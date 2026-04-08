import os
import shutil
from pathlib import Path
import argparse
from tqdm import tqdm

import xarray as xr
import numpy as np
import dask.diagnostics
import dask.array as da 
from ome_zarr.io import parse_url
from ome_zarr.reader import Reader
from multiview_stitcher import spatial_image_utils as si_utils
from multiview_stitcher import (
    fusion,
    io,
    msi_utils,
    ngff_utils,
    registration
)
from pylibCZIrw import czi as pyczi # to get mosaic shape from czi file

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--dataPath", help="The path to your data")
parser.add_argument("--extension", help="The extension of the files to be processed. Options: '.zarr', '.czi'", default='.zarr')
parser.add_argument("--metadataSubstring", type=str, help="The substring to use for metadata extraction", default='AcquisitionBlock')
parser.add_argument("--tilePruningMethod", type=str, help="The method to use for tile pruning before registration. Options are: 'keep_axis_aligned', 'alternating_pattern', 'shortest_paths_overlap_weighted', 'otsu_threshold_on_overlap'", default='keep_axis_aligned')
parser.add_argument("--overlapTolerance", type=float, help="Extended overlap tolerance between tiles in percentage (between 0 and 1)", default=0)
parser.add_argument("--keepIntermediateFiles", type=bool, help="Whether to keep intermediate OME-Zarr files created for each tile. Options: True, False", default=False)

args = parser.parse_args()

# get Source path
if args.dataPath is None:
    print("Please provide a data path")
    exit(1)
basedir = Path(args.dataPath)

# get the image reader based on the file extension
# if args.extension == '.czi':
#     # from bioio import BioImage
#     # import bioio_czi
#     from pylibCZIrw import czi as pyczi
# else:
#     from tifffile import imread


def get_unique_names(array, substring='.'):
    """
    Get unique names from a string array
    """
    try:
        unique_names = [f[:f.index(substring)] for f in array]
        unique_names = list(set(unique_names))
        unique_names.sort()
    except:
        unique_names = []
    return unique_names

def get_filename_from_tile_and_channel(data_path, tile):
    """
    This convenience function returns the filename given the tile and channel.
    """
    return data_path / f'{tile}'

def get_tile_grid_position_from_tile_index(tile_index, num_cols):
    """
    Function to get the grid position of a tile based on its index.
    Based on original work from https://github.com/multiview-stitcher/multiview-stitcher
    """
    return {
        'z': 0,
        'y': tile_index // num_cols,
        'x': tile_index % num_cols if (tile_index // num_cols) % 2 == 0 else num_cols - 1 - (tile_index % num_cols),
    }

def get_mosaic_shape_from_parent_file(data_path, file_name, name_substring, file_extension):
    """
    Get Mosaic tile shape from parent czi file.
    This function reads the metadata from the parent file and returns the shape of the mosaic.
    """
    parent_path = data_path.parent
    parent_filelist = os.listdir(parent_path)
    parent_filelist_filtered = [f for f in parent_filelist if f.endswith('.czi') and file_name in f]
    if len(parent_filelist_filtered) == 0:
        raise FileNotFoundError("No parent file found for %s in %s" % (file_name, parent_path))
    
    parent_filelist_filtered.sort()
    parent_file_path = parent_path / parent_filelist_filtered[0]
    print('Reading metadata from %s' % parent_filelist_filtered[0])
    
    with pyczi.open_czi(str(parent_file_path)) as czidoc:
        md_dic = czidoc.metadata
    
    # getting the metadata block corresponding to this mosaic when there is multiple acquisition blocks
    md_block = md_dic['ImageDocument']['Metadata']['Experiment']['ExperimentBlocks']['AcquisitionBlock']
    if name_substring != 'AcquisitionBlock':
        n_rows = int(md_block['SubDimensionSetups']['RegionsSetup']['SampleHolder']['TileRegions']['TileRegion']['Rows'])
        n_cols = int(md_block['SubDimensionSetups']['RegionsSetup']['SampleHolder']['TileRegions']['TileRegion']['Columns'])
    else:
        idx_start = parent_filelist_filtered[0].index(name_substring)
        offset = len(name_substring)
        if file_extension == '.czi':
            idx_end = parent_filelist_filtered[0].index('-', idx_start + offset, len(parent_filelist_filtered[0]))
        else:
            idx_1 = parent_filelist_filtered[0].index('_', idx_start + offset, len(parent_filelist_filtered[0]))
            idx_2 = parent_filelist_filtered[0].index('-', idx_start + offset, len(parent_filelist_filtered[0]))
            idx_end = min(idx_1, idx_2) if (idx_1 >=0 and idx_2 >=0) else max(idx_1, idx_2) # get the valid index when found
        block_index = int(parent_filelist_filtered[0][idx_start + offset : idx_end]) - 1 # 0-based index
        md_tileregions = md_block[block_index]['SubDimensionSetups']['RegionsSetup']['SampleHolder']['TileRegions']['TileRegion']
        if type(md_tileregions) is list: # when there are multiple tile regions it becomes a list, else it is a dict with TileRegion properties
            i = 0
            tile_used = False
            while i < len(md_tileregions) and tile_used is False:
                if md_tileregions[i]['IsUsedForAcquisition'] == 'true':
                    tile_used = True
                else:
                    i += 1
            n_rows = int(md_tileregions[i]['Rows'])
            n_cols = int(md_tileregions[i]['Columns'])
        else:
            n_rows = int(md_tileregions['Rows'])
            n_cols = int(md_tileregions['Columns'])
    return n_rows, n_cols


def tile_registration(data_array, overlap_tolerance, tile_pruning_method):
    """
    Wrapping function for tile stitching and registration. 
    """
    # Do channel alignment if needed
    perform_channel_alignment = False

    # Channel alignment is performed using a single tile.
    # Choose here which tile (index) to use.
    channel_alignment_tile_index = 0

    curr_transform_key = 'affine_metadata'
    if perform_channel_alignment:
        curr_transform_key = 'affine_metadata_ch_reg'

        channels = data_array[channel_alignment_tile_index]['scale0/image'].coords['c']

        # select chosen tiles for registration
        msims_ch_reg = [msi_utils.multiscale_sel_coords(data_array[5], {'c': ch})
                        for ch in channels]

        with dask.diagnostics.ProgressBar():
            params_c = registration.register(
                msims_ch_reg,
                registration_binning={'z': 2, 'y': 4, 'x': 4},
                reg_channel_index=0,
                transform_key='affine_metadata',
                pre_registration_pruning_method=None,
            )

        # assign channel coordinates to obtained parameters
        params_c = xr.concat(params_c, dim='c').assign_coords({'c': channels})

        # set obtained parameters for all tiles
        for msim in data_array:
            msi_utils.set_affine_transform(
                msim, params_c, transform_key=curr_transform_key, base_transform_key='affine_metadata')
    
    # do registration
    print('Performing registration...')
    with dask.diagnostics.ProgressBar():
        params = registration.register(
            data_array,
            registration_binning={'z': 1, 'y': 2, 'x': 2},
            reg_channel_index=0,
            transform_key=curr_transform_key,
            overlap_tolerance=overlap_tolerance,
            new_transform_key='affine_registered',
            pre_registration_pruning_method=tile_pruning_method,
        )
    
    # print obtained registration parameters
    for imsim, msim in enumerate(data_array):
        affine = np.array(msi_utils.get_transform_from_msim(msim, transform_key='affine_registered')[0])
        print('tile index %s \n %s' % (imsim, affine))
    
    return params, affine


def main(datapath='.', extension='.czi', metadata_substring='AcquisitionBlock', tile_pruning_method='keep_axis_aligned', overlap_tolerance=0, keep_intermediate_files=False):
    print('Processing folder: %s' % datapath)
    filelist = os.listdir(datapath)

    if extension.find('.') == -1:
        extension = '.' + extension
    
    filelist = [f for f in filelist if f.endswith(extension)]
    filelist.sort()
    print('Nr of czi files in dir: %i' % len(filelist))

    savedir = Path(str(datapath) + '/stitched_tiles_2d/')
    savedir.mkdir(parents=True, exist_ok=True)
    print('Saving output to: %s' % savedir)

    original_filenames = get_unique_names(filelist, substring='_tile')
    print("Nb of unique file names: %i" % len(original_filenames))

    for original_name in original_filenames:
        filelist_filtered = []
        for name in filelist:
            if name.find(original_name) >= 0:
                filelist_filtered.append(name)

        n_tiles = int(len(filelist_filtered))
        n_rows, n_columns = get_mosaic_shape_from_parent_file(datapath, original_name, metadata_substring, extension)
        print('Number of tiles: %i' % n_tiles)
        print('Mosaic shape: %i rows, %i columns' % (n_rows, n_columns))

        tile_file_indexes = [] # get tile file indexes from filelist
        for i in range(n_tiles):
            for file in filelist:
                if file.find(original_name) >= 0 and file.endswith('_tile' + str(i + 1).zfill(2) + extension):
                    tile_file_indexes.append(filelist.index(file))
        tile_file_indexes.sort()

        filelist_tiles = [filelist[i] for i in tile_file_indexes]
        print('\n '.join([f for f in filelist_tiles]))
        print('Tile grid indices:')
        print("\n".join([f"Tile {itile}: " + str(get_tile_grid_position_from_tile_index(itile, n_columns)) for itile, tile in enumerate(tile_file_indexes)]))

        # Getting image data voxel scales
        if extension == '.czi':
            # For CZI files, we need to read the first file to get the metadata
            file_path = str(datapath / filelist_tiles[0])
            # img = BioImage(
            #     file_path, 
            #     reader=bioio_czi.Reader, 
            #     reconstruct_mosaic=False,
            #     include_subblock_metadata=True,
            #     use_aicspylibczi=True,
            # )
            # scale = {'z': img.scale.Z, 'y': img.scale.Y, 'x': img.scale.X}

            with pyczi.open_czi(file_path) as czidoc:
                md_dic = czidoc.metadata
                tbd = czidoc.total_bounding_box
                pixelsize_x = float(md_dic['ImageDocument']['Metadata']['Scaling']['Items']['Distance'][0]['Value'])
                pixelsize_y = float(md_dic['ImageDocument']['Metadata']['Scaling']['Items']['Distance'][1]['Value'])
                pixelsize_z = float(md_dic['ImageDocument']['Metadata']['Scaling']['Items']['Distance'][2]['Value'])
            
            pixelsize_x = pixelsize_x / 10**-6  # convert scale from microns to meters
            pixelsize_y = pixelsize_y / 10**-6  # convert scale from microns to meters
            pixelsize_z = pixelsize_z / 10**-6  # convert scale from microns to meters

            scale = {'z': pixelsize_z, 'y': pixelsize_y, 'x': pixelsize_x}
        else:
            file_path = str(datapath / filelist_tiles[0])
            store = parse_url(file_path, mode="r")
            reader = Reader(store)
            nodes = list(reader())
            image_node = nodes[0]  # Get the first image

            # Get the pixel sizes (scales) of the first, raw, pyramid image
            scales = image_node.metadata['coordinateTransformations'][0][0]['scale']
            scale = {'z': scales[-3], 'y': scales[-2], 'x': scales[-1]}
        
        print('Voxel scales: %s' % scale)

        overlap = {
            'x': 0.1,
            'y': 0.1,
            # 'z': 0.1
        }
        if extension == '.czi':
            tile_shape = {
                'z': tbd['Z'][1],
                'y': tbd['Y'][1],
                'x': tbd['X'][1]
            }
        else:
            tile_shape = {
                'z': image_node.data[0].shape[-3], # to get the Z dimension
                'y': image_node.data[0].shape[-2], # to get the Y dimension
                'x': image_node.data[0].shape[-1] # to get the X dimension
            }
        print('Tile shape: %s' % tile_shape)

        translations = []
        for itile, tile in enumerate(tile_file_indexes):
            tile_grid_position = get_tile_grid_position_from_tile_index(itile, n_columns)
            translations.append(
                {
                    dim: tile_grid_position[dim] * (1 - (overlap[dim] if dim in overlap else 1)) * tile_shape[dim] * scale[dim]
                    for dim in scale
                }
            )

        print("Tile positions:")
        print("\n".join([f"Tile {itile}: " + str(t) for itile, t in enumerate(translations)]))

        # Read input tiles, convert to OME-Zarr files, then delete temporary files
        overwrite = True

        filelist_savenames = []
        if extension == '.czi':
            filelist_savenames = [f[:f.index(extension)].replace(' ', '_') + '.zarr' for f in filelist_tiles]
            print('Saving OME-Zarr files with names:')
            print('\n'.join([i for i in filelist_savenames]))

        msims = []
        zarr_paths = []
        for itile, tile in enumerate(tqdm(filelist_tiles)):

            # set save path for OME-Zarr files
            if extension == '.czi':
                zarr_path = os.path.join(os.path.dirname(get_filename_from_tile_and_channel(datapath, tile)), filelist_savenames[itile])
            else:
                zarr_path = os.path.join(os.path.dirname(get_filename_from_tile_and_channel(datapath, tile)), filelist_tiles[itile])

            # read tile image
            if os.path.exists(zarr_path):
                im_data = da.from_zarr(os.path.join(zarr_path, '0'))[0] # drop t axis automatically added
            else:
                file_path = str(datapath / tile)
                # img = BioImage(
                #     file_path, 
                #     reader=bioio_czi.Reader, 
                #     reconstruct_mosaic=False,
                #     include_subblock_metadata=True,
                #     use_aicspylibczi=True,
                # )
                # get data dimensions without T axis from metadata
                # im_data = img.get_image_data(img.dims.order[img.dims.order.index('T')+1:])

                with pyczi.open_czi(file_path) as cziimg:
                        tbd = cziimg.total_bounding_box
                        im_data = np.zeros((tbd['C'][1], tbd['Z'][1], tbd['Y'][1], tbd['X'][1]))

                        for t in range(tbd['T'][1]):
                            for c in range(tbd['C'][1]):
                                for z in range(tbd['Z'][1]):
                                    temp = cziimg.read(
                                        plane = {'C': c, "T": t, "Z": z},
                                        scene = 0,
                                    )
                                    im_data[c, z] = temp.squeeze()

            sim = si_utils.get_sim_from_array(
                im_data,
                dims=["c", "z", "y", "x"],
                scale=scale,
                translation=translations[itile],
                transform_key=io.METADATA_TRANSFORM_KEY,
                )
            
            if extension == '.czi':
                # write to OME-Zarr
                ngff_utils.write_sim_to_ome_zarr(sim, zarr_path, overwrite=overwrite)
                # replace sim with the sim read from the written OME-Zarr
                sim = ngff_utils.read_sim_from_ome_zarr(zarr_path)

            msim = msi_utils.get_msim_from_sim(sim)
            zarr_paths.append(zarr_path)

            msims.append(msim)

        try:
            params, affine = tile_registration(msims, overlap_tolerance, tile_pruning_method)
        except:
            print('Tile registration failed. Skipping this tile set.')
            print('====================')
            continue
        
        save_name = filelist_tiles[0][:filelist_tiles[0].index('_tile')] + '.zarr'
        
        print('Save name: %s' % save_name)
        output_filename = os.path.join(savedir, save_name)

        print('Fusing views...')
        fused = fusion.fuse(
            [msi_utils.get_sim_from_msim(msim) for msim in msims],
            transform_key='affine_registered',
            output_chunksize=256,
            )

        print('Fusing views and saving output to %s...', output_filename)
        with dask.diagnostics.ProgressBar():
            fused = ngff_utils.write_sim_to_ome_zarr(
                fused, output_filename, overwrite=True
            )
        
        print('Removing temporary files...')
        if extension == '.czi' and not keep_intermediate_files:
            for itile, tile in enumerate(tqdm(filelist_tiles)):
                zarr_path = os.path.join(os.path.dirname(get_filename_from_tile_and_channel(datapath, tile)), filelist_savenames[itile])
                if os.path.exists(zarr_path):
                    shutil.rmtree(zarr_path)

        print('====================')
    print('Done!')


if __name__ == '__main__':
    main(datapath=basedir, extension=args.extension, metadata_substring=args.metadataSubstring, tile_pruning_method=args.tilePruningMethod, overlap_tolerance=args.overlapTolerance, keep_intermediate_files=args.keepIntermediateFiles)
