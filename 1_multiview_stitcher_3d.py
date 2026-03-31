import os
import shutil
from pathlib import Path
import argparse
from tqdm import tqdm

import xarray as xr
import numpy as np
import dask.diagnostics
import dask.array as da 
from multiview_stitcher import spatial_image_utils as si_utils
from multiview_stitcher import (
    fusion,
    io,
    msi_utils,
    ngff_utils,
    registration
)

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--dataPath", help="The path to your data")
parser.add_argument("--extension", help="The extension of the files to be processed", default='.czi')
parser.add_argument("--overlapTolerance", type=float, help="Extended overlap tolerance between tiles in percentage (between 0 and 1)", default=0)

args = parser.parse_args()

# get Source path
if args.dataPath is None:
    print("Please provide a data path")
    exit(1)
basedir = Path(args.dataPath)

# get the image reader based on the file extension
if args.extension == '.czi':
    # from bioio import BioImage
    # import bioio_czi
    from pylibCZIrw import czi as pyczi
elif args.extension == '.zarr':
    from ome_zarr.io import parse_url
    from ome_zarr.reader import Reader
else:
    from tifffile import imread


def get_unique_names(array, substring='.'):
    """
    Get unique names from a string array
    Parameters
    ----------
    array : list of strings
        List of strings to be processed.
    substring : string
        Substring to identify unique names. Will use substring position to cut the string.
    Returns
    ----------
    unique_names : list of strings
        List of unique names.
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

def get_tile_grid_position_from_tile_index(tile_index, z_planes):
    """
    Function to get the grid position of a tile based on its index.
    Based on original work from https://github.com/multiview-stitcher/multiview-stitcher
    """
    return {
        'z': tile_index % z_planes,
        'y': 0,
        'x': 0
    }

def tile_registration(data_array, overlap_tolerance):
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
            pre_registration_pruning_method="keep_axis_aligned",
        )
    
    # print obtained registration parameters
    for imsim, msim in enumerate(data_array):
        affine = np.array(msi_utils.get_transform_from_msim(msim, transform_key='affine_registered')[0])
        print('tile index %s \n %s' % (imsim, affine))
    
    return params, affine


def main(datapath='.', extension='.czi', overlap_tolerance=0):
    print('Processing folder: %s' % datapath)
    filelist = os.listdir(datapath)

    filelist = [f for f in filelist if f.endswith(extension)]
    filelist.sort()
    print('Nr of %s files in dir: %i' % (extension, len(filelist)))

    savedir = Path(str(datapath) + '/stitched_tile_3d/')
    savedir.mkdir(parents=True, exist_ok=True)
    print('Saving output to: %s' % savedir)

    original_filenames = get_unique_names(filelist, substring='_sub')
    if len(original_filenames) == 0:
        original_filenames = get_unique_names(filelist, substring='_Sub')
    print("Nb of unique file names: %i" % len(original_filenames))

    for original_name in original_filenames:
        filelist_filtered = []
        for name in filelist:
            if name.find(original_name) >= 0:
                filelist_filtered.append(name)

        filelist_substacks = get_unique_names(filelist_filtered, substring='-Scene')
        n_positions = int(len(filelist_filtered) / len(filelist_substacks))
        n_substacks = len(filelist_substacks)

        for i in range(n_positions):
            substack_file_indexes = []
            if Path(datapath).stem == 'split_czi':
                for file in filelist:
                    if file.find(original_name) >= 0 and file.endswith('_tile' + str(i + 1).zfill(2) + extension):
                        substack_file_indexes.append(filelist.index(file))
            else:
                for file in filelist:
                    if file.find(original_name) >= 0 and file.endswith(extension):
                        substack_file_indexes.append(filelist.index(file))
            substack_file_indexes.sort()

            filelist_tiles = [filelist[i] for i in substack_file_indexes]
            print('\n '.join([x for x in filelist_tiles]))
            print('Tile grid indices:')
            print("\n".join([f"Tile {itile}: " + str(get_tile_grid_position_from_tile_index(itile, n_substacks)) for itile, tile in enumerate(substack_file_indexes)]))

            # Getting image data voxel scales
            file_path = str(datapath / filelist_tiles[0])

            if extension == '.czi':
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
                store = parse_url(file_path, mode="r")
                reader = Reader(store)
                nodes = list(reader())
                image_node = nodes[0]  # Get the first image

                # Get the pixel sizes (scales) of the first, raw, pyramid image
                scales = image_node.metadata['coordinateTransformations'][0][0]['scale']
                scale = {'z': scales[-3], 'y': scales[-2], 'x': scales[-1]}
            
            print('Voxel scales: %s' % scale)

            overlap = {
                # 'x': 0.1,
                # 'y': 0.1,
                'z': 0.1
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
            for itile, tile in enumerate(substack_file_indexes):
                tile_grid_position = get_tile_grid_position_from_tile_index(itile, n_substacks)
                translations.append(
                    {
                        dim: tile_grid_position[dim] * (1 - (overlap[dim] if dim in overlap else 1)) * tile_shape[dim] * scale[dim]
                        for dim in scale
                    }
                )

            print("Tile positions:")
            print("\n".join([f"Tile {itile}: " + str(t) for itile, t in enumerate(translations)]))

            # Read input tiles, convert to OME-Zarr files, then delete temporary files
            # remove spaces from filename
            filelist_savenames = [f[:f.index(extension)].replace(' ', '_') + '.zarr' for f in filelist_tiles]
            print('Saving OME-Zarr files with names:')
            print('\n'.join([i for i in filelist_savenames]))

            msims = []
            zarr_paths = []
            overwrite = True
            for itile, tile in enumerate(tqdm(filelist_tiles)):

                # set save path for OME-Zarr files
                zarr_path = os.path.join(os.path.dirname(get_filename_from_tile_and_channel(datapath, tile)), filelist_savenames[itile])
                print('OME-Zarr path: %s' % zarr_path)

                # read tile image
                # if os.path.exists(zarr_path) and not overwrite:
                if extension == '.zarr' and os.path.exists(zarr_path):
                    im_data = da.from_zarr(os.path.join(zarr_path, '0'))[0] # drop t axis automatically added
                else:
                    file_path = str(datapath / tile)
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

                if extension != '.zarr':
                    # write to OME-Zarr
                    ngff_utils.write_sim_to_ome_zarr(sim, zarr_path, overwrite=overwrite)
                    # replace sim with the sim read from the written OME-Zarr
                    sim = ngff_utils.read_sim_from_ome_zarr(zarr_path)

                msim = msi_utils.get_msim_from_sim(sim)
                zarr_paths.append(zarr_path)

                msims.append(msim)

            params, affine = tile_registration(msims, overlap_tolerance=overlap_tolerance)
            
            try:
                save_name = filelist_savenames[0][:filelist_savenames[0].index('_sub')] + '_tile'+ str(i + 1).zfill(2) + '.zarr'
            except:
                save_name = filelist_savenames[0][:filelist_savenames[0].index('_Sub')] + '_tile'+ str(i + 1).zfill(2) + '.zarr'
            
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
            
            if extension == '.czi':
                print('Removing temporary files...')
                for itile, tile in enumerate(tqdm(filelist_tiles)):
                    zarr_path = os.path.join(os.path.dirname(get_filename_from_tile_and_channel(datapath, tile)), filelist_savenames[itile])
                    if os.path.exists(zarr_path):
                        shutil.rmtree(zarr_path)
            
            print('====================')
    print('Done!')


if __name__ == '__main__':
    main(datapath=basedir, extension=args.extension, overlap_tolerance=args.overlapTolerance)