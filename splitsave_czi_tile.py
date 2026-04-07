import os
from pathlib import Path
import argparse
from tqdm import tqdm
import numpy as np

from bioio import BioImage
import bioio_czi
from pylibCZIrw import czi as pyczi

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--dataPath", help="The path to your data")
parser.add_argument("--extension", help="The extension of the files to be processed", default='.czi')

args = parser.parse_args()
print(args.dataPath)

if args.dataPath is None:
    print("Please provide a data path")
    exit(1)

basedir = Path(args.dataPath)

def main(datapath='.', extension='.czi'):
    filelist = os.listdir(datapath)

    filelist = [f for f in filelist if f.find(extension) > 0]
    filelist.sort()
    print('Nr of czi files in dir:', len(filelist))

    savedir = Path(str(basedir) + '/split_czi/')
    savedir.mkdir(parents=True, exist_ok=True)
    print('Saving output to:', savedir)

    for file in tqdm(filelist, desc='Processing files'):
        filepath = datapath / file
        filename_noext = file[:file.index(extension)]
        filename_noext = filename_noext.replace(' ', '_')
        print('Processing file:', filepath)

        img = BioImage(
            filepath,
            reader=bioio_czi.Reader,
            reconstruct_mosaic=False,
            include_subblock_metadata=True,
            use_aicspylibczi=True,
        )
        print('Image dimensions:', img.dims)

        if img.dims.order.find('M') < 0:
            print('Image is not mosaic, skipping')
            continue

        for scene in tqdm(img.scenes, desc='Processing scenes'):
            img.set_scene(scene)
            print('Current scene:', scene)

            n_tiles = img.dims['M'][0]

            scale_x = img.scale.X * 10**-6 # convert scale from BioIO to default CZI unit, which is meters
            scale_y = img.scale.Y * 10**-6
            scale_z = img.scale.Z * 10**-6

            for tile in tqdm(range(n_tiles), desc='Processing tiles'):
                img_data = img.get_image_dask_data(img.dims.order[1:], M=tile) #TODO better strategy to exclude dimension to split
                img_data_tile = img_data.compute()

                ch_names = {}
                for i in range(len(img.channel_names)):
                    ch_names[i] = img.channel_names[i]

                tile_save_path = str(savedir) + '/' + filename_noext + '_tile' + str(tile+1).zfill(2) + '.czi'
                print('Tile save path:', tile_save_path)

                if img_data_tile.shape[2] != img.dims['Z'][0]:
                    print('Warning: Z dimension of tile data does not match original image, filling with black slices')
                    # create new array with same shape as original image, but with tile data in the correct position and black slices for missing Z planes
                    new_shape = list(img_data_tile.shape)
                    new_shape[2] = img.dims['Z'][0]
                    img_data_tile_filled = np.zeros(tuple(new_shape), dtype=img_data_tile.dtype)
                    img_data_tile_filled[:, :, :img_data_tile.shape[2], :, :] = img_data_tile
                    img_data_tile = img_data_tile_filled
                
                with pyczi.create_czi(tile_save_path, exist_ok=True) as czidoc_w:
                    for t in range(img.dims['T'][0]):
                        for c in range(img.dims['C'][0]):
                            for z in range(img.dims['Z'][0]):
                                temp_image = img_data_tile[t][c][z]
                                czidoc_w.write(
                                    data=temp_image,
                                    plane={
                                        'T': t,
                                        'C': c,
                                        'Z': z,
                                    },
                                    compression_options = "zstd0:ExplicitLevel=2"
                                )
                    
                    czidoc_w.write_metadata(
                        filename_noext + '_tile' + str(tile+1).zfill(2),
                        channel_names=ch_names,
                        scale_x=scale_x,
                        scale_y=scale_y,
                        scale_z=scale_z,
                    )
                

if __name__ == '__main__':
    main(basedir, args.extension)