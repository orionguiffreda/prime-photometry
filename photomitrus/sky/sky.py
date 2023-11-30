"""
Subtracts sky
"""
import os
from astropy.io import fits
import argparse
#%%
def subtract_sky_and_normalize(science_data_directory, output_data_dir, sky):
    if not os.path.isdir(output_data_dir):
        os.makedirs(output_data_dir)
    image_fnames = [os.path.join(science_data_directory, f) for f in os.listdir(science_data_directory) if f.endswith('flat.fits')]
    image_fnames.sort()
    sky = fits.getdata(sky)
    for f in image_fnames:
        print(f)
        with fits.open(f) as hdul:
            image = hdul[0].data
            header = hdul[0].header
        reduced_image = (image-sky)
        output_fname = os.path.basename(f)
        output_fname = output_fname.replace('.flat.fits', '.sky.flat.fits')
        output_fname = os.path.join(output_data_dir, output_fname)
        print(output_fname)
        fits.HDUList(fits.PrimaryHDU(header=header, data=reduced_image)).writeto(output_fname, overwrite=True)

#%%

"""subtract_sky_and_normalize('/Users/orion/Desktop/PRIME/GRB/H_band/H_astrom_flat/', '/Users/orion/Desktop/PRIME/GRB/H_band/H_backs/',
                           '/Users/orion/Desktop/PRIME/GRB/H_band/H_sky/sky.2.5-4.fits')"""
#%%
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Subtracts sky from files in dir')
    parser.add_argument('paths', nargs='*', type=str, metavar='p', help='Put input path first, output path, then sky path')
    args = parser.parse_args()
    subtract_sky_and_normalize(args.paths[0],args.paths[1],args.paths[2])