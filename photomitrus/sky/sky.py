"""
Subtracts sky and divides out flat
"""
import os
from astropy.io import fits
import argparse
#%%
def subtract_sky_and_normalize(science_data_directory, output_data_dir, sky):
    if not os.path.isdir(output_data_dir):
        os.makedirs(output_data_dir)
    image_fnames = [os.path.join(science_data_directory, f) for f in os.listdir(science_data_directory) if f.endswith('fits.ramp') or f.endswith('ramp.new')]
    image_fnames.sort()
    sky = fits.getdata(sky)
    cropsky = sky[4:4092, 10:]  #remove if crop issue ever fixed
    for f in image_fnames:
        with fits.open(f) as hdul:
            image = hdul[0].data
            header = hdul[0].header
            cropimage = image[4:4092, 10:]  #remove if crop issue ever fixed
        reduced_image = (cropimage-cropsky)
        output_fname = os.path.basename(f)
        output_fname = output_fname.replace('.ramp.new', '.sky.flat.fits')
        output_fname = os.path.join(output_data_dir, output_fname)
        fits.HDUList(fits.PrimaryHDU(header=header, data=reduced_image)).writeto(output_fname, overwrite=True)
    print('Completed!')

#%%
def sky_flat_and_normalize(science_data_directory, output_data_dir, sky, flat):
    if not os.path.isdir(output_data_dir):
        os.makedirs(output_data_dir)
    image_fnames = [os.path.join(science_data_directory, f) for f in os.listdir(science_data_directory) if f.endswith('ramp.fits') or f.endswith('ramp.new')]
    image_fnames.sort()
    sky = fits.getdata(sky)
    flat = fits.getdata(flat)
    cropsky = sky[4:4092, 10:]  #remove if crop issue ever fixed
    for f in image_fnames:
        with fits.open(f) as hdul:
            image = hdul[0].data
            header = hdul[0].header
            cropimage = image[4:4092, 10:]  #remove if crop issue ever fixed
        sky_image = (cropimage-cropsky)
        reduced_image = sky_image / flat
        print(f + ' cropped, sky subtracted, and flat fielded!')
        output_fname = os.path.basename(f)
        output_fname = output_fname.replace('.ramp.new', '.sky.flat.fits')
        output_fname = os.path.join(output_data_dir, output_fname)
        print(output_fname)
        fits.HDUList(fits.PrimaryHDU(header=header, data=reduced_image)).writeto(output_fname, overwrite=True)
#%%

"""subtract_sky_and_normalize('/Users/orion/Desktop/PRIME/GRB/H_band/H_astrom_flat/', '/Users/orion/Desktop/PRIME/GRB/H_band/H_backs/',
                           '/Users/orion/Desktop/PRIME/GRB/H_band/H_sky/sky.2.5-4.fits')"""
#%%
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Crops and subtracts sky from files in dir, can also divide out flat')
    parser.add_argument('-flat', action='store_true', help='put optional arg if including flat fielding')
    parser.add_argument('-in_path', type=str, help='[str] Input imgs path (usually ramps w/ astrometry)')
    parser.add_argument('-out_path', type=str, help='[str] output sky sub image path')
    parser.add_argument('-sky_path', type=str, help='[str] input sky image path (for sky sub')
    parser.add_argument('-flat_path', type=str, help='[str] input flat path',default=None)
    args = parser.parse_args()
    if args.flat:
        sky_flat_and_normalize(args.in_path,args.out_path,args.sky_path,args.flat_path)
    else:
        subtract_sky_and_normalize(args.in_path,args.out_path,args.sky_path)