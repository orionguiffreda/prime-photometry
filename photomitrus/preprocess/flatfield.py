"""
Flat fields images
"""

from astropy.io import fits
import os
import argparse
#%%

def flatfield(science_data_directory, output_data_dir, flat):
    if not os.path.isdir(output_data_dir):
        os.makedirs(output_data_dir)
    image_fnames = [os.path.join(science_data_directory, f) for f in os.listdir(science_data_directory) if f.endswith('ramp.fits') or f.endswith('ramp.new')]
    image_fnames.sort()
    flat = fits.getdata(flat)
    print('Flat fielding imgs...')
    for f in image_fnames:
        with fits.open(f) as hdul:
            image = hdul[0].data
            header = hdul[0].header
            cropimage = image[4:4092, 4:4092]  #remove if crop issue ever fixed
        ff_image = cropimage/flat
        output_fname = os.path.basename(f)
        output_fname = output_fname.replace('.ramp.new', '.flat.fits')
        output_fname = os.path.join(output_data_dir, output_fname)
        fits.HDUList(fits.PrimaryHDU(header=header, data=ff_image)).writeto(output_fname, overwrite=True)
    print('Flat fielding completed!')

#%%
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Flatfields image data')
    parser.add_argument('-in_path', type=str, help='[str] input img path (usually astrom ramp images)')
    parser.add_argument('-out_path', type=str, help='[str] output flat fielded imgs path')
    parser.add_argument('-flat_path', type=str, help='[str] flat img path')
    args = parser.parse_args()

    flatfield(args.in_path,args.out_path,args.flat_path)