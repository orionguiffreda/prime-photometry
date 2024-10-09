"""
Flat fields images
"""

from astropy.io import fits
import os
import argparse
import warnings

warnings.filterwarnings("ignore", category=RuntimeWarning)
#%%

def flatfield(science_data, output_data_dir, flat):
    if not os.path.isdir(output_data_dir):
        os.makedirs(output_data_dir)
    flat = fits.getdata(flat)
    if os.path.isdir(science_data):
        image_fnames = [os.path.join(science_data, f) for f in os.listdir(science_data) if f.endswith('ramp.fits') or f.endswith('ramp.new')]
        image_fnames.sort()
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
    elif os.path.isfile(science_data):
        hdu = fits.open(science_data)
        image = hdu[0].data
        header = hdu[0].header
        cropimage = image[4:4092, 4:4092]
        ff_image = cropimage / flat
        output_fname = os.path.basename(science_data)
        output_fname = output_fname.replace('.ramp.new', '.flat.fits')
        output_fname = os.path.join(output_data_dir, output_fname)
        fits.HDUList(fits.PrimaryHDU(header=header, data=ff_image)).writeto(output_fname, overwrite=True)


#%%
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Flatfields image data')
    parser.add_argument('-in_path', type=str, help='[str] input img path or individual file (usually astrom ramp images)')
    parser.add_argument('-out_path', type=str, help='[str] output flat fielded imgs path')
    parser.add_argument('-flat_path', type=str, help='[str] flat img path')
    args = parser.parse_args()

    flatfield(args.in_path,args.out_path,args.flat_path)