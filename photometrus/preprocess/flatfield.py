"""
Flat fields images
"""

from astropy.io import fits
import os
import argparse
import warnings
import shutil

warnings.filterwarnings("ignore", category=RuntimeWarning)
#%%


def flatfield(science_data, output_data_dir, flat):
    if not os.path.isdir(output_data_dir):
        os.makedirs(output_data_dir)
    flat = fits.getdata(flat)
    if type(science_data) is list:
        print('Flat fielding initial ramps...')
        for imagepath in science_data:
            # copy initial ramp into flat dir
            if imagepath.endswith('.fz'):
                filename = imagepath[-23:]
                flatnewname = filename.replace('.fits.ramp.fz', '.flat.fits')
                flatnewpath = os.path.join(output_data_dir, flatnewname)
                shutil.copyfile(imagepath, flatnewpath)
                os.system('funpack -F %s' % flatnewpath)
            else:
                filename = imagepath[-20:]
                if filename.endswith('.ramp.fits'):
                    flatnewname = filename.replace('.ramp.fits', '.flat.fits')
                else:
                    flatnewname = filename.replace('.fits.ramp', '.flat.fits')
                flatnewpath = os.path.join(output_data_dir, flatnewname)
                shutil.copyfile(imagepath, flatnewpath)

            # flat-field copied ramp
            with fits.open(flatnewpath) as hdul:
                image = hdul[0].data
                header = hdul[0].header
                cropimage = image[4:4092, 4:4092]  # remove if crop issue ever fixed
            ff_image = cropimage/flat
            # write interim temp file
            tmp_path = flatnewpath+'.tmp'
            fits.HDUList(fits.PrimaryHDU(header=header, data=ff_image)).writeto(tmp_path, overwrite=True)
            # replace temp w/ flat field path
            os.replace(tmp_path, flatnewpath)
        print('Flat fielding completed!')

    elif os.path.isdir(science_data):
        image_fnames = [os.path.join(science_data, f) for f in os.listdir(science_data) if f.endswith('ramp.fits') or f.endswith('ramp.new')]
        image_fnames.sort()
        print('Flat fielding imgs...')
        for f in image_fnames:
            with fits.open(f) as hdul:
                image = hdul[0].data
                header = hdul[0].header
                cropimage = image[4:4092, 4:4092]  # remove if crop issue ever fixed
            ff_image = cropimage/flat
            output_fname = os.path.basename(f)
            if f.endswith('.ramp.fits'):
                output_fname = output_fname.replace('.ramp.fits', '.flat.fits')
            else:
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


def flat_field_cmd(in_path, out_path, flat_path):
    flatfield(in_path, out_path, flat_path)


#%%
def main():
    parser = argparse.ArgumentParser(description='Flatfields image data')
    parser.add_argument('-in_path', type=str, help='[str] input img path, individual file, or file list (usually ramp images)')
    parser.add_argument('-out_path', type=str, help='[str] output flat fielded imgs path')
    parser.add_argument('-flat_path', type=str, help='[str] flat img path')
    args = parser.parse_args()

    flatfield(args.in_path, args.out_path, args.flat_path)


if __name__ == "__main__":
    main()