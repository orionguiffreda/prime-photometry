"""
Subtracts sky
"""
import os
from astropy.io import fits
import subprocess
import argparse
from settings import gen_config_file_name
import numpy as np
#%%
def subtract_sky_and_normalize(science_data_directory, output_data_dir, sky):
    if not os.path.isdir(output_data_dir):
        os.makedirs(output_data_dir)
    image_fnames = [os.path.join(science_data_directory, f) for f in os.listdir(science_data_directory) if f.endswith('ramp.new')]
    image_fnames.sort()
    sky = fits.getdata(sky)
    cropsky = sky [4:4092, 4:4092]  #remove if crop issue ever fixed
    for f in image_fnames:
        with fits.open(f) as hdul:
            image = hdul[0].data
            header = hdul[0].header
            cropimage = image [4:4092, 4:4092]  #remove if crop issue ever fixed
            CRPIX1 = (header['CRPIX1'])         #changing ref pixels to work w/ cropped imgs
            CRPIX2 = (header['CRPIX2'])
            header.set('CRPIX1', value=CRPIX1 - 4)
            header.set('CRPIX2', value=CRPIX2 - 4)
        reduced_image = (cropimage-cropsky*np.nanmedian(cropimage))
        output_fname = os.path.basename(f)
        output_fname = output_fname.replace('.ramp.new', '.sky.flat.fits')
        output_fname = os.path.join(output_data_dir, output_fname)
        fits.HDUList(fits.PrimaryHDU(header=header, data=reduced_image)).writeto(output_fname, overwrite=True)
    print('Completed!')

#%%
def sky_flat_and_normalize(science_data_directory, output_data_dir, sky):
    if not os.path.isdir(output_data_dir):
        os.makedirs(output_data_dir)
    image_fnames = [os.path.join(science_data_directory, f) for f in os.listdir(science_data_directory) if f.endswith('flat.fits')]
    image_fnames.sort()
    sky = fits.getdata(sky)
    cropsky = sky #[4:4092, 4:4092]  #remove if crop issue ever fixed
    for f in image_fnames:
        with fits.open(f) as hdul:
            image = hdul[0].data
            header = hdul[0].header
            cropimage = image #[4:4092, 4:4092]  #remove if crop issue ever fixed
            CRPIX1 = (header['CRPIX1'])  # changing ref pixels to work w/ cropped imgs
            CRPIX2 = (header['CRPIX2'])
            header.set('CRPIX1', value=CRPIX1 - 4)
            header.set('CRPIX2', value=CRPIX2 - 4)
        reduced_image = (cropimage-cropsky*np.nanmedian(cropimage))
        output_fname = os.path.basename(f)
        output_fname = output_fname.replace('.flat.fits', '.sky.flat.fits')
        output_fname = os.path.join(output_data_dir, output_fname)
        fits.HDUList(fits.PrimaryHDU(header=header, data=reduced_image)).writeto(output_fname, overwrite=True)
    print('Sky sub on FF imgs completed!')

#%%
def sexback(imgdir,outdir):
    os.chdir(str(imgdir))
    sx = gen_config_file_name('sexback.config')
    ap = gen_config_file_name('astrom.param')
    for f in sorted(os.listdir(str(imgdir))):
        if f.endswith('flat.fits'):
            output_fname = f.replace('.flat.fits', '.sky.flat.fits')
            pre = os.path.splitext(output_fname)[0]
            ext = os.path.splitext(output_fname)[1]
            com = ["sex ", imgdir + f, ' -c '+sx, " -CATALOG_NAME " + outdir + pre + '.cat', ' -PARAMETERS_NAME '+ap,
                   ' -CHECKIMAGE_NAME  '+ outdir + output_fname]
            s0 = ''
            com = s0.join(com)
            out = subprocess.Popen([com], shell=True,  stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            out.wait()
            print(output_fname + ' back subbed!')
    print('sxtrctr back sub complete!')
#%%

"""subtract_sky_and_normalize('/Users/orion/Desktop/PRIME/GRB/H_band/H_astrom_flat/', '/Users/orion/Desktop/PRIME/GRB/H_band/H_backs/',
                           '/Users/orion/Desktop/PRIME/GRB/H_band/H_sky/sky.2.5-4.fits')"""
#%%
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Crops and subtracts sky from files in dir, can also divide out flat')
    parser.add_argument('-flat', action='store_true', help='put optional arg if included flat fielding previously')
    parser.add_argument('-sex', action='store_true', help='optional arg to use sxtrctr background sub instead, '
                                                          'outputs .cats and sky subbed imgs')
    parser.add_argument('-in_path', type=str, help='[str] Input imgs path (usually ramps w/ astrometry)')
    parser.add_argument('-out_path', type=str, help='[str] output sky sub image path')
    parser.add_argument('-sky_path', type=str, help='[str] input sky image path (for sky sub')
    args = parser.parse_args()
    if args.flat:
        sky_flat_and_normalize(args.in_path,args.out_path,args.sky_path)
    elif args.sex:
        sexback(args.in_path, args.out_path)
    else:
        subtract_sky_and_normalize(args.in_path,args.out_path,args.sky_path)