"""
Subtracts sky
"""
import os
from astropy.io import fits
import subprocess
import argparse
import numpy as np
import threading
import sys

from photomitrus.settings import gen_config_file_name
from photomitrus.sky.gen_sky import checkplot

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
    cropsky = fits.getdata(sky)
    hdr_sky = fits.getheader(sky)
    try:
        airmass_sky = hdr_sky['AIRMASS']
    except KeyError:
        print('No Airmass value in sky header')
        airmass_sky = 0
    for f in image_fnames:
        with fits.open(f) as hdul:
            image = hdul[0].data
            header = hdul[0].header
            cropimage = image
            # CRPIX1 = (header['CRPIX1'])  # changing ref pixels to work w/ cropped imgs (NOT NECESSARY IF NO ASTROM)
            # CRPIX2 = (header['CRPIX2'])
            # header.set('CRPIX1', value=CRPIX1 - 4)
            # header.set('CRPIX2', value=CRPIX2 - 4)
            try:
                airmass_sci = header['AIRMASS']
            except KeyError:
                print('No Airmass value in sci header')
                airmass_sci = 0
        reduced_image = (cropimage-cropsky*np.nanmedian(cropimage))
        # if not airmass_sky or not airmass_sci:
        #     print('Airmasses not found in either sci or sky, defaulting to normal scaling')
        # else:
        #     airmass_ratio = airmass_sci / airmass_sky
        #     print('Airmass ratio =',airmass_ratio)
        #     reduced_image = (cropimage - (cropsky * np.nanmedian(cropimage) * airmass_ratio))
        output_fname = os.path.basename(f)
        output_fname = output_fname.replace('.flat.fits', '.sky.flat.fits')
        output_fname = os.path.join(output_data_dir, output_fname)
        fits.HDUList(fits.PrimaryHDU(header=header, data=reduced_image)).writeto(output_fname, overwrite=True)
    print('Sky sub on FF imgs completed!')

#%%


def sexback(imgdir,outdir):
    print('Using sextractor background subtraction...')
    os.chdir(str(imgdir))
    sx = gen_config_file_name('bulge_new.config')
    ap = gen_config_file_name('tempsource.param')

    def sxbackcmd(imgpath, sx, catpath, ap, outpath, backpath, first=False):
        if first:
            command = ('sex %s -c %s -CATALOG_NAME %s -PARAMETERS_NAME %s -CHECKIMAGE_TYPE -BACKGROUND,BACKGROUND'
                       ' -CHECKIMAGE_NAME %s,%s'
                       % (imgpath, sx, catpath, ap, outpath, backpath))
            rval = subprocess.run(command.split(), check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        else:
            command = ('sex %s -c %s -CATALOG_TYPE NONE -PARAMETERS_NAME %s -CHECKIMAGE_TYPE -BACKGROUND,BACKGROUND'
                       ' -CHECKIMAGE_NAME %s,%s'
                       % (imgpath, sx, ap, outpath, backpath))
            rval = subprocess.run(command.split(), check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    threads = []
    files = [f for f in sorted(os.listdir(str(imgdir))) if f.endswith('flat.fits')]
    if files:
        output_fname = files[0].replace('.flat.fits', '.sky.flat.fits')
        pre = os.path.splitext(output_fname)[0]
        imgpath = os.path.join(imgdir, files[0])
        catpath = os.path.join(outdir, pre+'.cat')
        outpath = os.path.join(outdir, output_fname)
        backname = output_fname.replace('.sky.flat.fits', '.sky.flat.back.fits')
        backpath = os.path.join(outdir, backname)
        sxbackcmd(imgpath, sx, catpath, ap, outpath, backpath, first=True)
        checkplot(output_directory=outdir, save_name=backname)
        for f in files[1:]:
            output_fname = f.replace('.flat.fits', '.sky.flat.fits')
            pre = os.path.splitext(output_fname)[0]
            imgpath = os.path.join(imgdir, f)
            catpath = os.path.join(outdir, pre+'.cat')
            outpath = os.path.join(outdir, output_fname)
            backname = output_fname.replace('.sky.flat.fits', '.sky.flat.back.fits')
            backpath = os.path.join(outdir, backname)
            first = False

            thread = threading.Thread(target=sxbackcmd, args=(imgpath, sx, catpath, ap, outpath, backpath, first), name=f)
            thread.start()
            threads.append(thread)
        for thread in threads:
            thread.join()
    else:
        sys.exit('No applicable flat-fielded .flat.fits files found!')
    print('Sxtrctr back sub complete!')

#%%


def sky_sub(in_path, out_path, sky_path=None, no_flat=False, sex=False):
    if no_flat:
        subtract_sky_and_normalize(in_path,out_path,sky_path)
    elif sex:
        sexback(in_path, out_path)
    else:
        sky_flat_and_normalize(in_path, out_path, sky_path)


def main():
    parser = argparse.ArgumentParser(description='Crops and subtracts sky from files in dir, can also divide out flat')
    parser.add_argument('-no_flat', action='store_true', help='put optional arg if you DIDNT flat field previously')
    parser.add_argument('-sex', action='store_true', help='optional arg to use sxtrctr background sub instead, '
                                                          'outputs .cats and sky subbed imgs')
    parser.add_argument('-in_path', type=str, help='[str] Input imgs path (usually ramps w/ astrometry)')
    parser.add_argument('-out_path', type=str, help='[str] output sky sub image path')
    parser.add_argument('-sky_path', type=str, help='[str] input sky image path (for sky sub')
    args, unknown = parser.parse_known_args()

    sky_sub(args.in_path, args.out_path, args.sky_path, args.no_flat, args.sex)


if __name__ == "__main__":
    main()