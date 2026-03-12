"""
Stacks astrometrically calibrated files using swarp
"""
#%%
import os
import shutil
import argparse
import subprocess
import sys
from astropy.io import fits
import astropy.units as u
from astropy.coordinates import SkyCoord
import numpy as np
import fnmatch
import matplotlib.pyplot as plt

# sys.path.insert(0,'C:\PycharmProjects\prime-photometry\photometrus')
from photometrus.astrom.astrometry import sextract, scamp
from photometrus.photometry.photometry import photometry
from photometrus.settings import gen_config_file_name, auto_bulge_detect, gen_mask_file_name
from photometrus.utils.utils import combine_header_and_fits

#%%


def badpixmask(parent, subpath, chip):
    otherdir = os.path.join(parent, 'temp/')
    exists = os.path.exists(otherdir)
    if not exists:
        os.mkdir(otherdir)
    if exists:
        print(otherdir,' exists!')
    mask = gen_mask_file_name('badpixmask_c%i.fits' % chip)
    badmask = fits.getdata(mask)
    badmask = badmask.astype(bool)
    imgdir = subpath
    print('moving imgs to temp dir...')
    for i in sorted(os.listdir(imgdir)):
        if i.endswith('.sky.flat.fits'):
            shutil.move(os.path.join(imgdir, i), os.path.join(otherdir, i))
    print('Applying bad pixel mask...')
    for i in sorted(os.listdir(otherdir)):
        img = fits.open(os.path.join(otherdir, i))
        hdr = img[0].header
        data = img[0].data
        data[~badmask] = np.nan
        fits.writeto(os.path.join(imgdir, i), data, hdr)
    return otherdir


def astromfin(directory, chip):
    if not chip:
        raise ValueError('Specify a chip when using this functionality!')
    print('Re-running astrometry on swarped image! Running sextractor...')
    catpath = swarp_sx(directory, chip)
    print('Applying 4th order scamp fit to stacked image...')
    scamp(directory, swarpcat=catpath)
    print('Combining scamp .head and stacked image...')
    swarp_missfits(directory, chip)
    try:
        os.remove(catpath)
    except Exception as e:
        print(f"Error removing file: {catpath} - {e}")
    print('Absolute astrometry complete!')


def astrom_check(imgdir):
    """
    Check to determine if all images are w/in the same area in the sky (2x dither rad) before attempting stacking.

    Parameters
    ----------
    imgdir: str
        Directory where input images are stored
    """

    image_fnames = [os.path.join(imgdir, f) for f in sorted(os.listdir(imgdir)) if
                    f.endswith('.flat.fits') or f.endswith('.flat.new')]

    image_hdrs = [fits.getheader(img) for img in image_fnames]

    image_ras = [(float(hdr['CRVAL1'])) for hdr in image_hdrs]
    image_decs = [(float(hdr['CRVAL2'])) for hdr in image_hdrs]
    image_coords = SkyCoord(ra=image_ras, dec=image_decs, frame='icrs', unit='degree')

    med_ra = np.nanmedian(image_ras)
    med_dec = np.nanmedian(image_decs)
    med_coords = SkyCoord(ra=med_ra, dec=med_dec, frame='icrs', unit='degree')

    med_dith_rad = np.nanmedian([(float(hdr['DITHRAD'])) for hdr in image_hdrs])
    acc_radius = (2 * med_dith_rad) * u.arcsec

    seps = image_coords.separation(med_coords)

    acc_mask = seps < acc_radius

    if not acc_mask.all():
        print(' Pre-stacking astrom check shows 1 or more images are *NOT* w/in acceptable area!')
        bad_idxs = np.where(~acc_mask)[0]
        bad_imgs = [item for index, item in enumerate(image_fnames) if index in bad_idxs]
        bad_img_names = [os.path.split(img)[1] for img in bad_imgs]
        print(f' Recommend checking quality / astrometry on offending images: {bad_img_names}')

        if len(bad_imgs) > len(image_fnames) // 2:
            raise Exception(f'*WARNING* {len(bad_imgs)}/{len(image_fnames)} images (> 1/2 total # of images) are NOT '
                            f'in acceptable area, examine images!  Is there an issue with image acquisition tracking '
                            f'or astrometry?')

        print(' Renaming offending images to avoid stacking issues...')
        for img_name in bad_imgs:
            os.rename(img_name, img_name.replace('.flat.','.flat.EXCL.'))


def swarp(imgdir, finout):
    print('SWARP Stacking!')
    image_fnames = [os.path.join(imgdir, f) for f in os.listdir(imgdir) if f.endswith('.flat.fits') or f.endswith('.flat.new')]
    image_fnames.sort()
    header = fits.getheader(image_fnames[-1])
    filter1 = header.get('FILTER1', 'unknown')
    filter2 = header.get('FILTER2', 'unknown')
    ext = os.path.splitext(image_fnames[-1])[1]
    if ext == '.new':
        save_name = 'coadd.{}-{}.{}-{}.C{}.fits'.format(filter1, filter2, image_fnames[0][-23:-15],
                                                        image_fnames[-1][-23:-15], image_fnames[0][-14])
        print(save_name)
        weight_name = 'weight.{}-{}.{}-{}.C{}.fits'.format(filter1, filter2, image_fnames[0][-23:-15],
                                                        image_fnames[-1][-23:-15], image_fnames[0][-14])
    else:
        save_name = 'coadd.{}-{}.{}-{}.C{}.fits'.format(filter1, filter2, image_fnames[0][-24:-16],
                                                        image_fnames[-1][-24:-16], image_fnames[0][-15])
        print(save_name)
        weight_name = 'weight.{}-{}.{}-{}.C{}.fits'.format(filter1, filter2, image_fnames[0][-24:-16],
                                                        image_fnames[-1][-24:-16], image_fnames[0][-15])

    os.chdir(str(finout))
    #save_name = 'coaddastr.fits'
    #weight_name = 'coaddastrweight.fits'

    sw = gen_config_file_name('default.swarp')
    com = ["swarp ", os.path.join(imgdir, '*.flat'+ext), ' -c '+sw
           , ' -IMAGEOUT_NAME '+ save_name, ' -WEIGHTOUT_NAME '+weight_name]
    s0 = ''
    com = s0.join(com)
    out = subprocess.Popen([com], shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    out.wait()
    print('Co-added image created, all done!')


def swarp_increm(imgdir, finout, im_num=5):
    batch_size = im_num
    files = [os.path.join(imgdir, f) for f in sorted(os.listdir(imgdir)) if f.endswith('.flat.new') or f.endswith('.flat.fits')]
    total_files = len(files)

    stack_files = [os.path.join(finout, f) for f in sorted(os.listdir(finout)) if f.startswith('coadd') and f.endswith('.fits')]
    if len(stack_files) < 1 - round(total_files / batch_size):

        for i in range(0, total_files, batch_size):
            batch = files[:i + batch_size]
            print(batch)
            print('images taken = ', len(batch))

            image_fnames = batch
            header = fits.getheader(image_fnames[-1])
            filter1 = header.get('FILTER1', 'unknown')
            filter2 = header.get('FILTER2', 'unknown')
            ext = os.path.splitext(image_fnames[-1])[1]
            if ext == '.new':
                save_name = 'coadd.{}-{}.{}-{}.C{}.fits'.format(filter1, filter2, image_fnames[0][-23:-15],
                                                                image_fnames[-1][-23:-15], image_fnames[0][-14])
                print(save_name)
                weight_name = 'weight.{}-{}.{}-{}.C{}.fits'.format(filter1, filter2, image_fnames[0][-23:-15],
                                                                   image_fnames[-1][-23:-15], image_fnames[0][-14])
            elif ext == '.fits':
                save_name = 'coadd.{}-{}.{}-{}.C{}.fits'.format(filter1, filter2, image_fnames[0][-24:-16],
                                                                image_fnames[-1][-24:-16], image_fnames[0][-15])
                print(save_name)
                weight_name = 'weight.{}-{}.{}-{}.C{}.fits'.format(filter1, filter2, image_fnames[0][-24:-16],
                                                                   image_fnames[-1][-24:-16], image_fnames[0][-15])
            os.chdir(str(finout))
            sw = gen_config_file_name('default.swarp')
            #save_name = 'coaddastr.fits'
            #weight_name = 'coaddastrweight.fits'
            image_list = (',').join(image_fnames)
            com = ["swarp ", image_list, ' -c '+sw
                   , ' -IMAGEOUT_NAME '+save_name, ' -WEIGHTOUT_NAME '+weight_name]
            s0 = ''
            com = s0.join(com)
            out = subprocess.Popen([com], shell=True)
            out.wait()
            print('Co-added image created, all done!')

    stack_files = [os.path.join(finout, f) for f in sorted(os.listdir(finout)) if f.startswith('coadd') and f.endswith('.fits')]
    total_stacks = len(stack_files)
    band = fits.getheader(stack_files[0])['FILTER2']
    if len(band) > 1:
        band = 'Z'

    ecsv_files = [f for f in sorted(os.listdir(finout)) if f.startswith('coadd') and f.endswith('.ecsv')]

    lim_mags = []
    for stack in stack_files:
        if len(ecsv_files) < 1 - round(total_files / batch_size):
            print(f'\nRunning photometry on {stack}\n')
            photometry(full_filename=stack, band=band)
        hdr = fits.getheader(stack)
        lim_mag = hdr['lim_mag_auto']
        lim_mags.append(lim_mag)

    exptime = fits.getheader(files[0])['EXPTIMEC']
    single_stack_exp = exptime * im_num
    max_exp = (total_stacks * single_stack_exp)
    exptimes_axis = np.arange(single_stack_exp, max_exp + single_stack_exp, single_stack_exp)

    print('\nGenerating incremental lim mag plot!')
    plt.figure(1, figsize=(10,8))
    plt.plot(exptimes_axis, lim_mags, 'ro')
    plt.grid()
    # plt.yscale('log')
    # plt.xscale('log')
    plt.xticks(exptimes_axis, rotation=45, ha='right')
    plt.xlabel('Exposure Time (s)')
    plt.ylabel(f'{band}MAG_AUTO Limiting Magnitude (AB)')
    plt.title('PRIME Limiting Magnitude Evolution')
    plt.savefig(os.path.join(finout, 'lim_mag_increm.png'), dpi=200)
    plt.clf()
    print('Generated!')

    plt.close('all')


def swarp_alt(imgdir, imout):
    image_fnames = [os.path.join(imgdir, f) for f in os.listdir(imgdir) if f.endswith('.flat.new') or f.endswith('.flat.fits')]
    image_fnames.sort()
    image_fnames_1 = image_fnames[::2]
    image_fnames_2 = image_fnames[1::2]
    print('1st stack = ', len(image_fnames_1))
    print('2nd stack = ', len(image_fnames_2))

    header = fits.getheader(image_fnames[-1])
    filter1 = header.get('FILTER1', 'unknown')
    filter2 = header.get('FILTER2', 'unknown')
    ext = os.path.splitext(image_fnames_1[-1])[1]
    if ext == '.new':
        save_name1 = 'coadd.{}-{}.{}-{}.C{}.fits'.format(filter1, filter2, image_fnames_1[0][-23:-15],
                                                        image_fnames_1[-1][-23:-15], image_fnames_1[0][-14])
        print(save_name1)
        weight_name1 = 'weight.{}-{}.{}-{}.C{}.fits'.format(filter1, filter2, image_fnames_1[0][-23:-15],
                                                           image_fnames_1[-1][-23:-15], image_fnames_1[0][-14])

        save_name2 = 'coadd.{}-{}.{}-{}.C{}.fits'.format(filter1, filter2, image_fnames_2[0][-23:-15],
                                                        image_fnames_2[-1][-23:-15], image_fnames_2[0][-14])
        print(save_name2)
        weight_name2 = 'weight.{}-{}.{}-{}.C{}.fits'.format(filter1, filter2, image_fnames_2[0][-23:-15],
                                                           image_fnames_2[-1][-23:-15], image_fnames_2[0][-14])
    elif ext == '.fits':
        save_name1 = 'coadd.{}-{}.{}-{}.C{}.fits'.format(filter1, filter2, image_fnames_1[0][-24:-16],
                                                        image_fnames_1[-1][-24:-16], image_fnames_1[0][-15])
        print(save_name1)
        weight_name1 = 'weight.{}-{}.{}-{}.C{}.fits'.format(filter1, filter2, image_fnames_1[0][-24:-16],
                                                           image_fnames_1[-1][-24:-16], image_fnames_1[0][-15])

        save_name2 = 'coadd.{}-{}.{}-{}.C{}.fits'.format(filter1, filter2, image_fnames_2[0][-24:-16],
                                                        image_fnames_2[-1][-24:-16], image_fnames_2[0][-15])
        print(save_name2)
        weight_name2 = 'weight.{}-{}.{}-{}.C{}.fits'.format(filter1, filter2, image_fnames_2[0][-24:-16],
                                                           image_fnames_2[-1][-24:-16], image_fnames_2[0][-15])
    os.chdir(str(imout))
    sw = gen_config_file_name('default.swarp')
    # save_name = 'coaddastr.fits'
    # weight_name = 'coaddastrweight.fits'
    image_list1 = (',').join(image_fnames_1)
    image_list2 = (',').join(image_fnames_2)
    com = ["swarp ", image_list1, ' -c ' + sw
        , ' -IMAGEOUT_NAME ' + save_name1, ' -WEIGHTOUT_NAME ' + weight_name1]
    com = ''.join(com)
    out = subprocess.Popen([com], shell=True)
    out.wait()
    print('First set co-added image created, all done!')

    com2 = ["swarp ", image_list2, ' -c ' + sw
        , ' -IMAGEOUT_NAME ' + save_name2, ' -WEIGHTOUT_NAME ' + weight_name2]
    com2 = ''.join(com2)
    out2 = subprocess.Popen([com2], shell=True)
    out2.wait()
    print('2nd set co-added image created, all done!')


def swarp_sx(imgpath, chip):
    # imgpath can be either full file path to stack or directory where stack is in

    if os.path.isdir(imgpath):
        os.chdir(str(imgpath))
        bulge = auto_bulge_detect(imgpath)
        if bulge:
            sx = gen_config_file_name('bulge_new.config')
            ap = gen_config_file_name('tempsource.param')
        else:
            sx = gen_config_file_name('sex.config')
            ap = gen_config_file_name('astrom_coadd.param')

        stackimg = [f for f in sorted(os.listdir(imgpath)) if fnmatch.fnmatch(f, 'coadd.*.C%i.fits' % chip)]
        coaddimg = stackimg[0]
        catname = coaddimg.replace('.fits', '.cat')
        catpath = os.path.join(imgpath, catname)
        weightname = 'weight' + coaddimg[5:]
    elif os.path.isfile(imgpath):
        coaddimgdir, coaddimg = os.path.split(imgpath)
        os.chdir(str(coaddimgdir))
        bulge = auto_bulge_detect(coaddimgdir)
        if bulge:
            sx = gen_config_file_name('bulge_new.config')
            ap = gen_config_file_name('tempsource.param')
        else:
            sx = gen_config_file_name('sex.config')
            ap = gen_config_file_name('astrom_coadd.param')

        catname = coaddimg.replace('.fits', '.cat')
        catpath = os.path.join(coaddimgdir, catname)
        weightname = 'weight' + coaddimg[5:]
    else:
        print('Must specify a directory to stacked image or full file path to stacked image! '
              'Cannot continue with stack image sextraction!')
        catname = catpath = weightname = None
    if catname:
        if os.path.isfile(weightname):
            print('Including weight map!')
            command = ('sex %s -c %s -CATALOG_NAME %s -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE %s -PARAMETERS_NAME %s' %
                       (coaddimg, sx, catname, weightname, ap))
            # print('Executing command: %s' % command)
            subprocess.run(command.split(), check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        else:
            command = ('sex %s -c %s -CATALOG_NAME %s -PARAMETERS_NAME %s' %
                       (coaddimg, sx, catname, ap))
            # print('Executing command: %s' % command)
            subprocess.run(command.split(), check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return catpath


def swarp_missfits(imgpath, chip):
    # imgpath can be either full file path to stack or directory where stack is in
    if os.path.isdir(imgpath):
        stackimg = [f for f in sorted(os.listdir(imgpath)) if fnmatch.fnmatch(f, 'coadd.*.C%i.fits' % chip)]
        stackhdr = [f for f in sorted(os.listdir(imgpath)) if fnmatch.fnmatch(f, 'coadd.*.C%i.head' % chip)]
        combine_header_and_fits(stackhdr[0], stackimg[0], remove_header_file=True)
    else:
        stackimg = imgpath
        stackhdr = imgpath.replace('.fits', '.head')
        combine_header_and_fits(stackhdr, stackimg, remove_header_file=True)


#%%


def stack(subpath, stackpath, chip, num=5, no_astrom=False, astrom_only=False, increm=False, alt=False, mosaic=False):
    # if args.mask:
        # otherdir = badpixmask(args.parent,args.sub,args.chip)
        # print('removing temp dir...')
        # shutil.rmtree(otherdir)
    if not mosaic:
        astrom_check(subpath)

    if no_astrom:
        swarp(subpath, stackpath)
    elif astrom_only:
        astromfin(stackpath, chip)
    elif increm:
        swarp_increm(subpath, stackpath, num)
    elif alt:
        swarp_alt(subpath, stackpath)
        astromfin(stackpath, chip)
    else:
        if not chip:
            raise ValueError('Remember to specify chip number using default stacking behavior!  It is required for '
                             'absolute astrometry check!')
        swarp(subpath, stackpath)
        astromfin(stackpath, chip)


def main():
    parser = argparse.ArgumentParser(description='Runs swarp to stack imgs, then reruns astrometry for improved wcs *NOTE* if using -mask flag, do not keyboard interrupt')
    # parser.add_argument('-mask', action='store_true', help='optional flag, use if you want to utilize a bad pixel mask')
    parser.add_argument('-no_astrom', action='store_true', help='optional flag, use if you just want the swarped image, not the image with improved astrometry')
    parser.add_argument('-astrom_only', action='store_true',
                        help='optional flag, use if you already have the swarped image, but want improved astrometry')
    parser.add_argument('-increm', action='store_true',
                        help='optional flag, use if you want to generate a stacked img from increments of images, ex. 5 stack, then 10 stack, etc.')
    parser.add_argument('-alt', action='store_true',
                        help='create 2 stacked images from 1 set of data, alternating images used')
    parser.add_argument('-mosaic', action='store_true',
                        help='Use this flag if you are attempting to make a large mosaic, will disable the default '
                             'image location screening')
    parser.add_argument('-sub', type=str, help='[str] Processed images path')
    parser.add_argument('-stack', type=str, help='[str] Output stacked image path')
    parser.add_argument('-num', type=int, help='*USE ONLY W/ -INCREM* [int] # of imgs to increment by', default=5)
    #parser.add_argument('-parent', type=str, help='*USE ONLY W/ -MASK FLAG* [str] Parent directory where all img folders are stored', default=None)
    parser.add_argument('-chip', type=int, help='*USE ONLY W/O -no_astrom FLAG* [int] Detector chip number', default=None)
    args, unknown = parser.parse_known_args()

    stack(args.sub, args.stack, args.chip, args.num, args.no_astrom, args.astrom_only, args.increm, args.alt, args.mosaic)


if __name__ == "__main__":
    main()
