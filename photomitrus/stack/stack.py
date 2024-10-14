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
import numpy as np

# sys.path.insert(0,'C:\PycharmProjects\prime-photometry\photomitrus')
from photomitrus.settings import gen_config_file_name
from photomitrus.settings import gen_mask_file_name

#%%


def badpixmask(parent, subpath, chip):
    otherdir = parent+'temp/'
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
            shutil.move(imgdir+i,otherdir+i)
    print('Applying bad pixel mask...')
    for i in sorted(os.listdir(otherdir)):
        img = fits.open(otherdir + i)
        hdr = img[0].header
        data = img[0].data
        data[~badmask] = np.nan
        fits.writeto(imgdir + i, data, hdr)
    return otherdir


def swarp(imgdir, finout):
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
    com = ["swarp ", imgdir + '*.flat'+ext, ' -c '+sw
           , ' -IMAGEOUT_NAME '+ save_name, ' -WEIGHTOUT_NAME '+weight_name]
    s0 = ''
    com = s0.join(com)
    out = subprocess.Popen([com], shell=True)
    out.wait()
    print('Co-added image created, all done!')


def swarp_increm(imgdir, finout,im_num):
    batch_size = im_num
    files = [os.path.join(imgdir, f) for f in sorted(os.listdir(imgdir)) if f.endswith('.flat.new') or f.endswith('.flat.fits')]
    total_files = len(files)

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
                                                            image_fnames[-1][-23:-15], image_fnames[0][-13])
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
               , ' -IMAGEOUT_NAME '+ save_name, ' -WEIGHTOUT_NAME '+weight_name]
        s0 = ''
        com = s0.join(com)
        out = subprocess.Popen([com], shell=True)
        out.wait()
        print('Co-added image created, all done!')


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


def astromfin(directory,chip):
    print('Re-running astrometry on swarped image!')
    import fnmatch
    instack = sorted(os.listdir(dir))
    stackimg = []
    for f in instack:
        if fnmatch.fnmatch(f, 'coadd.*.C%i.fits' % chip):
            stackimg.append(f)
    for f in stackimg:
        pre = os.path.splitext(f)[0]
        try:
            command = ('solve-field '
                       '--backend-config /home/prime/miniconda3/pkgs/astrometry-0.94-py39h33f06bc_5/share/astrometry/astrometry.cfg '
                       '--scale-units arcsecperpix --scale-low 0.45 --scale-high 0.55 --no-verify -U none --axy none '
                       '-S none -M none -R none -B none -O -p -z 4 -D %s %s') % (
                          directory, directory + f)
            print('Executing command: %s' % command)
            subprocess.run(command.split(), check=True)
        except subprocess.CalledProcessError as err:
            print('Could not run with exit error %s' % err)
        # renaming
        if os.path.isfile(directory+pre+'.wcs'):
            os.remove(directory+f)
            os.rename(directory+pre+'.new',directory+f)
            print('final image generated, original stack discarded!')
        else:
            print('New astrometry on stacked image failed, defaulting to original...')


#%%


def stack(subpath, stackpath, chip, num=5, no_astrom=False, astrom_only=False, increm=False, alt=False):
    # if args.mask:
        # otherdir = badpixmask(args.parent,args.sub,args.chip)
        # print('removing temp dir...')
        # shutil.rmtree(otherdir)
    if no_astrom:
        swarp(subpath, stackpath)
    elif astrom_only:
        astromfin(stackpath,chip)
    elif increm:
        swarp_increm(subpath,stackpath,num)
        astromfin(stackpath,chip)
    elif alt:
        swarp_alt(subpath,stackpath)
        astromfin(stackpath,chip)
    else:
        swarp(subpath,stackpath)
        astromfin(stackpath,chip)


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
    parser.add_argument('-sub', type=str, help='[str] Processed images path')
    parser.add_argument('-stack', type=str, help='[str] Output stacked image path')
    parser.add_argument('-num', type=int, help='*USE ONLY W/ -INCREM* [int] # of imgs to increment by')
    #parser.add_argument('-parent', type=str, help='*USE ONLY W/ -MASK FLAG* [str] Parent directory where all img folders are stored', default=None)
    parser.add_argument('-chip', type=int, help='*USE ONLY W/O -no_astrom FLAG* [int] Detector chip number', default=None)
    args, unknown = parser.parse_known_args()

    stack(args.sub, args.stack, args.chip, args.num, args.no_astrom, args.astrom_only, args.increm, args.alt)


if __name__ == "__main__":
    main()
