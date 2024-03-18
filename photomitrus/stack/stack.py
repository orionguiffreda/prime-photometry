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

sys.path.insert(0,'C:\PycharmProjects\prime-photometry\photomitrus')
from settings import gen_config_file_name
from settings import gen_mask_file_name

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
    image_fnames = [os.path.join(imgdir, f) for f in os.listdir(imgdir) if f.endswith('.flat.fits')]
    image_fnames.sort()
    header = fits.getheader(image_fnames[-1])
    filter1 = header.get('FILTER1', 'unknown')
    filter2 = header.get('FILTER2', 'unknown')
    save_name = 'coadd.{}-{}.{}-{}.C{}.fits'.format(filter1, filter2, image_fnames[0][-24:-16],
                                                    image_fnames[-1][-24:-16], image_fnames[0][-15])
    print(save_name)
    weight_name = 'weight.{}-{}.{}-{}.C{}.fits'.format(filter1, filter2, image_fnames[0][-24:-16],
                                                    image_fnames[-1][-24:-16], image_fnames[0][-15])
    os.chdir(str(finout))
    sw = gen_config_file_name('default.swarp')
    #save_name = 'coaddback.fits'
    #weight_name = 'coaddbackweight.fits'
    com = ["swarp ", imgdir + '*.flat.fits', ' -c '+sw
           , ' -IMAGEOUT_NAME '+ save_name, ' -WEIGHTOUT_NAME '+weight_name]
    s0 = ''
    com = s0.join(com)
    out = subprocess.Popen([com], shell=True)
    out.wait()
    print('Co-added image created, all done!')


#%%
"""imgdir='/mnt/d/PRIME_photometry_test_files/GRBastrom/'
image_fnames = [os.path.join(imgdir, f) for f in os.listdir(imgdir) if f.endswith('.ramp.new')]
image_fnames.sort()
header = fits.getheader(image_fnames[-1])
filter1 = header.get('FILTER1', 'unknown')
filter2 = header.get('FILTER2', 'unknown')
save_name = 'coadd.{}-{}.{}-{}.C{}.fits'.format(filter1, filter2, image_fnames[0][-19:-11], image_fnames[-1][-19:-11], image_fnames[0][-10])
print(save_name)"""

#%%
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Runs swarp to stack imgs, if using -mask flag, do not keyboard interupt')
    parser.add_argument('-mask', action='store_true', help='optional flag, use if you want to utilize a bad pixel mask')
    parser.add_argument('-sub', type=str, help='[str] Processed images path')
    parser.add_argument('-stack', type=str, help='[str] Output stacked image path')
    parser.add_argument('-parent', type=str, help='[str] Parent directory where all img folders are stored *USE ONLY W/ -MASK FLAG*',
                        default=None)
    parser.add_argument('-chip', type=int, help='[int] Detector chip number *USE ONLY W/ -MASK FLAG*',
                        default=None)
    args = parser.parse_args()
    if args.mask:
        otherdir = badpixmask(args.parent,args.sub,args.chip)
        print('removing temp dir...')
        shutil.rmtree(otherdir)
    swarp(args.sub,args.stack)

#%%
"""swarp('/Users/orion/Desktop/PRIME/GRB/H_band/H_backs/','/Users/orion/Desktop/PRIME/GRB/H_band/H_swarp/')"""