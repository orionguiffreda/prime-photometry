"""
Stacks astrometrically calibrated files using swarp
"""
#%%
import os
import argparse
import subprocess
import sys
from astropy.io import fits

sys.path.insert(0,'C:\PycharmProjects\prime-photometry\photomitrus')
from settings import gen_config_file_name

#%%
def swarp(imgdir, finout):
    image_fnames = [os.path.join(imgdir, f) for f in os.listdir(imgdir) if f.endswith('.sky.flat.fits')]
    image_fnames.sort()
    header = fits.getheader(image_fnames[-1])
    filter1 = header.get('FILTER1', 'unknown')
    filter2 = header.get('FILTER2', 'unknown')
    save_name = 'coadd.{}-{}.{}-{}.C{}.fits'.format(filter1, filter2, image_fnames[0][-24:-16],
                                                    image_fnames[-1][-24:-16], image_fnames[0][-15])
    print(save_name)
    weight_name =  'weight.{}-{}.{}-{}.C{}.fits'.format(filter1, filter2, image_fnames[0][-24:-16],
                                                    image_fnames[-1][-24:-16], image_fnames[0][-15])
    os.chdir(str(finout))
    sw = gen_config_file_name('default.swarp')
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
    parser = argparse.ArgumentParser(description='Runs swarp to stack imgs')
    parser.add_argument('paths', type=str, metavar='p', help='Put input file path first and output path last')
    args = parser.parse_args()
    swarp(args.paths[0],args.paths[1])

#%%
"""swarp('/Users/orion/Desktop/PRIME/GRB/H_band/H_backs/','/Users/orion/Desktop/PRIME/GRB/H_band/H_swarp/')"""