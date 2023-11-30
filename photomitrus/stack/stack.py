"""
Stacks astrometrically calibrated files using swarp
"""
#%%
import os
import argparse
import subprocess
from astropy.io import fits
#%%
def swarp(imgdir, finout):
    image_fnames = [os.path.join(imgdir, f) for f in os.listdir(imgdir) if f.endswith('.sky.flat.fits')]
    image_fnames.sort()
    header = fits.getheader(image_fnames[-1])
    filter1 = header.get('FILTER1', 'unknown')
    filter2 = header.get('FILTER2', 'unknown')
    save_name = 'coadd.{}-{}.{}-{}.C{}.fits'.format(filter1, filter2, image_fnames[0][-28:-20], image_fnames[-1][-28:-20], image_fnames[0][-19])
    print(save_name)
    os.chdir(str(finout))
    com = ["swarp ", imgdir + '*.flat.fits', ' -c /Users/orion/Desktop/PRIME/default.swarp'
           , ' -IMAGEOUT_NAME '+ save_name] #change config file path, also
    s0 = ''
    com = s0.join(com)
    out = subprocess.Popen([com], shell=True)
    out.wait()
    print('Co-added image created, all done!')

#%%
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Runs swarp to stack imgs')
    parser.add_argument('paths', nargs='*', type=str, metavar='p', help='Put input file path first and output path last')
    args = parser.parse_args()
    swarp(args.paths[0],args.paths[1])

#%%
"""swarp('/Users/orion/Desktop/PRIME/GRB/H_band/H_backs/','/Users/orion/Desktop/PRIME/GRB/H_band/H_swarp/')"""