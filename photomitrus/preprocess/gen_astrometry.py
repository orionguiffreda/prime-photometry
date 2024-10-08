"""
Runs astrom.net on specified file
"""
#%%
import subprocess
import argparse
import sys
from astropy.io import fits
sys.path.insert(0,'C:\PycharmProjects\prime-photometry\photomitrus')
from settings import flist
import os

def astrom(outpath,inlist):
    #inpath = flist()
    #dir = '/mnt/d/PRIME_photometry_test_files/UKIDSS/C4/'
    for f in inlist:
        com = ["solve-field ", '--backend-config /home/prime/miniconda3/pkgs/astrometry-0.94-py39h33f06bc_5/share/astrometry/astrometry.cfg'
            ,' -U none',' --no-verify', ' --axy none', ' -S none', ' -M none', ' -R none', ' -B none', ' -O',
            ' -p', ' -z 4', ' -D ' + outpath, ' ' + f] #can change downsample w/ z, or change what outputs here
        s0 = ''
        com = s0.join(com)
        out = subprocess.Popen([com], shell=True)   #include stdout and stderr to suppress astrom.net output
        out.wait()
    print('fields solved, all done!')

def astromdir(outpath,dir):
    #inpath = flist()
    #dir = '/mnt/d/PRIME_photometry_test_files/UKIDSS/C4/'
    inpath = sorted(os.listdir(dir))
    for f in inpath:
        com = ["solve-field ", '--backend-config /home/prime/miniconda3/pkgs/astrometry-0.94-py39h33f06bc_5/share/astrometry/astrometry.cfg'
            ,' --scale-units arcsecperpix', ' --scale-low 0.45',' --scale-high 0.55',' --no-verify',' -U none', ' --axy none', ' -S none', ' -M none', ' -R none', ' -B none', ' -O',
            ' -p', ' -z 4', ' -D ' + outpath, ' ' + dir + f] #can change downsample w/ z, or change what outputs here
        s0 = ''
        com = s0.join(com)
        out = subprocess.Popen([com], shell=True)   #include stdout and stderr to suppress astrom.net output
        out.wait()
    print('fields solved, all done!')

# astromdir for hard to solve fields, change ra and dec currently
def astromdirhard(outpath,dir,rad,ds):
    #inpath = flist()
    #dir = '/mnt/d/PRIME_photometry_test_files/UKIDSS/C4/'
    inpath = sorted(os.listdir(dir))
    for f in inpath:
        if f.endswith('.fits'):
            img = fits.open(dir+f)
            hdr = img[0].header
            ra = hdr['RA-D']
            dec = hdr['DEC-D']
            com = ["solve-field ", '--backend-config /home/prime/miniconda3/pkgs/astrometry-0.94-py39h33f06bc_5/share/astrometry/astrometry.cfg'
                ,' --scale-units arcsecperpix', ' --scale-low 0.45',' --scale-high 0.55',' --ra ' + str(ra),' --dec ' + str(dec),' --radius '+str(rad),' --cpulimit 60',' -U none',
                   ' --axy list.axy',' -S none',' -M none', ' -R none', ' -B none', ' -O',' -p',' -z '+str(ds), ' -D ' + outpath, ' ' + dir + f] #can change downsample w/ z, or change what outputs here
            s0 = ''
            com = s0.join(com)
            out = subprocess.Popen([com], shell=True)   #include stdout and stderr to suppress astrom.net output
            out.wait()
    print('fields solved, all done!')
#%%
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Runs astrom.net on specified file')
    parser.add_argument('-list', action='store_true', help='if downloading from a list of file paths, use this optional arg')
    parser.add_argument('-hard', action='store_true',
                        help='if cant solve with normal method, use this optional arg')
    parser.add_argument('-output', type=str, help='[str] output path for ramps w/ astrom.')
    parser.add_argument('-input', type=str, help='[str] input path for ramps (or list if using that)')
    parser.add_argument('-rad', type=str, help='*USE ONLY WITH -HARD* radius in deg for astrometry searching, default = 1',default=1)
    parser.add_argument('-ds', type=str, help='*USE ONLY WITH -HARD* amount of downsample, default = 4',default=4)
    args = parser.parse_args()
    if args.list:
        astrom(args.output,args.input)
    elif args.hard:
        astromdirhard(args.output,args.input,args.rad,args.ds)
    else:
        astromdir(args.output,args.input)
