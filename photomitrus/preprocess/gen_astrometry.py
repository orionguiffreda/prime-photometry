"""
Runs astrom.net on specified file
"""
#%%
import subprocess
import argparse
import sys
sys.path.insert(0,'C:\PycharmProjects\prime-photometry\photomitrus')
from settings import flist
import os

def astrom(outpath,inlist):
    #inpath = flist()
    #dir = '/mnt/d/PRIME_photometry_test_files/UKIDSS/C4/'
    for f in inlist:
        com = ["solve-field ", '--backend-config /home/prime/miniconda3/pkgs/astrometry-0.94-py39h33f06bc_5/share/astrometry/astrometry.cfg'
            ,' -U none', ' --axy none', ' -S none', ' -M none', ' -R none', ' -B none', ' -O',
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
            ,' --scale-units arcsecperpix', ' --scale-low 0.45',' --scale-high 0.55',' -U none', ' --axy none', ' -S none', ' -M none', ' -R none', ' -B none', ' -O',
            ' -p', ' -z 4', ' -D ' + outpath, ' ' + dir + f] #can change downsample w/ z, or change what outputs here
        s0 = ''
        com = s0.join(com)
        out = subprocess.Popen([com], shell=True)   #include stdout and stderr to suppress astrom.net output
        out.wait()
    print('fields solved, all done!')
#%%
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Runs astrom.net on specified file')
    parser.add_argument('files', nargs='*', type=str, metavar='f', help='Put output path then dir or input list (if applicable)')
    parser.add_argument('-list', action='store_true', help='if downloading from a list of file paths, use this optional arg')
    parser.add_argument('-dir', action='store_true', help='if downloading directly from a dir, use this optional arg,'
                                                          ' include the desired dir after the output path')
    args = parser.parse_args()
    if args.list:
        astrom(args.files[0],args.files[1])
    if args.dir:
        astromdir(args.files[0],args.files[1])
