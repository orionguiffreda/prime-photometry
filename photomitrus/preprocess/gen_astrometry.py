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

def astrom(outpath):
    inpath = flist()
    #dir = '/mnt/d/PRIME_photometry_test_files/C4_GRBsub/'
    #inpath = sorted(os.listdir(dir))
    for f in inpath:
        com = ["solve-field ", '--backend-config /home/prime/miniconda3/pkgs/astrometry-0.94-py39h33f06bc_5/share/astrometry/astrometry.cfg'
            ,' -U none', ' --axy none', ' -S none', ' -M none', ' -R none', ' -B none', ' -O',
            ' -p', ' -z 4', ' -D ' + outpath, ' ' + f] #can change downsample w/ z, or change what outputs here
        s0 = ''
        com = s0.join(com)
        out = subprocess.Popen([com], shell=True)   #include stdout and stderr to suppress astrom.net output
        out.wait()
    print('fields solved, all done!')

#%%
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Runs astrom.net on specified file')
    parser.add_argument('files', nargs='*', type=str, metavar='f', help='Put output path')
    args = parser.parse_args()
    astrom(args.files[0])
