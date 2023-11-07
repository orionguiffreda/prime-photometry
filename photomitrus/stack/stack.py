"""
Stacks astrometrically calibrated files using swarp
"""
#%%
import os
import argparse
import subprocess

def swarp(indir, finout):
    os.chdir(str(finout))
    com = ["swarp ", indir + '*.new', ' -c /Users/orion/Desktop/PRIME/default.swarp'] #change config file path, also
    s0 = ''
    com = s0.join(com)
    out = subprocess.Popen([com], shell=True)
    out.wait()
    print('Co-added image created, all done!')

#%%
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Runs swarp to stack imgs')
    parser.add_argument('files', nargs='*', type=str, metavar='f', help='Put input file first and output path last')
    args = parser.parse_args()
    swarp(args.files[0],args.files[1])