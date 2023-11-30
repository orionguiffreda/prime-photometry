"""
Runs astrometry.net on specified file
"""
#%%
import subprocess
import argparse
def astrom(inpath, outpath):
    com = ["solve-field ", ' -U none', ' --axy none', ' -S none', ' -M none', ' -R none', ' -B none', ' -O',
           ' -p', ' -z 2', ' -D ' + outpath, ' ' + inpath + '*.fits'] #can change downsample w/ z, or change what outputs here
    s0 = ''
    com = s0.join(com)
    out = subprocess.Popen([com], shell=True)   #include stdout and stderr to suppress astrometry.net output
    out.wait()
    print('fields solved, all done!')

#%%
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Runs astrometry.net on specified file')
    parser.add_argument('files', nargs='*', type=str, metavar='f', help='Put input file first and output path last')
    args = parser.parse_args()
    astrom(args.files[0],args.files[1])