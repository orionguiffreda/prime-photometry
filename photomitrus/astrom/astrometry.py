"""
Runs improved astrom on a given file using sextractor and scamp
"""

import os
import sys
import argparse
import subprocess
import sys

sys.path.insert(0,'/Users/orion/Desktop/PRIME/prime-photometry/photomitrus/')
from settings import gen_config_file_name

sys.path.insert(0,'C:\PycharmProjects\prime-photometry\photomitrus')
from settings import gen_config_file_name

#%%
def sex(imgdir):
    os.chdir(str(imgdir))
    sx = gen_config_file_name('sex.config')
    for f in sorted(os.listdir(str(imgdir))):
        if f.endswith('.fits'):
            pre = os.path.splitext(f)[0]
            ext = os.path.splitext(f)[1]
            com = ["sex ", imgdir + pre + ext, ' -c '+sx,
                   " -CATALOG_NAME " + pre + '.cat']
            s0 = ''
            com = s0.join(com)
            out = subprocess.Popen([com], shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            out.wait()
            print(pre + ext + ' sextracted!')

def scamp(imgdir):
    os.chdir(imgdir)
    sc = gen_config_file_name('scamp.conf')
    for f in sorted(os.listdir(str(imgdir))):  # loop to generate scamp .head headers
        if f.endswith('.cat'):
            pre = os.path.splitext(f)[0]
            ext = os.path.splitext(f)[1]
            com = ["scamp ", imgdir + pre + ext, ' -c '+sc]
            s0 = ''
            com = s0.join(com)
            out = subprocess.Popen([com], shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            out.wait()
            print(pre + ext + ' scamped!')

#%%
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='runs sextractor and scamp on input imgs; generates LDAC .cat and .head files')
    parser.add_argument('-sex', action='store_true', help='if you want to run JUST sextractor')
    parser.add_argument('-scamp', action='store_true', help='if you want to run JUST scamp')
    parser.add_argument('-all', action='store_true', help='if you want to run both')
    parser.add_argument('-path', type=str, help='[str] Images path (currently just dumps .cat '
                                                                       '& .head files in same path)')
    args = parser.parse_args()
    if args.sex:
        sex(args.path)
    if args.scamp:
        scamp(args.path)
    if args.all:
        sex(args.path)
        scamp(args.path)