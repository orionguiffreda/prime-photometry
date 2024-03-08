import os
import sys
import subprocess
import argparse
sys.path.insert(0,'/mnt/c/PycharmProjects/prime-photometry/photomitrus/')
from settings import makedirs

#%% directory creation
def makedirectories(parentdir,chip):
    print('generating directories for astrom, sub, sky, and stacked imgs...')
    astromdir, skydir, subdir, stackdir = makedirs(parentdir,chip)
    return astromdir, skydir, subdir, stackdir

#%% mflat creation
"""
def mflat(flatdir, chip):
    print('generating master flats...')
    try:
        command = 'python ./preprocess/gen_flat.py -dir %s -chip %s' % (flatdir, chip)
        print('Executing command: %s' % command)
        rval = subprocess.run(command.split(), check=True)
    except subprocess.CalledProcessError as err:
        print('Could not run with exit error %s'%err)
"""
#%% initial astrometry
def initastrom(astrompath, parentdir, chip):
    os.chdir('/mnt/c/PycharmProjects/prime-photometry/photomitrus/')
    ramppath = parentdir + 'C%i/' % (chip)
    print('running initial astrometry on raw imgs...')
    try:
        command = 'python ./preprocess/gen_astrometry.py -output %s -input %s' % (astrompath, ramppath)
        print('Executing command: %s' % command)
        rval = subprocess.run(command.split(), check=True)
    except subprocess.CalledProcessError as err:
        print('Could not run with exit error %s'%err)

#%% sky gen
def sky(astrompath, skypath, sigma):
    os.chdir('/mnt/c/PycharmProjects/prime-photometry/photomitrus/')
    print('generating sky...')
    try:
        command = 'python ./sky/gen_sky.py -in_path %s -sky_path %s -sigma %i' % (astrompath, skypath, sigma)
        print('Executing command: %s' % command)
        rval = subprocess.run(command.split(), check=True)
    except subprocess.CalledProcessError as err:
        print('Could not run with exit error %s' % err)

#%% sky sub
def skysub(astrompath, subpath, skypath, chip):
    os.chdir('/mnt/c/PycharmProjects/prime-photometry/photomitrus/')
    import fnmatch
    for file in os.listdir(skypath):
        if fnmatch.fnmatch(file, '*.C{}.*'.format(chip)):
            skyfile = file
    print('cropping and subtracting sky...')
    try:
        command = 'python ./sky/sky.py -in_path %s -out_path %s -sky_path %s%s' % (astrompath, subpath, skypath, skyfile)
        print('Executing command: %s' % command)
        rval = subprocess.run(command.split(), check=True)
    except subprocess.CalledProcessError as err:
        print('Could not run with exit error %s' % err)

#%% better astrometry
def astrometry(subpath):
    print('using SXTRCTR and SCAMP to generate better astrometry...')
    try:
        command = 'python ./astrom/astrometry.py -all -path %s' % (subpath)
        print('Executing command: %s' % command)
        rval = subprocess.run(command.split(), check=True)
    except subprocess.CalledProcessError as err:
        print('Could not run with exit error %s' % err)

#%% stacking
def stack(subpath, stackpath):
    print('stacking all images using SWARP...')
    try:
        command = 'python ./stack/stack.py -sub %s -stack %s' % (subpath, stackpath)
        print('Executing command: %s' % command)
        rval = subprocess.run(command.split(), check=True)
    except subprocess.CalledProcessError as err:
        print('Could not run with exit error %s' % err)

#%%
defaults = dict(sigma=4)

#%%
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Automation of the backbone of pipeline, currently processes 1 chip at a time')
    parser.add_argument('-skygen_start',  action='store_true', help='optional flag, starts pipeline at sky gen step')
    parser.add_argument('-skysub_start', action='store_true', help='optional flag, starts pipeline at sky sub step')
    parser.add_argument('-astrom_start', action='store_true', help='optional flag, starts pipeline at sxtrctr / scamp step')
    parser.add_argument('-stack_start', action='store_true',help='optional flag, starts pipeline at swarp step')
    parser.add_argument('-parent', type=str, help='[str], parent directory of outputs, should include folders of the chips ramp data')
    parser.add_argument('-chip', type=int, help='[int], number of detector')
    parser.add_argument('-sigma', type=int, help='[int], sigma value for sky sub sigma clipping, default = 4', default=defaults["sigma"])
    args = parser.parse_args()

    if args.skygen_start:
        astromdir, skydir, subdir, stackdir = makedirectories(args.parent, args.chip)
        sky(astromdir, skydir, args.sigma)
        skysub(astromdir, subdir, skydir, args.chip)
        astrometry(subdir)
        stack(subdir, stackdir)
    elif args.skysub_start:
        astromdir, skydir, subdir, stackdir = makedirectories(args.parent, args.chip)
        skysub(astromdir, subdir, skydir, args.chip)
        astrometry(subdir)
        stack(subdir, stackdir)
    elif args.astrom_start:
        astromdir, skydir, subdir, stackdir = makedirectories(args.parent, args.chip)
        astrometry(subdir)
        stack(subdir, stackdir)
    else:
        astromdir, skydir, subdir, stackdir = makedirectories(args.parent, args.chip)
        initastrom(astromdir, args.parent, args.chip)
        sky(astromdir, skydir, args.sigma)
        skysub(astromdir, subdir, skydir, args.chip)
        astrometry(subdir)
        stack(subdir, stackdir)
