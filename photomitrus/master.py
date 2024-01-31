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

#%% initial astrometry
def initastrom(astrompath, ramppath):
    os.chdir('/mnt/c/PycharmProjects/prime-photometry/photomitrus/')
    print('running initial astrometry on raw imgs...')
    try:
        command = 'python ./preprocess/gen_astrometry.py -dir %s %s' % (astrompath, ramppath)
        print('Executing command: %s' % command)
        rval = subprocess.run(command.split(), check=True)
    except subprocess.CalledProcessError as err:
        print('Could not run with exit error %s'%err)

#%% sky gen
def sky(astrompath, skypath, sigma):
    print('generating sky...')
    try:
        command = 'python ./sky/gen_sky.py %s %s %i' % (astrompath, skypath, sigma)
        print('Executing command: %s' % command)
        rval = subprocess.run(command.split(), check=True)
    except subprocess.CalledProcessError as err:
        print('Could not run with exit error %s' % err)

#%% sky sub
def skysub(astrompath, subpath, skypath, chip):
    import fnmatch
    for file in os.listdir(skypath):
        if fnmatch.fnmatch(file, '*.C{}.*'.format(chip)):
            skyfile = file
    print('cropping and subtracting sky...')
    try:
        command = 'python ./sky/sky.py %s %s %s%s' % (astrompath, subpath, skypath, skyfile)
        print('Executing command: %s' % command)
        rval = subprocess.run(command.split(), check=True)
    except subprocess.CalledProcessError as err:
        print('Could not run with exit error %s' % err)

#%% better astrometry
def astrometry(subpath):
    print('using SXTRCTR and SCAMP to generate better astrometry...')
    try:
        command = 'python ./astrom/astrometry.py -all %s' % (subpath)
        print('Executing command: %s' % command)
        rval = subprocess.run(command.split(), check=True)
    except subprocess.CalledProcessError as err:
        print('Could not run with exit error %s' % err)

#%% stacking
def stack(subpath, stackpath):
    print('stacking all images using SWARP...')
    try:
        command = 'python ./stack/stack.py %s %s' % (subpath, stackpath)
        print('Executing command: %s' % command)
        rval = subprocess.run(command.split(), check=True)
    except subprocess.CalledProcessError as err:
        print('Could not run with exit error %s' % err)

#%%
defaults = dict(sigma=4)

#%%
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Automation of the backbone of pipeline, currently processes 1 chip at a time')
    parser.add_argument('-parent', type=str, help='[str], parent directory of outputs, should include astrom, sky, sub, and stack dirs')
    parser.add_argument('-ramp',type=str,  help='[str], directory of ramp files')
    parser.add_argument('-chip', type=int, help='[int], number of detector')
    parser.add_argument('-sigma', type=int, help='[int], sigma value for sky sub sigma clipping, default = 4', default=defaults["sigma"])
    args = parser.parse_args()

    astromdir, skydir, subdir, stackdir = makedirectories(args.parent,args.chip)
    initastrom(astromdir,args.ramp)
    sky(astromdir,skydir,args.sigma)
    skysub(astromdir,subdir,skydir,args.chip)
    astrometry(subdir)
    stack(subdir,stackdir)
