import os
import sys
import subprocess
import argparse
import fnmatch
# sys.path.insert(0,'/mnt/c/PycharmProjects/prime-photometry/photomitrus/')

from photomitrus.settings import makedirs
from photomitrus.settings import makedirsFF
from photomitrus.settings import gen_pipeline_file_name
from photomitrus.settings import gen_mflat_file_name


#%% directory creation
def makedirectories(parentdir,chip):
    print('generating directories for astrom, sub, sky, and stacked imgs...')
    astromdir, skydir, subdir, stackdir = makedirs(parentdir,chip)
    return astromdir, skydir, subdir, stackdir


def makedirectoriesFF(parentdir,chip):
    print('generating FF directory')
    FFdir = makedirsFF(parentdir,chip)
    return FFdir

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


def initastrom(astrompath, parentdir, chip=None):
    os.chdir(gen_pipeline_file_name())
    if chip:
        ramppath = parentdir + 'C%i/' % (chip)
        print('running initial astrometry on ramp imgs...')
        try:
            command = 'python ./preprocess/gen_astrometry.py -hard -output %s -input %s -rad 1 -ds 3' % (astrompath, ramppath)
            print('Executing command: %s' % command)
            subprocess.run(command.split(), check=True)
        except subprocess.CalledProcessError as err:
            print('Could not run with exit error %s'%err)
    if not chip:
        print('running astrometry.net on subbed imgs...')
        try:
            command = 'python ./preprocess/gen_astrometry.py -hard -output %s -input %s -rad 1 -ds 3' % (astrompath, astrompath)
            print('Executing command: %s' % command)
            subprocess.run(command.split(), check=True)
        except subprocess.CalledProcessError as err:
            print('Could not run with exit error %s'%err)

#%% angle astrometry


def astrom_angle(astrompath, parentdir, chip, rot_val=None):
    os.chdir(gen_pipeline_file_name())
    ramppath = parentdir + 'C%i/' % (chip)
    print('running initial astrometry on ramp imgs...')
    if rot_val:
        try:
            command = 'python ./preprocess/astromangle_new.py -input %s -output %s -rot_val %s' % (ramppath, astrompath, rot_val)
            print('Executing command: %s' % command)
            subprocess.run(command.split(), check=True)
        except subprocess.CalledProcessError as err:
            print('Could not run with exit error %s'%err)
    else:
        try:
            command = 'python ./preprocess/astromangle_new.py -input %s -output %s' % (ramppath, astrompath)
            print('Executing command: %s' % command)
            subprocess.run(command.split(), check=True)
        except subprocess.CalledProcessError as err:
            print('Could not run with exit error %s' % err)

#%% flat fielding


def flatfield(astrompath,FFpath,band,chip):
    os.chdir(gen_pipeline_file_name())
    print('using master flat to flat field ramp imgs..')
    flatpath = gen_mflat_file_name(band,chip)
    try:
        command = 'python ./preprocess/flatfield.py -in_path %s -out_path %s -flat_path %s' % (astrompath,FFpath,flatpath)
        print('Executing command: %s' % command)
        subprocess.run(command.split(), check=True)
    except subprocess.CalledProcessError as err:
        print('Could not run with exit error %s'%err)

#%% sky gen


def sky(astrompath, skypath, sigma):
    os.chdir(gen_pipeline_file_name())
    print('generating sky...')
    FFstring = '_FF'
    if FFstring in astrompath:
        try:
            command = 'python ./sky/gen_sky.py -flat -in_path %s -sky_path %s -sigma %i' % (astrompath, skypath, sigma)
            print('Executing command: %s' % command)
            subprocess.run(command.split(), check=True)
        except subprocess.CalledProcessError as err:
            print('Could not run with exit error %s' % err)
    else:
        try:
            command = 'python ./sky/gen_sky.py -in_path %s -sky_path %s -sigma %i' % (astrompath, skypath, sigma)
            print('Executing command: %s' % command)
            subprocess.run(command.split(), check=True)
        except subprocess.CalledProcessError as err:
            print('Could not run with exit error %s' % err)

#%% sky sub


def skysub(astrompath, subpath, skypath, chip):
    os.chdir(gen_pipeline_file_name())
    import fnmatch
    for file in os.listdir(skypath):
        if fnmatch.fnmatch(file, '*.C{}.*'.format(chip)):
            skyfile = file
    print('cropping and subtracting sky...')
    FFstring = '_FF'
    if FFstring in astrompath:
        try:
            command = 'python ./sky/sky.py -flat -in_path %s -out_path %s -sky_path %s%s' % (astrompath, subpath, skypath, skyfile)
            print('Executing command: %s' % command)
            subprocess.run(command.split(), check=True)
        except subprocess.CalledProcessError as err:
            print('Could not run with exit error %s' % err)
    else:
        try:
            command = 'python ./sky/sky.py -in_path %s -out_path %s -sky_path %s%s' % (astrompath, subpath, skypath, skyfile)
            print('Executing command: %s' % command)
            subprocess.run(command.split(), check=True)
        except subprocess.CalledProcessError as err:
            print('Could not run with exit error %s' % err)
#%%sex sky sub


def sexskysub(astrompath,subpath):
    os.chdir(gen_pipeline_file_name())
    try:
        command = 'python ./sky/sky.py -sex -in_path %s -out_path %s' % (astrompath, subpath)
        print('Executing command: %s' % command)
        subprocess.run(command.split(), check=True)
    except subprocess.CalledProcessError as err:
        print('Could not run with exit error %s' % err)

#%% astrometry shift


def shift(subpath, band):
    os.chdir(gen_pipeline_file_name())
    print('Shifting astrometry...')
    all_fits = [f for f in sorted(os.listdir(subpath)) if f.endswith('.flat.fits')]
    if len(all_fits) >= 100:
        imgname = all_fits[0]
        try:
            # change back to large once fully tested
            command = 'python ./astrom/astrom_shift.py -remove -segment -pipeline -dir %s -imagename %s -band %s' % (subpath, imgname, band)
            print('Executing command: %s' % command)
            subprocess.run(command.split(), check=True)
        except subprocess.CalledProcessError as err:
            print('Could not run with exit error %s' % err)
    else:
        imgname = all_fits[0]
        try:
            command = 'python ./astrom/astrom_shift.py -remove -segment -pipeline -dir %s -imagename %s -band %s' % (subpath, imgname, band)
            print('Executing command: %s' % command)
            subprocess.run(command.split(), check=True)
        except subprocess.CalledProcessError as err:
            print('Could not run with exit error %s' % err)

#%% better astrometry


def astrometry(subpath,sex=None):
    os.chdir(gen_pipeline_file_name())
    print('using SXTRCTR and SCAMP to generate better astrometry...')
   # if sex:
   #     try:
    #        command = 'python ./astrom/astrometry.py -scamp -path %s' % (subpath)
    #        print('Executing command: %s' % command)
    #        rval = subprocess.run(command.split(), check=True)
    #    except subprocess.CalledProcessError as err:
    #        print('Could not run with exit error %s' % err)
    #else:
    try:
        command = 'python ./astrom/astrometry.py -all -path %s' % (subpath)
        print('Executing command: %s' % command)
        subprocess.run(command.split(), check=True)
    except subprocess.CalledProcessError as err:
        print('Could not run with exit error %s' % err)

#%% stacking


def stack(subpath, stackpath,chip):
    os.chdir(gen_pipeline_file_name())
    print('stacking all images using SWARP...')
    try:
        command = 'python ./stack/stack.py -sub %s -stack %s -chip %i' % (subpath, stackpath, chip)
        print('Executing command: %s' % command)
        subprocess.run(command.split(), check=True)
    except subprocess.CalledProcessError as err:
        print('Could not run with exit error %s' % err)

#%% packing compression


def fpack(stackpath,chip):
    os.chdir(gen_pipeline_file_name())
    print('compressing stacked images with fpack...')
    instack = sorted(os.listdir(stackpath))
    stackimg = []
    for f in instack:
        if fnmatch.fnmatch(f, 'coadd.Open-*.C%i.fits' % chip):
            stackimg.append(f)
    for f in stackimg:
        try:
            command = 'fpack -D -Y %s%s' % (stackpath,f)
            subprocess.run(command.split(), check=True)
        except subprocess.CalledProcessError as err:
            print('Could not run with exit error %s' % err)
    print('Stack compressed!')

#%% astromnet refining


def astromnet_refine(subdir):
    initastrom(subdir, subdir)
    oldlist = [f for f in sorted(os.listdir(subdir)) if f.endswith('.flat.fits')]
    newlist = [j for j in sorted(os.listdir(subdir)) if j.endswith('.flat.new')]
    errnum = len(oldlist) * 0.2
    if len(newlist) < len(oldlist) - errnum:
        print(
            'Not enough success with new astrometry! (%i fields) Continuing with initial astrometry...' % len(
                newlist))
        for f in newlist:
            os.remove(subdir + f)
    else:
        for f in oldlist:
            os.remove(subdir + f)
        print('%i files refined! removed old fits files' % len(newlist))

#%%


defaults = dict(sigma=4)


def master(
        parentdir, chip, band, sigma, rot_val=None, no_ff=False, no_shift=False, sex=False, compress=False,
        net_refine=False
):

    if no_ff:
        astromdir, skydir, subdir, stackdir = makedirectories(parentdir, chip)
        astrom_angle(astromdir, parentdir, chip, rot_val)
        sky(astromdir, skydir, sigma)
        if sex:
            sexskysub(astromdir, subdir)
        else:
            skysub(astromdir, subdir, skydir, chip)
        astrometry(subdir)
        stack(subdir, stackdir, chip)
    else:
        astromdir, skydir, subdir, stackdir = makedirectories(parentdir, chip)
        FFdir = makedirectoriesFF(parentdir, chip)
        astrom_angle(astromdir, parentdir, chip, rot_val)
        flatfield(astromdir, FFdir, band, chip)
        if sex:
            pass
        else:
            sky(FFdir, skydir, sigma)
        if sex:
            sexskysub(FFdir, subdir)
        else:
            skysub(FFdir, subdir, skydir, chip)
        if net_refine:
            astromnet_refine(subdir)
        else:
            if not no_shift:
                shift(subdir, band)
            else:
                pass
        astrometry(subdir)
        stack(subdir, stackdir, chip)
        if compress:
            fpack(stackdir, chip)


def main():
    parser = argparse.ArgumentParser(description='Automation of the backbone of pipeline, currently processes 1 chip at a time')
    parser.add_argument('-parent', type=str, help='[str], parent directory of outputs, should include folders of the chips ramp data')
    parser.add_argument('-chip', type=int, help='[int], number of detector')
    parser.add_argument('-band', type=str, help='*NOT NECESSARY UNLESS USING -FF* [str], band of images, ex. "J"',default=None)
    parser.add_argument('-sigma', type=int, help='[int], sigma value for sky sub sigma clipping, default = 4', default=defaults["sigma"])
    parser.add_argument('-rot_val', type=float, help='[float] optional, put in your rot angle in deg,'
                                                     ' if you had a non-default rotation angle in your obs'
                                                     ' (default = 48 deg or 172800")', default=None)
    parser.add_argument('-no_FF', action='store_true', help='optional flag, does not use flat fielding in pipeline')
    parser.add_argument('-sex', action='store_true', help='optional flag, to utlize sextractor background subtraction instead, do not currently use!')
    parser.add_argument('-compress', action='store_true',help='optional flag, use fpack to compress stacked images')
    parser.add_argument('-no_shift', action='store_true', help='optional flag, STOPS use of astrometric shifting script')
    # parser.add_argument('-skygen_start',  action='store_true', help='optional flag, starts pipeline at sky gen step')
    # parser.add_argument('-skysub_start', action='store_true', help='optional flag, starts pipeline at sky sub step')
    # parser.add_argument('-astrom_start', action='store_true', help='optional flag, starts pipeline at sxtrctr / scamp step')
    # parser.add_argument('-stack_start', action='store_true',help='optional flag, starts pipeline at swarp step')
    parser.add_argument('-net_refine', action='store_true', help='optional flag, used to automatically refine astrometry using '
                                                             'astrometry.net')
    args, unknown = parser.parse_known_args()

    master(args.parent, args.chip, args.band, args.sigma, args.rot_val, args.no_FF, args.no_shift, args.sex,
           args.compress, args.net_refine)


if __name__ == "__main__":
    main()
