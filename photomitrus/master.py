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
from photomitrus.preprocess import astromangle_new
from photomitrus.preprocess import gen_astrometry
from photomitrus.preprocess import flatfield
from photomitrus.sky import gen_sky
from photomitrus.sky import sky_sub
from photomitrus.astrom import astrom_shift
from photomitrus.astrom import astrometry
from photomitrus.stack import stack


# %% directory creation
def makedirectories(parentdir, chip):
    print('generating directories for astrom, sub, sky, and stacked imgs...')
    astromdir, skydir, subdir, stackdir = makedirs(parentdir, chip)
    return astromdir, skydir, subdir, stackdir


def makedirectoriesFF(parentdir, chip):
    print('generating FF directory')
    FFdir = makedirsFF(parentdir, chip)
    return FFdir


# %% mflat creation
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


# %% initial astrometry


def initastrom(astrompath, parentdir, chip=None):
    os.chdir(gen_pipeline_file_name())
    if chip:
        ramppath = os.path.join(parentdir, 'C%i' % chip)
        print('running initial astrometry on ramp imgs...')

        print('\nEquivalent argparse cmd: python ./preprocess/gen_astrometry.py -output %s -input %s' % (
        astrompath, ramppath))

        gen_astrometry.gen_astrom(output=astrompath, input=ramppath)
    if not chip:
        print('running astrometry.net on subbed imgs...')

        print('\nEquivalent argparse cmd: python ./preprocess/gen_astrometry.py -output %s -input %s' % (
        astrompath, astrompath))

        gen_astrometry.gen_astrom(output=astrompath, input=astrompath)


# %% angle astrometry


def astrom_angle(astrompath, parentdir, chip, rot_val=48):
    os.chdir(gen_pipeline_file_name())
    ramppath = os.path.join(parentdir,  'C%i' % chip)
    print('running initial astrometry on ramp imgs...')

    if not rot_val:
        placeholder = 48
    else:
        placeholder = rot_val

    print(
        '\nEquivalent argparse cmd: python ./preprocess/astromangle_new.py -input %s -output %s -rot_val %s' % (ramppath,
                                                                                                              astrompath,
                                                                                                              placeholder))

    astromangle_new.astrom_angle(input_dir=ramppath, output_dir=astrompath, rot_val=rot_val)


# %% flat fielding


def flatfielding(astrompath, FFpath, band, chip):
    os.chdir(gen_pipeline_file_name())
    print('using master flat to flat field ramp imgs..')
    flatpath = gen_mflat_file_name(band, chip)

    print('\nEquivalent argparse cmd: python ./preprocess/flatfield.py -in_path %s -out_path %s'
          ' -flat_path %s' % (astrompath, FFpath, flatpath))

    flatfield.flat_field_cmd(in_path=astrompath, out_path=FFpath, flat_path=flatpath)


# %% sky gen


def sky(astrompath, skypath, sigma, chip):
    os.chdir(gen_pipeline_file_name())
    filelist = [f for f in os.listdir(skypath) if f.endswith('.C{}.fits'.format(chip))]
    if filelist:
        print('Previous sky %s found! Skipping sky gen..\n' % filelist[0])
        pass
    else:
        print('generating sky...')
        FFstring = '_FF'
        if FFstring in astrompath:
            no_flat = False
        else:
            no_flat = True

        print('\nEquivalent argparse cmd: python ./sky/gen_sky.py -in_path %s -sky_path %s -sigma %s ' % (astrompath, skypath, sigma))

        gen_sky.sky_gen(in_path=astrompath, sky_path=skypath, sigma=sigma, no_flat=no_flat)


# %% sky sub


def skysub(astrompath, subpath, skypath, chip, sky_override_path=None):
    os.chdir(gen_pipeline_file_name())
    if sky_override_path:
        skyfilepath = sky_override_path
    else:
        for file in os.listdir(skypath):
            if file.endswith('.C{}.fits'.format(chip)):
                skyfile = file
                skyfilepath = os.path.join(skypath, skyfile)
    print('cropping and subtracting sky...')
    FFstring = '_FF'
    if FFstring in astrompath:
        no_flat = False
    else:
        no_flat = True

    print('\nEquivalent argparse cmd: python ./sky/sky_sub.py -in_path %s -out_path %s -sky_path %s' %
          (astrompath, subpath, skyfilepath))

    sky_sub.sky_sub(in_path=astrompath, out_path=subpath, sky_path=skyfilepath, no_flat=no_flat)

# %%sex sky sub


def sexskysub(astrompath, subpath):
    # currently w/o ff?
    os.chdir(gen_pipeline_file_name())

    print('Using sextractor sky sub...')
    sky.sky_sub(in_path=astrompath, out_path=subpath, sex=True)


# %% astrometry shift


def shift(subpath, band):
    os.chdir(gen_pipeline_file_name())
    print('Shifting astrometry...')
    all_fits = [f for f in sorted(os.listdir(subpath)) if f.endswith('.flat.fits')]
    if len(all_fits) >= 100:
        imgname = all_fits[0]

        print('\nEquivalent argparse cmd: python ./astrom/astrom_shift.py -dir %s -imagename %s -band %s' %
              (subpath, imgname, band))

        astrom_shift.shift(directory=subpath, imagename=imgname, band=band)
    else:
        imgname = all_fits[0]

        print('\nEquivalent argparse cmd: python ./astrom/astrom_shift.py -dir %s -imagename %s -band %s' %
              (subpath, imgname, band))

        astrom_shift.shift(directory=subpath, imagename=imgname, band=band)


# %% better astrometry


def astromatic_astrometry(subpath, sex=None):
    os.chdir(gen_pipeline_file_name())
    print('using SXTRCTR and SCAMP to generate better astrometry...')
    # if sex:
    #     try:
    #        command = 'python ./astrom/astrometry.py -scamp -path %s' % (subpath)
    #        print('Executing command: %s' % command)
    #        rval = subprocess.run(command.split(), check=True)
    #    except subprocess.CalledProcessError as err:
    #        print('Could not run with exit error %s' % err)
    # else:

    print('\nEquivalent argparse cmd: python ./astrom/astrometry.py -path %s' % subpath)

    astrometry.astrometry(path=subpath)

# %% stacking


def stacking(subpath, stackpath, chip):
    os.chdir(gen_pipeline_file_name())
    print('stacking all images using SWARP...')

    print('\nEquivalent argparse cmd: python ./stack/stack.py -sub %s -stack %s -chip %i' % (subpath, stackpath, chip))

    stack.stack(subpath=subpath, stackpath=stackpath, chip=chip)


# %% packing compression


def fpack(stackpath, chip):
    os.chdir(gen_pipeline_file_name())
    print('compressing stacked images with fpack...')
    instack = sorted(os.listdir(stackpath))
    stackimg = []
    for f in instack:
        if fnmatch.fnmatch(f, 'coadd.Open-*.C%i.fits' % chip):
            stackimg.append(f)
    for f in stackimg:
        try:
            command = 'fpack -D -Y %s%s' % (stackpath, f)
            subprocess.run(command.split(), check=True)
        except subprocess.CalledProcessError as err:
            print('Could not run with exit error %s' % err)
    print('Stack compressed!')


# %% astromnet refining


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


# %%


defaults = dict(sigma=4)


def master(
        parentdir, chip, band, sigma=4, rot_val=None, no_ff=False, no_shift=False, sex=False, compress=False,
        net_refine=False, sky_override=None
):
    if no_ff:
        astromdir, skydir, subdir, stackdir = makedirectories(parentdir, chip)
        astrom_angle(astromdir, parentdir, chip, rot_val)
        sky(astromdir, skydir, sigma, chip)
        if sex:
            sexskysub(astromdir, subdir)
        else:
            skysub(astromdir, subdir, skydir, chip, sky_override)
        astromatic_astrometry(subdir)
        stacking(subdir, stackdir, chip)
    else:
        astromdir, skydir, subdir, stackdir = makedirectories(parentdir, chip)
        FFdir = makedirectoriesFF(parentdir, chip)
        astrom_angle(astromdir, parentdir, chip, rot_val)
        flatfielding(astromdir, FFdir, band, chip)
        if sex or sky_override:
            pass
        else:
            sky(FFdir, skydir, sigma, chip)
        if sex:
            sexskysub(FFdir, subdir)
        else:
            skysub(FFdir, subdir, skydir, chip, sky_override)
        if net_refine:
            astromnet_refine(subdir)
        else:
            if not no_shift:
                shift(subdir, band)
            else:
                pass
        astromatic_astrometry(subdir)
        stacking(subdir, stackdir, chip)
        if compress:
            fpack(stackdir, chip)


def main():
    parser = argparse.ArgumentParser(
        description='Automation of the backbone of pipeline, currently processes 1 chip at a time')
    parser.add_argument('-parent', type=str,
                        help='[str], parent directory of outputs, should include folders of the chips ramp data')
    parser.add_argument('-chip', type=int, help='[int], number of detector')
    parser.add_argument('-band', type=str, help='*NOT NECESSARY UNLESS USING -FF* [str], band of images, ex. "J"',
                        default=None)
    parser.add_argument('-sigma', type=int, help='[int], sigma value for sky sub sigma clipping, default = 4',
                        default=defaults["sigma"])
    parser.add_argument('-rot_val', type=float, help='[float] optional, put in your rot angle in deg,'
                                                     ' if you had a non-default rotation angle in your obs'
                                                     ' (default = 48 deg or 172800")', default=None)
    parser.add_argument('-no_FF', action='store_true', help='optional flag, does not use flat fielding in pipeline')
    parser.add_argument('-sex', action='store_true',
                        help='optional flag, to utlize sextractor background subtraction instead, do not currently use!')
    parser.add_argument('-compress', action='store_true', help='optional flag, use fpack to compress stacked images')
    parser.add_argument('-no_shift', action='store_true',
                        help='optional flag, STOPS use of astrometric shifting script')
    # parser.add_argument('-skygen_start',  action='store_true', help='optional flag, starts pipeline at sky gen step')
    # parser.add_argument('-skysub_start', action='store_true', help='optional flag, starts pipeline at sky sub step')
    # parser.add_argument('-astrom_start', action='store_true', help='optional flag, starts pipeline at sxtrctr / scamp step')
    # parser.add_argument('-stack_start', action='store_true',help='optional flag, starts pipeline at swarp step')
    parser.add_argument('-net_refine', action='store_true',
                        help='optional flag, used to automatically refine astrometry using '
                             'astrometry.net')
    parser.add_argument('-sky_override', type=str, help='[str], Optional path to specify already generated '
                                                        'sky to use in sky sub, skipping sky gen. Input full file path.',
                        default=None)
    args, unknown = parser.parse_known_args()

    master(args.parent, args.chip, args.band, args.sigma, args.rot_val, args.no_FF, args.no_shift, args.sex,
           args.compress, args.net_refine, args.sky_override)


if __name__ == "__main__":
    main()
