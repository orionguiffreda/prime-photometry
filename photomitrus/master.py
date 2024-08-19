import os
import sys
import subprocess
import argparse
import fnmatch
#sys.path.insert(0,'/mnt/c/PycharmProjects/prime-photometry/photomitrus/')
from settings import makedirs
from settings import makedirsFF
from settings import gen_pipeline_file_name
from settings import gen_mflat_file_name


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
            rval = subprocess.run(command.split(), check=True)
        except subprocess.CalledProcessError as err:
            print('Could not run with exit error %s'%err)
    if not chip:
        print('running astrometry.net on subbed imgs...')
        try:
            command = 'python ./preprocess/gen_astrometry.py -hard -output %s -input %s -rad 1 -ds 3' % (astrompath, astrompath)
            print('Executing command: %s' % command)
            rval = subprocess.run(command.split(), check=True)
        except subprocess.CalledProcessError as err:
            print('Could not run with exit error %s'%err)

#%% angle astrometry
def astrom_angle(astrompath,parentdir, chip):
    os.chdir(gen_pipeline_file_name())
    ramppath = parentdir + 'C%i/' % (chip)
    print('running initial astrometry on ramp imgs...')
    try:
        command = 'python ./preprocess/astromangle_new.py -input %s -output %s' % (ramppath, astrompath)
        print('Executing command: %s' % command)
        rval = subprocess.run(command.split(), check=True)
    except subprocess.CalledProcessError as err:
        print('Could not run with exit error %s'%err)

#%% flat fielding
def flatfield(astrompath,FFpath,filter,chip):
    os.chdir(gen_pipeline_file_name())
    print('using master flat to flat field ramp imgs..')
    flatpath = gen_mflat_file_name(filter,chip)
    try:
        command = 'python ./preprocess/flatfield.py -in_path %s -out_path %s -flat_path %s' % (astrompath,FFpath,flatpath)
        print('Executing command: %s' % command)
        rval = subprocess.run(command.split(), check=True)
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
            rval = subprocess.run(command.split(), check=True)
        except subprocess.CalledProcessError as err:
            print('Could not run with exit error %s' % err)
    else:
        try:
            command = 'python ./sky/gen_sky.py -in_path %s -sky_path %s -sigma %i' % (astrompath, skypath, sigma)
            print('Executing command: %s' % command)
            rval = subprocess.run(command.split(), check=True)
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
            rval = subprocess.run(command.split(), check=True)
        except subprocess.CalledProcessError as err:
            print('Could not run with exit error %s' % err)
    else:
        try:
            command = 'python ./sky/sky.py -in_path %s -out_path %s -sky_path %s%s' % (astrompath, subpath, skypath, skyfile)
            print('Executing command: %s' % command)
            rval = subprocess.run(command.split(), check=True)
        except subprocess.CalledProcessError as err:
            print('Could not run with exit error %s' % err)
#%%sex sky sub
def sexskysub(astrompath,subpath):
    os.chdir(gen_pipeline_file_name())
    try:
        command = 'python ./sky/sky.py -sex -in_path %s -out_path %s' % (astrompath, subpath)
        print('Executing command: %s' % command)
        rval = subprocess.run(command.split(), check=True)
    except subprocess.CalledProcessError as err:
        print('Could not run with exit error %s' % err)

#%% astrometry shift
def shift(subpath, filter):
    os.chdir(gen_pipeline_file_name())
    print('Shifting astrometry...')
    all_fits = [f for f in sorted(os.listdir(subpath)) if f.endswith('.flat.fits')]
    imgname = all_fits[0]
    try:
        command = 'python ./astrom/astrom_shift.py -remove -pipeline -dir %s -imagename %s -filter %s' % (subpath, imgname, filter)
        print('Executing command: %s' % command)
        rval = subprocess.run(command.split(), check=True)
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
        rval = subprocess.run(command.split(), check=True)
    except subprocess.CalledProcessError as err:
        print('Could not run with exit error %s' % err)

#%% stacking
def stack(subpath, stackpath,chip):
    os.chdir(gen_pipeline_file_name())
    print('stacking all images using SWARP...')
    try:
        command = 'python ./stack/stack.py -sub %s -stack %s -chip %i' % (subpath, stackpath, chip)
        print('Executing command: %s' % command)
        rval = subprocess.run(command.split(), check=True)
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
            #print('Executing command: %s' % command)
            rval = subprocess.run(command.split(), check=True)
        except subprocess.CalledProcessError as err:
            print('Could not run with exit error %s' % err)
    print('Stack compressed!')

#%%
defaults = dict(sigma=4)

#%%
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Automation of the backbone of pipeline, currently processes 1 chip at a time')
    parser.add_argument('-FF', action='store_true', help='optional flag, changes the default pipeline to include a flat fielding step, should improve photometry')
    parser.add_argument('-angle', action='store_true', help='optional flag, forgoes astrometry.net, utilizes center and corner positions to create wcs')
    parser.add_argument('-sex', action='store_true', help='optional flag, to utlize sextractor background subtraction instead, do not currently use!')
    parser.add_argument('-fpack', action='store_true',help='optional flag, use fpack to compress stacked images')
    parser.add_argument('-skygen_start',  action='store_true', help='optional flag, starts pipeline at sky gen step')
    parser.add_argument('-skysub_start', action='store_true', help='optional flag, starts pipeline at sky sub step')
    parser.add_argument('-astrom_start', action='store_true', help='optional flag, starts pipeline at sxtrctr / scamp step')
    parser.add_argument('-stack_start', action='store_true',help='optional flag, starts pipeline at swarp step')
    parser.add_argument('-qual_check', action='store_true', help='optional flag, puts an input after stacking, if the stacked img is '
                                                                 'subpar, will rerun astrometry.net on proc. images and restack')
    parser.add_argument('-qual_check_only', action='store_true', help='optional flag, include both this and -qual_check if'
                                                                      ' you want to quality check after the pipeline is run already')
    parser.add_argument('-refine', action='store_true', help='optional flag, used to automatically refine astrometry using '
                                                             'the same method as "-qual_check"')
    parser.add_argument('-parent', type=str, help='[str], parent directory of outputs, should include folders of the chips ramp data')
    parser.add_argument('-chip', type=int, help='[int], number of detector')
    parser.add_argument('-filter', type=str, help='*NOT NECESSARY UNLESS USING -FF* [str], filter of images, ex. "J"',default=None)
    parser.add_argument('-sigma', type=int, help='[int], sigma value for sky sub sigma clipping, default = 4', default=defaults["sigma"])
    args = parser.parse_args()

    if not args.qual_check_only:
        if args.FF:
            if args.skygen_start:
                astromdir, skydir, subdir, stackdir = makedirectories(args.parent, args.chip)
                FFdir = makedirectoriesFF(args.parent, args.chip)
                #flatfield(astromdir, FFdir, args.filter, args.chip)
                sky(FFdir, skydir, args.sigma)
                skysub(FFdir, subdir, skydir, args.chip)
                if args.refine:
                    initastrom(subdir, subdir)
                    oldlist = [f for f in sorted(os.listdir(subdir)) if f.endswith('.flat.fits')]
                    newlist = [j for j in sorted(os.listdir(subdir)) if j.endswith('.flat.new')]
                    if len(oldlist) <= 15:
                        errnum = 3
                    else:
                        errnum = 6
                    if len(newlist) < len(oldlist)-errnum:
                        print('Not enough success with new astrometry! Continuing with initial astrometry...')
                        for f in newlist:
                            os.remove(subdir + f)
                    else:
                        for f in oldlist:
                            os.remove(subdir + f)
                    print('%i files refined! removed old fits files' % len(newlist))
                    astrometry(subdir)
                    stack(subdir, stackdir,args.chip)
                    if args.fpack:
                        fpack(stackdir,args.chip)
                else:
                    astrometry(subdir)
                    stack(subdir, stackdir,args.chip)
                    if args.fpack:
                        fpack(stackdir, args.chip)
            elif args.skysub_start:
                astromdir, skydir, subdir, stackdir = makedirectories(args.parent, args.chip)
                FFdir = makedirectoriesFF(args.parent, args.chip)
                if args.sex:
                    sexskysub(FFdir,subdir)
                    astrometry(subdir, args.sex)
                else:
                    skysub(FFdir, subdir, skydir, args.chip)
                    astrometry(subdir)
                stack(subdir, stackdir,args.chip)
            elif args.astrom_start:
                astromdir, skydir, subdir, stackdir = makedirectories(args.parent, args.chip)
                FFdir = makedirectoriesFF(args.parent, args.chip)
                astrometry(subdir,args.sex)
                stack(subdir, stackdir,args.chip)
            elif args.stack_start:
                astromdir, skydir, subdir, stackdir = makedirectories(args.parent, args.chip)
                FFdir = makedirectoriesFF(args.parent, args.chip)
                stack(subdir, stackdir,args.chip)
                if args.fpack:
                    fpack(stackdir, args.chip)
            else:
                astromdir, skydir, subdir, stackdir = makedirectories(args.parent, args.chip)
                FFdir = makedirectoriesFF(args.parent, args.chip)
                if args.angle:
                    astrom_angle(astromdir, args.parent, args.chip)
                else:
                    initastrom(astromdir, args.parent, args.chip)
                flatfield(astromdir, FFdir, args.filter, args.chip)
                if args.sex:
                    pass
                else:
                    sky(FFdir, skydir, args.sigma)
                if args.sex:
                    sexskysub(FFdir,subdir)
                else:
                    skysub(FFdir, subdir, skydir, args.chip)
                if args.refine:
                    initastrom(subdir, subdir)
                    oldlist = [f for f in sorted(os.listdir(subdir)) if f.endswith('.flat.fits')]
                    newlist = [j for j in sorted(os.listdir(subdir)) if j.endswith('.flat.new')]
                    if len(oldlist) <= 15:
                        errnum = 3
                    else:
                        errnum = 5
                    if len(newlist) < len(oldlist)-errnum:
                        print('Not enough success with new astrometry! (%i fields) Continuing with initial astrometry...' % len(newlist))
                        for f in newlist:
                            os.remove(subdir + f)
                    else:
                        for f in oldlist:
                            os.remove(subdir + f)
                        print('%i files refined! removed old fits files' % len(newlist))
                    astrometry(subdir)
                    stack(subdir, stackdir,args.chip)
                    if args.fpack:
                        fpack(stackdir,args.chip)
                else:
                    shift(subdir,args.filter)
                    astrometry(subdir,args.sex)
                    stack(subdir, stackdir,args.chip)
                    if args.fpack:
                        fpack(stackdir,args.chip)
        if not args.FF:
            if args.skygen_start:
                astromdir, skydir, subdir, stackdir = makedirectories(args.parent, args.chip)
                sky(astromdir, skydir, args.sigma)
                skysub(astromdir, subdir, skydir, args.chip)
                astrometry(subdir)
                stack(subdir, stackdir,args.chip)
            elif args.skysub_start:
                astromdir, skydir, subdir, stackdir = makedirectories(args.parent, args.chip)
                skysub(astromdir, subdir, skydir, args.chip)
                astrometry(subdir)
                stack(subdir, stackdir,args.chip)
            elif args.astrom_start:
                astromdir, skydir, subdir, stackdir = makedirectories(args.parent, args.chip)
                astrometry(subdir)
                stack(subdir, stackdir,args.chip)
            elif args.stack_start:
                astromdir, skydir, subdir, stackdir = makedirectories(args.parent, args.chip)
                stack(subdir, stackdir,args.chip)
            else:
                astromdir, skydir, subdir, stackdir = makedirectories(args.parent, args.chip)
                if args.angle:
                    astrom_angle(astromdir, args.parent, args.chip)
                else:
                    initastrom(astromdir, args.parent, args.chip)
                sky(astromdir, skydir, args.sigma)
                if args.sex:
                    sexskysub(astromdir,subdir)
                else:
                    skysub(astromdir, subdir, skydir, args.chip)
                astrometry(subdir)
                stack(subdir, stackdir,args.chip)

        if args.qual_check:
            quality = input("Is stacked image good quality? (Input 'Y' or 'N'): ")
            if quality == 'Y':
                print('Pipeline complete! Enjoy your stacked image!')
            elif quality == 'N':
                print('Understood, rerunning astrometry.net on processed images!')
                for f in sorted(os.listdir(subdir)):
                    if f.endswith('.head') or f.endswith('.cat'):
                        os.remove(subdir+f)
                print('removed old scamp and sex files')
                initastrom(subdir, subdir)
                oldlist = [f for f in sorted(os.listdir(subdir)) if f.endswith('.flat.fits')]
                newlist = [j for j in sorted(os.listdir(subdir)) if j.endswith('.flat.new')]
                if len(newlist) < len(oldlist) - 5:
                    print('Not enough success with new astrometry, %i files refined.. Breaking...' % len(newlist))
                    sys.exit()
                else:
                    for f in oldlist:
                        os.remove(subdir + f)
                print('%i files refined! removed old fits files' % len(newlist))
                astrometry(subdir)
                stack(subdir, stackdir, args.chip)

    else:
        if args.qual_check:
            if args.qual_check_only:
                astromdir, skydir, subdir, stackdir = makedirectories(args.parent, args.chip)
                FFdir = makedirectoriesFF(args.parent, args.chip)
            quality = input("Is stacked image good quality? (Input 'Y' or 'N'): ")
            if quality == 'Y':
                print('Pipeline complete! Enjoy your stacked image!')
            elif quality == 'N':
                print('Understood, rerunning astrometry.net on processed images!')
                for f in sorted(os.listdir(subdir)):
                    if f.endswith('.head') or f.endswith('.cat'):
                        os.remove(subdir+f)
                print('removed old scamp and sex files')
                initastrom(subdir, subdir)
                oldlist = [f for f in sorted(os.listdir(subdir)) if f.endswith('.flat.fits')]
                newlist = [j for j in sorted(os.listdir(subdir)) if j.endswith('.flat.new')]
                if len(newlist) < len(oldlist) - 5:
                    print('Not enough success with new astrometry, %i files refined.. Breaking...' % len(newlist))
                    sys.exit()
                else:
                    for f in oldlist:
                        os.remove(subdir + f)
                print('%i files refined! removed old fits files' % len(newlist))
                astrometry(subdir)
                stack(subdir, stackdir, args.chip)
