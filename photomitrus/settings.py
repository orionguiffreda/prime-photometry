"""
Settings for pipeline
    - Directory locations
    - Detector specific settings
    - Filter specific settings
    - Catalog settings
"""
# %% Config File Names
import os

# base_dir = os.path.dirname(__file__)
# base_dir = os.path.dirname(os.path.abspath('__file__'))

PIPELINE_DEFAULT_DIR = '/mnt/d/PRIME_photometry_test_files/'

PHOTOMETRY_MAG_LOWER_LIMIT = 12.5
PHOTOMETRY_MAG_UPPER_LIMIT = 21
PHOTOMETRY_QUERY_WIDTH = 48  # in arcmin, should cover the whole chip even with 45 deg rotation
PHOTOMETRY_QUERY_CATALOGS = {'VHS': ['J', 'II/367/'], 'VIKING': ['J', 'II/343/viking2'], '2MASS': ['J', 'II/246/'],
'Skymapper': ['Z', 'II/379/smssdr4'],
                             'DES_Z': ['Z', 'II/371/des_dr2'], 'SDSS': ['Z', 'V/154/sdss16'], 'UKIDSS': ['Y', 'II/319/las9'],
                             'DES_Y': ['Y', 'II/371/des_dr2']
                             }


def gen_config_file_name(filename):
    base_dir = os.path.dirname(__file__)  # os.path.abspath('__file__')
    return os.path.join(base_dir, 'configs', filename)


def gen_mask_file_name(filename):
    base_dir = os.path.dirname(__file__)  # os.path.abspath('__file__')
    return os.path.join(base_dir, 'weightmaps', filename)


def gen_pipeline_file_name():
    base_dir = os.path.dirname(os.path.realpath(__file__))
    return base_dir


def gen_master_name():
    base_dir = os.path.dirname(__file__)
    return os.path.join(base_dir, 'master.py')


def gen_mflat_file_name(band, chip):
    base_dir = os.path.dirname(os.path.realpath(__file__))
    flat_dir = os.path.join(base_dir, 'mflats')
    mflat_list = [
        f for f in sorted(os.listdir(flat_dir)) if f.endswith('.fits') if '.%s.' % band in f if 'C%s' % chip in f]

    # get latest mflat
    mflat_list = sorted(mflat_list, reverse=True)
    filename = mflat_list[0]
    return os.path.join(flat_dir, filename)


# %%
os.chdir('/mnt/c/PycharmProjects/prime-photometry/photomitrus/')
test = gen_mflat_file_name('J', 1)

# %%
import pandas as pd

# Settings for list of directory+filenames
object = 'field4057'
filter = 'H'
chip = 4


def flist(Object=object, Filter=filter, Chip=chip):
    log = pd.read_csv('/mnt/d/PRIME_photometry_test_files/ramp_fit_log_2023-12-10.clean.dat',
                      delimiter=' ')  # reads in csv
    if Filter == 'Z':
        fnames = log['filename'][
            log['CHIP'] == Chip & log['OBJNAME'].str.contains(str(Object)) & ~log['OBJNAME'].str.contains('test') & log[
                'FILTER1'].str.contains('Z')
            & log['Open'].str.contains(str(Filter))]
        fnames = fnames.tolist()
        dir = ('/mnt/d/PRIME_photometry_test_files/C{}/'.format(Chip))
        fullnames = [dir + x for x in sorted(fnames)]
    else:
        # screening log w/ constraints
        fnames = log['filename'][log['OBJNAME'].str.contains(str(Object)) & ~log['OBJNAME'].str.contains('test') & log[
            'FILTER1'].str.contains('Open')
                                 & log['FILTER2'].str.contains(str(Filter)) & log['OBSERVER'].str.contains('NASA')]
        fnames = fnames.tolist()
        # adding path
        dir = ('/mnt/d/PRIME_photometry_test_files/C{}/'.format(Chip))
        fullnames = [dir + x for x in sorted(fnames)]
        # fullnames = fnames
    fullnames = [f.replace('fits.ramp', 'ramp.fits') for f in fullnames]
    fullnames = [f.replace('C1', 'C%s' % chip) for f in fullnames]
    if not fullnames:
        print('No files found for specified fields!')
    return fullnames


# a = flist()
# %% make directories
from pathlib import Path


def makedirs(dir, chip):
    os.chdir(dir)
    sky = os.path.join(dir, 'sky')
    stack = os.path.join(dir, 'stack')
    a = os.path.join(dir, 'C%i_astrom' % chip)
    sub = os.path.join(dir, 'C%i_sub' % chip)
    directories = (a, sky, sub, stack)
    for directory in directories:
        if os.path.exists(directory):
            print(directory + ' already exists!')
        else:
            os.mkdir(directory)
            print(directory)
    # skyexists = os.path.exists('sky')
    # if not skyexists:
    #     os.mkdir(sky)
    #     print(sky)
    # if skyexists:
    #     print(sky + ' already exists!')
    # stackexists = os.path.exists(stack)
    # if not stackexists:
    #     os.mkdir(stack)
    #     print(stack)
    # if stackexists:
    #     print(dir + 'stack already exists!')
    # aexists = os.path.exists(a)
    # if not aexists:
    #     os.mkdir(a)
    #     print(dir + a)
    # if aexists:
    #     print(dir + a+' already exists!')
    # subexists = os.path.exists(sub)
    # if not subexists:
    #     os.mkdir(sub)
    #     print(dir + sub)
    # if subexists:
    #     print(dir + sub+' already exists!')
    # dirnames = os.listdir('.')
    # astromlist = [i for i in dirnames if i.endswith(a)]
    # astrom = ' '.join(astromlist)
    # sublist = [i for i in dirnames if i.endswith(sub)]
    # sub = ' '.join(sublist)
    # skylist = [i for i in dirnames if i.endswith('sky')]
    # sky = ' '.join(skylist)
    # stacklist = [i for i in dirnames if i.endswith('stack')]
    # stack = ' '.join(stacklist)
    # return dir + astrom +'/', dir + sky +'/', dir + sub +'/', dir + stack +'/'
    return directories


# flat field directory creation
def makedirsFF(dir, chip):
    # os.chdir(dir)
    FF = os.path.join(dir, 'C%i_FF' % chip)
    skyexists = os.path.exists(FF)
    if not skyexists:
        os.mkdir(FF)
        print(FF)
    if skyexists:
        print(FF + ' already exists!')
    # dirnames = os.listdir('.')
    # FFdir = [i for i in dirnames if i.endswith(FF)]
    # FFname = ' '.join(FFdir)
    return FF


# %%
# astrom, sky, sub, stack = makedirs('/mnt/d/PRIME_photometry_test_files/GRB240205B/J_Band/flats/',1)
# FF = makedirsFF('/mnt/d/PRIME_photometry_test_files/xrf_3_20240318/J_Band/',2)

# %% for multiple object fields
"""
object = ['field4037']
filter = 'J'
chip = 1
def flist(Object=object, Filter=filter, Chip=chip):
    log = pd.read_csv('/mnt/d/PRIME_photometry_test_files/ramp_fit_log_2023-12-10.clean.dat', delimiter=' ') #reads in csv
    if Filter == 'Z':
        fnames = []
        for f in Object:
            fobjnames = log['filename'][log['CHIP'] == Chip & log['OBJNAME'].str.contains(str(Object)) & log['FILTER1'].str.contains('Z')
            & log['Open'].str.contains(str(Filter))]
            fobjnames = fobjnames.tolist()
            fnames.extend(fobjnames)
        dir = ('/mnt/d/PRIME_photometry_test_files/C{}/'.format(Chip))
        fullnames = [dir + x for x in sorted(fnames)]
    else:
        #screening log w/ constraints
        fnames = []
        for f in Object:
            fobjnames = log['filename'][log['CHIP'] == Chip & log['OBJNAME'].str.contains(str(f)) & log['FILTER1'].str.contains('Open')
                 & log['FILTER2'].str.contains(str(Filter))]
            fobjnames = fobjnames.tolist()
            fnames.extend(fobjnames)
            #adding path
        dir = ('/mnt/d/PRIME_photometry_test_files/C{}/'.format(Chip))
        fullnames = [dir + x for x in sorted(fnames)]
    if not fullnames:
        print('No files found for specified fields!')
    return fullnames

a = flist()
"""
