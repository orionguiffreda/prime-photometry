"""
Settings for pipeline
    - Directory locations
    - Detector specific settings
    - Filter specific settings
    - Catalog settings
"""
# %% Config File Names
import os
import pandas as pd
from datetime import datetime

# base_dir = os.path.dirname(__file__)
# base_dir = os.path.dirname(os.path.abspath('__file__'))

PIPELINE_DEFAULT_DIR = '/mnt/photometry/'

PHOTOMETRY_MAG_LOWER_LIMIT = 12.5
PHOTOMETRY_MAG_UPPER_LIMIT = 23
PHOTOMETRY_QUERY_WIDTH = 48  # in arcmin, should cover the whole chip even with 45 deg rotation
PHOTOMETRY_QUERY_CATALOGS = {'VHS': ['J', 'II/367/'], 'VIKING': ['J', 'II/343/viking2'], 'VVV_J': ['J', 'II/348/vvv2'],
                             '2MASS': ['J', 'II/246/'],
                            'VVV_Y': ['J', 'II/348/vvv2'], 'Skymapper': ['Z', 'II/379/smssdr4'],'SDSS': ['Z', 'V/154/sdss16'],'DES_Z': ['Z', 'II/371/des_dr2'],
                             'VVV_Z': ['J', 'II/348/vvv2'], 'UKIDSS': ['Y', 'II/319/las9'],'PanSTARRS': ['Y', 'II/349/ps1'],'DES_Y': ['Y', 'II/371/des_dr2']
                             }
# AB mag surveys: DES, Skymapper, SDSS, & PSTARRS
AB_OFFSET_DICT = {'J': 0.94, 'H': 1.38, 'Y': 0.62, 'Z': 0.52}

PHOTOMETRY_LIM_MAGS = {
    'VHS': 19.5, 'VIKING': 21.2, 'VVV': 20.2, '2MASS': 15.5, 'Skymapper': 22, 'DES_Z': 23.1, 'SDSS': 23,
                       'UKIDSS': 20.2, 'DES_Y': 21.7, 'PanSTARRS': 21.4
                       }

# chip avg zp (for bulge field astrom)
# 23.88
# 21.836
GB_QUERY_CATALOGS = {'VVV': ['J', 'II/348/vvv2'], '2MASS': ['J', 'II/246/'],
                     'Skymapper': ['Z', 'II/379/smssdr4']}

CHIP_ZPS = {'Z': [23.88], 'Y': [23.88], 'J': [23.88], 'H': [24.111]}
# 24.111

def bulge_checker(case):
    bulge_list = ['bulge', 'gb', 'gp', 'plane']

    bulge_check = [name for name in bulge_list if name in case.lower()]
    if bulge_check:
        return True
    else:
        return False


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


def gen_mflat_file_name(band, chip, date=None):
    base_dir = os.path.dirname(os.path.realpath(__file__))
    flat_dir = os.path.join(base_dir, 'mflats')
    mflat_list = [
        f for f in sorted(os.listdir(flat_dir)) if f.endswith('.fits') if '.%s.' % band in f if 'C%s' % chip in f
        if len(f) <= 24]

    if date:
        # get mflat closest to obs date, if there is a tie, it picks the earlier one to be safe
        def closest_file(file_list, target_date):
            target = datetime.strptime(target_date, "%Y%m%d")

            def extract_date(file):
                return datetime.strptime(file.split(".")[2], "%Y%m%d")

            return min(file_list, key=lambda f: (abs((extract_date(f) - target).days), extract_date(f)))

        filename = closest_file(mflat_list, date)
        print('Getting mflat closest to given date: ', filename)
    else:
        # if no date given, get latest mflat
        mflat_list = sorted(mflat_list, reverse=True)
        filename = mflat_list[0]
        print('No date, getting latest mflat: ', filename)
    return os.path.join(flat_dir, filename)


def mflat_checker(date):
    base_dir = os.path.dirname(os.path.realpath(__file__))
    flat_dir = os.path.join(base_dir, 'mflats')
    # flat_dir = '/home/alex/PycharmProjects/prime-photometry/photomitrus/mflats/'
    mflat_list = [
        f for f in sorted(os.listdir(flat_dir)) if f.endswith('.fits') if '.%s.' % date in f]
    if not mflat_list:
        print('\nNo master flat currently generated for this date!')
        return False
    else:
        print('\nMaster flats exist for this date!')
    return True


#%%

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

# directory creation


def makedirs(dir, chip):
    os.chdir(dir)
    sky = os.path.join(dir, 'sky')
    stack = os.path.join(dir, 'stack')
    a = os.path.join(dir, 'C%i_astrom' % chip)
    FF = os.path.join(dir, 'C%i_FF' % chip)
    sub = os.path.join(dir, 'C%i_sub' % chip)
    directories = (a, FF, sky, sub, stack)
    for directory in directories:
        if os.path.exists(directory):
            print(directory + ' already exists!')
        else:
            os.mkdir(directory)
            print(directory)
    return directories


# flat field directory creation, probably outdated


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
