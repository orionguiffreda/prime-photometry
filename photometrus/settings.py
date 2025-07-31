"""
Settings for pipeline
    - Directory locations
    - Detector specific settings
    - Filter specific settings
    - Catalog settings
"""
# %% Config File Names
import os
from pathlib import Path

import pandas as pd
from datetime import datetime
from astropy.io import fits
import json5


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
GET_DATA_SETTINGS = dict(
    replace_list = [
        ('all sky grid', 'all_sky_grid'),
        ('no_grid  ', 'no_grid '),
        ('no_grid_test  ', 'no_grid_test '),
        ('field10554_ test', 'field10554_test'),
        ('fucus test', 'focus_test'),
        ('galactic plane', 'galactic_plane'),
        ('standard star', 'standard_star'),
        ('rotator test', 'rotator_test'),
        ('SKY M7', 'SKY_M7'),
        ('LVC S230522n', 'LVC_S230522n'),
        ('Dec test', 'Dec_test'),
        ('Galctic plane', 'Galactic_plane'),
        ('Glactic plane', 'Galactic_plane'),
        ('Galactic plane', 'Galactic_plane'),
        ('Galactic Plane', 'Galactic_plane'),
        ('Vignetting Test', 'Vignetting_Test'),
        ('SPIS J0539-0059', 'SPIS_J0539-0059'),
        ('Frost Check', 'Frost_Check'),
        ('GainTest LED', 'GainTest_LED'),
    ]
)


def bulge_checker(case):
    bulge_list = ['bulge', 'gb', 'gp', 'plane']

    bulge_check = [name for name in bulge_list if name in case.lower()]
    if bulge_check:
        return True
    else:
        return False


def auto_bulge_detect(directory):
    fits_files = [f for f in os.listdir(directory) if f.endswith('.fits') or f.endswith('.new')]
    first_fits = os.path.join(directory, fits_files[0])
    hdr = fits.getheader(first_fits)

    case = hdr['OBJTYPE']
    bulge = bulge_checker(case)
    if bulge:
        print('Bulge field detected! Switching to bulge setup if not already specified!')

    return bulge


def gen_user_prime_dir():
    home_dir = Path.home()
    prime_dir = os.path.join(home_dir, '.prime')
    if not os.path.isdir(prime_dir):
        os.makedirs(prime_dir)
    return prime_dir


def gen_config_dir():
    return os.path.join(gen_user_prime_dir(), 'configs')


def gen_pipeline_file_name():
    base_dir = os.path.dirname(os.path.realpath(__file__))
    return base_dir


def gen_config_file_name(filename):
    # base_dir = os.path.dirname(__file__)  # os.path.abspath('__file__')
    # base_dir = gen_pipeline_file_name()
    return os.path.join(gen_config_dir(), filename)


def load_settings(settings_file='photometrus.json5'):
    settings_file_abs = os.path.join(gen_config_dir(), settings_file)
    if not os.path.isfile(settings_file_abs):
        settings_file_abs = os.path.join(gen_pipeline_file_name(), 'configs', settings_file)
    with open(settings_file_abs, 'r') as f:
        json_data = json5.load(f)
    return json_data


def update_settings(settings_file='photometrus.json5'):
    settings = load_settings(settings_file)
    global PIPELINE_DEFAULT_DIR
    global PHOTOMETRY_MAG_LOWER_LIMIT
    global PHOTOMETRY_MAG_UPPER_LIMIT
    global PHOTOMETRY_QUERY_WIDTH
    global PHOTOMETRY_QUERY_CATALOGS
    global AB_OFFSET_DICT
    global PHOTOMETRY_LIM_MAGS
    global GB_QUERY_CATALOGS
    global CHIP_ZPS
    global GET_DATA_SETTINGS

    PIPELINE_DEFAULT_DIR = settings['PIPELINE_DEFAULT_DIR']
    PHOTOMETRY_MAG_LOWER_LIMIT = settings['PHOTOMETRY_MAG_LOWER_LIMIT']
    PHOTOMETRY_MAG_UPPER_LIMIT = settings['PHOTOMETRY_MAG_UPPER_LIMIT']
    PHOTOMETRY_QUERY_WIDTH = settings['PHOTOMETRY_QUERY_WIDTH']
    PHOTOMETRY_QUERY_CATALOGS = settings['PHOTOMETRY_QUERY_CATALOGS']
    AB_OFFSET_DICT = settings['AB_OFFSET_DICT']
    PHOTOMETRY_LIM_MAGS = settings['PHOTOMETRY_LIM_MAGS']
    GB_QUERY_CATALOGS = settings['GB_QUERY_CATALOGS']
    CHIP_ZPS = settings['CHIP_ZPS']
    GET_DATA_SETTINGS = settings['GET_DATA_SETTINGS']


update_settings()


def gen_mask_file_name(filename):
    base_dir = os.path.dirname(__file__)  # os.path.abspath('__file__')
    return os.path.join(base_dir, 'weightmaps', filename)


def gen_master_name():
    base_dir = os.path.dirname(__file__)
    return os.path.join(base_dir, 'master.py')


def gen_flat_dir():
    base_dir = gen_user_prime_dir()
    flat_dir = os.path.join(base_dir, 'mflats')
    if not os.path.exists(flat_dir):
        os.makedirs(os.path.join(flat_dir, 'auxiliary_mflats'))
    return flat_dir


def gen_mflat_file_name(band, chip, date=None):
    flat_dir = gen_flat_dir()
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
    base_dir = gen_pipeline_file_name()
    flat_dir = gen_flat_dir()
    # flat_dir = '/home/alex/PycharmProjects/prime-photometry/photometrus/mflats/'
    mflat_list = [
        f for f in sorted(os.listdir(flat_dir)) if f.endswith('.fits') if '.%s.' % date in f]
    if not mflat_list:
        print('\nNo master flat currently generated for this date!')
        return False
    else:
        print('\nMaster flats exist for this date!')
    return True


#%%
"""
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
"""

# a = flist()
# %% make directories

# directory creation





# flat field directory creation, probably outdated

"""
def makedirsFF(basedir, chip):
    # os.chdir(basedir)
    FF = os.path.join(basedir, 'C%i_FF' % chip)
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
"""
