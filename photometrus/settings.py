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
from astroquery.vizier import Vizier, Conf
import warnings

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

ASTROM_QUERY_CATALOGS = {
    'VIRAC': ['J,H','II/364/virac2'],
    'GAIA': ['J,H,Y,Z','I/350/gaiaedr3']
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

WEIGHT_SIGMA_THRESHOLDS = {'C1': 2.0, 'C2': 2.0, 'C3': 2.0, 'C4': 1.0}

PRIME_FILTERS_DICT = {
    'Z': {
        'lambda': 890,
        'fwhm': 80
    },
    'Y': {
        'lambda': 1022.5,
        'fwhm': 105
    },
    'J': {
        'lambda': 1250,
        'fwhm': 160
    },
    'H': {
        'lambda': 1635,
        'fwhm': 290
    }
}

VIZIER_MIRRORS = [
    "vizier.cds.unistra.fr",
    "vizier.cfa.harvard.edu",
    "vizier.idia.ac.za",
    "vizier.nao.ac.jp",
    "vizier.iucaa.in",
    "vizier.inasan.ru",
    "vizier.china-vo.org"
]

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
    ],
    temperature_date_dict = {
        '2023-06-15T23:00:00': {
            '117': ['20240414', '20240904'],
            '122': ['20250122'],
            '111.5': ['20250220'],
            '114': ['20250320']
        }
    },
    ramp_cal_directory = '/nfs/home/prime/hamada/nlc_result/'
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
    clean = [f for f in fits_files if 'weight' not in f and 'PSF' not in f and 'GRB' not in f]
    if len(clean) == 0:
        print('No applicable fields to check, defaulting to non-bulge!')
        return False
    else:
        coadds = [f for f in clean if 'coadd' in f]
        if coadds:
            first_fits = os.path.join(directory, coadds[0])
        else:
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
    global WEIGHT_SIGMA_THRESHOLDS
    global PRIME_FILTERS_DICT
    global VIZIER_MIRRORS
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
    WEIGHT_SIGMA_THRESHOLDS = settings['WEIGHT_SIGMA_THRESHOLDS']
    PRIME_FILTERS_DICT = settings['PRIME_FILTERS_DICT']
    VIZIER_MIRRORS = settings['VIZIER_MIRRORS']
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


def gen_mflat_file_name(band, chip, date=None, sflat=False):
    flat_dir = gen_flat_dir()

    if sflat:
        name = 'sflat'
        length = 33
    else:
        name = 'mflat'
        length = 24

    mflat_list = [
        f for f in sorted(os.listdir(flat_dir)) if f.startswith(name) if f.endswith('.fits') if '.%s.' % band in f
        if 'C%s' % chip in f if len(f) <= length]

    if date:
        # get mflat closest to obs date, if there is a tie, it picks the earlier one to be safe
        def closest_file(file_list, target_date):
            target = datetime.strptime(target_date, "%Y%m%d")

            def extract_date(file):
                if sflat:
                    half = file.split("-")[0]
                    return datetime.strptime(half.split(".")[2], "%Y%m%d")
                else:
                    return datetime.strptime(file.split(".")[2], "%Y%m%d")

            return min(file_list, key=lambda f: (abs((extract_date(f) - target).days), extract_date(f)))

        filename = closest_file(mflat_list, date)
        print(f'Getting {name} closest to given date: ', filename)
    else:
        # if no date given, get latest mflat
        mflat_list = sorted(mflat_list, reverse=True)
        filename = mflat_list[0]
        print(f'No date, getting latest {name}: ', filename)
    return os.path.join(flat_dir, filename)


def mflat_checker(date, band=None, sflat=False):
    base_dir = gen_pipeline_file_name()
    flat_dir = gen_flat_dir()
    # flat_dir = '/home/alex/PycharmProjects/prime-photometry/photometrus/mflats/'
    if band:
        mflat_list = [
            f for f in sorted(os.listdir(flat_dir)) if f.endswith('.fits') if '.%s.%s.' % (band,date) in f]
    else:
        mflat_list = [
            f for f in sorted(os.listdir(flat_dir)) if f.endswith('.fits') if '.%s.' % date in f]
    if not mflat_list:
        if not sflat:
            print('\nNo master flat currently generated for this date!')
        return False
    else:
        if not sflat:
            print('\nMaster flats exist for this date!')
        return True


def gen_sflat_file_name(band, chip, date=None):
    flat_dir = gen_flat_dir()
    sflat_list = [
        f for f in sorted(os.listdir(flat_dir)) if f.endswith('.fits') if '.%s.' % band in f if 'C%s' % chip in f
        if f.startswith('sflat.')]

    sflat_list = sorted(sflat_list, reverse=True)
    filename = sflat_list[0]
    print('Getting latest sflat: ', filename)
    return os.path.join(flat_dir, filename)


def get_weight_thresh(chip):
    chosen_thresh = WEIGHT_SIGMA_THRESHOLDS['C%i' % chip]
    return chosen_thresh


# vizier mirror automated check
def set_vizier_mirror():
    """
    Checks vizier mirrors and auto picks working one, then updates astroquery's config.
    """

    mirrors = VIZIER_MIRRORS

    test_cat = 'II/246/'

    for url in mirrors:
        old = Conf.server
        try:
            Conf.server = url

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                Vizier(columns=["*"], row_limit=1).query_constraints(catalog=test_cat)

            print(f'Vizier mirror: {url}')
            break

        except Exception as e:
            print(f' {url} mirror check failed: {e}')
            Conf.server = old
    else:
        raise RuntimeError('No working Vizier mirror found!')