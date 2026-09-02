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
from astroquery.vizier import Vizier, conf
from astropy.table import Table
from astroquery.utils import TableList
import warnings
import calendar
import psycopg2
import math

import pandas as pd
from datetime import datetime
from astropy.io import fits
import json5


# base_dir = os.path.dirname(__file__)
# base_dir = os.path.dirname(os.path.abspath('__file__'))

PIPELINE_DEFAULT_DIR = '/mnt/photometry/'
FLAT_DEFAULT_DIR = '/mnt/windows_fits_data/fits_data_drive_backup/flat_storage_dir/'

PHOTOMETRY_MAG_LOWER_LIMIT = 12.5
PHOTOMETRY_MAG_UPPER_LIMIT = 23
PHOTOMETRY_QUERY_WIDTH = 48  # in arcmin, should cover the whole chip even with 45 deg rotation
PHOTOMETRY_QUERY_CATALOGS = {'VHS': ['J', 'II/367/'], 'VIKING': ['J', 'II/382/viking4'], 'VVV_J': ['J', 'II/348/vvv2'],
                             '2MASS': ['J', 'II/246/'],
                            'VVV_Y': ['J', 'II/348/vvv2'], 'Skymapper': ['Z', 'II/379/smssdr4'],'SDSS': ['Z', 'V/154/sdss16'],'DES_Z': ['Z', 'II/371/des_dr2'],
                             'VVV_Z': ['J', 'II/348/vvv2'], 'UKIDSS': ['Y', 'II/319/las9'],'PanSTARRS': ['Y', 'II/349/ps1'],'DES_Y': ['Y', 'II/371/des_dr2']
                             }

# PHOTOMETRY_QUERY_FUNCTIONS = {
#     'Z': ['viking_query', 'vvv_query', 'vhs_query', 'panstarrs_query', 'sdss_query', 'des_query'],
#     'Y': ['viking_query', 'vvv_query', 'vhs_query', 'ukidss_query', 'panstarrs_query', 'des_query'],
#     'J': ['viking_query', 'vvv_query', 'vhs_query', 'ukidss_query', 'twomass_query'],
#     'H': ['viking_query', 'vvv_query', 'vhs_query', 'ukidss_query', 'twomass_query']
# }

ASTROM_QUERY_CATALOGS = {
    'VIRAC': ['J,H','II/364/virac2'],
    'GAIA': ['J,H,Y,Z','I/350/gaiaedr3']
}

# AB mag surveys: DES, Skymapper, SDSS, & PSTARRS
AB_OFFSET_DICT = {'J': 0.94, 'H': 1.38, 'Y': 0.62, 'Z': 0.52}

PHOTOMETRY_LIM_MAGS = {
    'VHS': 19.5, 'VIKING': 21.2, 'VVV': 20.2, '2MASS': 15.5, 'Skymapper': 22, 'DES_Z': 23.1, 'SDSS': 23,
    'LAS': 21.0, 'UHS': 20.5, 'GPS': 21.0, 'GCS': 20.6, 'DES_Y': 21.7, 'PanSTARRS': 21.4
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

# keep magtype name *BELOW 4 CHARS*
MAGTYPES = {
    "AUTO": "AUTO",
    "PSF": "PSF",
    "APER": "APER"
}

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
    """Generates filepath for configs in .prime/config/, unless an existing filepath is specified"""
    # base_dir = os.path.dirname(__file__)  # os.path.abspath('__file__')
    # base_dir = gen_pipeline_file_name()
    if os.path.isfile(filename):
        return filename
    else:
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
    global FLAT_DEFAULT_DIR
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
    FLAT_DEFAULT_DIR = settings['FLAT_DEFAULT_DIR']
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
    weightmap_directory = os.path.join(gen_user_prime_dir(), 'weightmaps')
    return os.path.join(weightmap_directory, filename)


def gen_master_name():
    base_dir = os.path.dirname(__file__)
    return os.path.join(base_dir, 'master.py')


def gen_flat_dir():
    flat_dir = FLAT_DEFAULT_DIR
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
                    # return datetime.strptime(half.split(".")[2], "%Y%m%d")
                    given_date = datetime.strptime(half.split(".")[2], "%Y%m%d")
                    days_in_month = calendar.monthrange(given_date.year, given_date.month)[1]
                    target_date = given_date.replace(day=days_in_month // 2)    # middle of month
                    return target_date
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
        old = conf.server
        try:
            conf.server = url

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                v = Vizier(columns=["*"],
                           column_filters={'Jmag': '>12.5'},
                           row_limit=1)
                result = v.query_constraints(catalog=test_cat)
            print(f'Vizier mirror: {url}')
            break

        except Exception as e:
            print(f' {url} mirror check failed: \n{e}\n')
            conf.server = old
    else:
        raise RuntimeError('No working Vizier mirror found!')

    return url


# local 2mass query functions

def ang_convert(ang):
    """
    Converts input to appropriate degree equivalent for query (ex. 33m -> 33 arcmin to degree).  Use just a number or 'd' for
    degrees, m for arcmin, and s for arcsec
    """

    if isinstance(ang, (float, int)):
        arcconvert = ang
    elif not any(char.isalpha() for char in ang) or 'd' in ang:
        if 'd' in ang:
            ang = ang.split('d')[0]
        arcconvert = float(ang)
    elif 'm' in ang:
        angsplit = ang.split('m')
        rad_f = float(angsplit[0])
        arcconvert = rad_f / 60
    elif 's' in ang:
        angsplit = ang.split('s')
        rad_f = float(angsplit[0])
        arcconvert = rad_f / 3600
    else:
        print('Only arcsec, arcmin, and deg are supported! Default = deg')
        raise Exception('Use supported units.')
    return arcconvert


def local_query_box(ra_center, dec_center, dbname, tablename,
                    width, height=None,
                    columns=None,
                    column_filters=None):

    """
    Replicates astroquery's query_region box query functionality using the local psql db.

    Parameters
    ----------
    ra_center: float
        RA coordinate for query
    dec_center: float
        Dec coordinate for query
    width: float, int, or str
        full width of box query. If float or int or str w/ 'd' (ex. '2d'), degrees are assumed, if str w/ 'm' (ex. 33m),
        arcmin assumed, & if str w/ 's' (ex. '4s'), arcsec assumed
    height: float, int, or str
        optional, full height of box query.  If not specified, will be the same as width
    columns: list of str
        optional, specify specific cols to return w/ query. Same format as query_region
    column_filters: dict
        optional, apply filters to specific cols for query. Same format as query_region
    """

    ra_center = ra_center % 360

    if dec_center > 90 or dec_center < -90:
        raise ValueError("Declination must be between -90 and +90 degrees")

    if height is None:
        height = width

    print(f' Local query box size: {height, width}')
    height = ang_convert(height)
    width = ang_convert(width)

    POLE_THRESHOLD = 0.5  # degrees from pole to switch query strategy (box to radial)
    near_south_pole = dec_center - height / 2 <= -90 + POLE_THRESHOLD
    near_north_pole = dec_center + height / 2 >= 90 - POLE_THRESHOLD

    if near_south_pole or near_north_pole:
        # Use radial query near poles
        radius = math.sqrt((height / 2) ** 2 + (width / 2) ** 2)
        print(f" Near pole detected, switching to radial query with radius {radius:.4f} deg")
        query_type = 'radial'
        query_params = (ra_center, dec_center, radius)
    else:
        delta_ra = (width / 2) / math.cos(math.radians(dec_center))
        delta_dec = height / 2
        corners = [
            (ra_center - delta_ra, dec_center + delta_dec),  # TL
            (ra_center + delta_ra, dec_center + delta_dec),  # TR
            (ra_center + delta_ra, dec_center - delta_dec),  # BR
            (ra_center - delta_ra, dec_center - delta_dec),  # BL
        ]
        corners = [(ra % 360, dec) for ra, dec in corners]
        query_type = 'poly'
        query_params = corners

    conn = psycopg2.connect(service=dbname)

    cur = conn.cursor()

    colstr = "*" if columns is None else ",".join(columns)

    if query_type == 'radial':
        where = [f"""
           q3c_radial_query(
                ra, dec, {query_params[0]}, {query_params[1]}, {query_params[2]}
           )
           """]
    else:
        where = [f"""
           q3c_poly_query(
               ra, dec,
               ARRAY[
                   {query_params[0][0]}, {query_params[0][1]},
                   {query_params[1][0]}, {query_params[1][1]},
                   {query_params[2][0]}, {query_params[2][1]},
                   {query_params[3][0]}, {query_params[3][1]}
               ]::double precision[]
           )
           """]

    if column_filters:
        for col, expr in column_filters.items():
            if expr == '!= null':
                adj_expr = 'IS NOT NULL'
            else:
                adj_expr = expr
            where.append(f"{col} {adj_expr}")

    if columns:
        print(' Specific cols specified, removing rows where column vals = None')
        for col in columns:
            where.append(f"{col} IS NOT NULL")

    sql = f"""
        SELECT {colstr}
        FROM {tablename}
        WHERE {" AND ".join(where)}
    """

    # print(f' Constraints: {where}')

    cur.execute(sql)
    rows = cur.fetchall()
    names = [d[0] for d in cur.description]
    conn.close()

    tbl = Table(rows=rows, names=names) if rows else Table(names=names)

    return TableList([("local_query", tbl)])
