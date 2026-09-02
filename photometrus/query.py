"""
Query functions: for individual catalogs & automatic query
"""

import numpy as np
import pandas as pd
import astropy.units as u
from astroquery.vizier import Vizier
from astropy.coordinates import Angle, SkyCoord
from requests.exceptions import ConnectionError, Timeout
from astroquery.exceptions import RemoteServiceError
from astropy.io import ascii
from astropy.io import fits
from astropy.table import Table, Column
import psycopg2

from photometrus.settings import (PHOTOMETRY_MAG_LOWER_LIMIT, PHOTOMETRY_MAG_UPPER_LIMIT, AB_OFFSET_DICT,
                                  PHOTOMETRY_QUERY_WIDTH, local_query_box, gen_config_file_name)

from photometrus.utils.utils import convert_table_to_ldac
from photometrus.utils.defaults import PROCESSING_DEFAULTS as defaults


# ALL CATALOG QUERY FUNCTIONS
"""
Below is the format of an example Vizier query, if you want to add one.  It should have the following parameters:

coords: SkyCoord
    SkyCoord object of input center coordinates for query.  Normal usage with the main query fctn at the bottom of this 
    file usually provides this as: coords = SkyCoord(ra=raImage * u.degree, dec=decImage * u.degree, frame='fk5')
frame_long_str: str
    String providing the coordinate longitude. Normal usage with the main query fctn will provide this, usually as:
    frame_long_str = 'RA: %.4f' % raImage
frame_lat_str: str
    String providing the coordinate latitude. Normal usage with the main query fctn will provide this, usually as:
    frame_long_str = 'DEC: %.4f' % decImage
band: str
    Filter of observation
width: float
    Full width of query box in arcmin
mag_low_cutoff: float
    Bright end mag cutoff to avoid saturated sources (specified in settings.py, usually = 12.5)
mag_high_cutoff: float
    Dim end mag cutoff to handle very deep surveys beyond PRIME's limit (specified in settings.py, usually = 23)
chosen_frame: str
    String denoting coordinate system frame.  Normal usage with the main query fctn will provide this, usually as 'fk5'
    
*EXCEPTIONS TO THIS ARE THE ASTROMETRIC-BASED GAIA QUERIES*

Below is the formatting for a query function.  Note the following:

* The Vizier object should have 4 columns: [RA, DEC, Mag, MagErr].  While a 5th column is supported for Vista surveys, 
(Mclass), photometry.py currently only supports this specific column with those surveys to prune galaxies for zp calc.

* Column filters should feature a mag cutoff of bright sources, as in the example below.  If the survey is especially 
deep, it should also have a dim mag cutoff.  See the skymapper or sdss query fctns for formatting in that case.  Column
filters should also have whatever error flag constraints you deem necessary to prune unwanted sources.

* Your function should return two parameters: Q (the query object from query_region), and survey_name (name of the survey
you're adding).  The survey_name should match the name of the fctn, ex: survey_name='VHS' for the vhs_query() fctn.

* Finally, make sure to add your query to the PHOTOMETRY_QUERY_FUNCTIONS dict near the bottom of the file.  Add it for 
the applicable bands the survey provides.  The order of the survey fctns in the dict denotes hierarchy of attempts, ex. 
for J band VIKING will be tried first, then VVV, and so on.


def example_query(coords, frame_long_str, frame_lat_str, band, width, mag_low_cutoff,
              mag_high_cutoff=PHOTOMETRY_MAG_UPPER_LIMIT, chosen_frame='fk5'):

    survey_name = 'name_of_input_survey' (i.e. 'VHS')
    catNum = 'Vizier_catalog_number_for_survey' (i.e. 'II/367' for VHS)
    print('\nQuerying Vizier %s around %s, %s, boxwidth %.2f arcmin, mag lim of %s - %s'
          % (catNum, frame_long_str, frame_lat_str, width, mag_low_cutoff, mag_high_cutoff))
    try:
        v = Vizier(columns=['RAJ2000', 'DEJ2000', '%sap3' % band, 'e_%sap3' % band],
                   column_filters={
                                "%sap3" % band: f">{mag_low_cutoff:f}",
                                "%sperrbits" % band: '<128'},
                   row_limit=-1)
        Q = v.query_region(SkyCoord(coords, unit=(u.deg, u.deg)), width=str(width) + 'm'
                           , catalog=catNum, cache=False, frame=chosen_frame)
    except (RemoteServiceError, ConnectionError, Timeout) as e:
        print(f'Error: {e}')
        print(
            'Error in Vizier query')
    return Q, survey_name              
"""


def gaia_vizier_query(coords, width, chosen_frame='fk5'):
    """Vizier GAIA query function"""

    if chosen_frame != 'fk5':
        raise Exception('Frame other than fk5 detected! *WARNING* Local query currently doesnt '
                        'support galactic coords!')

    survey_name = 'GAIA'
    catNum = 'I/350/gaiaedr3'
    try:
        v = Vizier(columns=['RA_ICRS', 'DE_ICRS', 'RPmag'],
                   column_filters={"Dup": "<1", "Nd": ">6"},
                   row_limit=-1)
        Q = v.query_region(SkyCoord(coords, unit=(u.deg, u.deg)), width=str(width) + 'm'
                           , catalog=catNum, cache=False, frame=chosen_frame)
    except (RemoteServiceError, ConnectionError, Timeout) as e:
        print(f'Error: {e}')
        print(
            'Error in Vizier query.  Is Vizier down?')

    return Q, survey_name


def gaia_query(coords, width, chosen_frame='fk5'):
    """Local GAIA query function"""

    if chosen_frame != 'fk5':
        raise Exception('Frame other than fk5 detected! *WARNING* Local query currently doesnt '
                        'support galactic coords!')

    ra = coords.ra.deg
    dec = coords.dec.deg

    survey_name = 'GAIA'
    print('Querying GAIA EDR3 locally')
    try:
        Q = local_query_box(ra_center=ra,
                            dec_center=dec,
                            dbname='prime_vhs_local',
                            tablename='gaia_sources',
                            width=str(width) + 'm',
                            columns=['ra', 'dec', 'ra_error', 'dec_error', 'phot_rp_mean_mag', 'ref_epoch'],
                            column_filters={
                                "duplicated_source": "IS FALSE"
                            }
                            )
    except (psycopg2.ProgrammingError, psycopg2.OperationalError) as e:
        print(f'Error: {e}')
        print(f'Local GAIA query unsuccessful! Cannot continue!')
        print('Attempting Vizier query as backup!')
        Q, survey_name = gaia_vizier_query(coords, width, chosen_frame=chosen_frame)

    return Q, survey_name


def twomass_vizier_query(coords, frame_long_str, frame_lat_str, band, width, mag_low_cutoff,
                  mag_high_cutoff=PHOTOMETRY_MAG_UPPER_LIMIT, chosen_frame='fk5', errbits=None):
    """Vizier 2mass query function"""

    survey_name = '2MASS'
    catNum = 'II/246'
    print('\nQuerying Vizier %s around %s, %s, boxwidth %.2f arcmin, mag lim of %s - %s'
          % (catNum, frame_long_str, frame_lat_str, width, mag_low_cutoff, mag_high_cutoff))

    column_filters = {
        "%smag" % band: f">{mag_low_cutoff:f}",
        "Nd": ">6",
    }

    # errbits constraint
    if errbits is not None:
        column_filters["Cflg"] = f"{errbits}"

    try:
        v = Vizier(columns=['RAJ2000', 'DEJ2000', '%smag' % band, 'e_%smag' % band],
                   column_filters=column_filters,
                   row_limit=-1)
        Q = v.query_region(SkyCoord(coords, unit=(u.deg, u.deg)), width=str(width) + 'm'
                           , catalog=catNum, cache=False, frame=chosen_frame)
    except (RemoteServiceError, ConnectionError, Timeout) as e:
        print(f'Error: {e}')
        print(
            'Error in Vizier query.  Is Vizier down?')

    return Q, survey_name


def twomass_query(coords, frame_long_str, frame_lat_str, band, width, mag_low_cutoff,
                  mag_high_cutoff=PHOTOMETRY_MAG_UPPER_LIMIT, chosen_frame='fk5', errbits=None):
    """Local 2mass query function"""

    if chosen_frame != 'fk5':
        raise Exception('Frame other than fk5 detected! *WARNING* Local query currently doesnt '
                        'support galactic coords!')

    ra = coords.ra.deg
    dec = coords.dec.deg

    column_filters = {
        f"{band.lower()}mag": f">{mag_low_cutoff:f}",
    }

    # errbits constraint
    if errbits is not None:
        column_filters["cc_flg"] = f"{errbits}"

    survey_name = '2MASS'
    print('\nVizier catalogs exhausted, switching to local 2MASS query...')
    print('Local 2MASS Query around %s, %s, boxwidth %.2f arcmin, mag lim of %s - %s'
          % (frame_long_str, frame_lat_str, width, mag_low_cutoff, mag_high_cutoff))
    try:
        Q = local_query_box(ra_center=ra,
                            dec_center=dec,
                            dbname='prime_vhs_local',
                            tablename='twomass_sources',
                            width=str(width) + 'm',
                            columns=["ra", "dec", f"{band.lower()}mag", f"e_{band.lower()}mag"],
                            column_filters=column_filters
                            )
    except (psycopg2.ProgrammingError, psycopg2.OperationalError) as e:
        print(f'Error: {e}')
        print(f'Local 2MASS query unsuccessful!  Is the field in Y or Z band? '
                        ' If so, and there are no other surveys, 2MASS does not have'
                        ' these filters!  Cannot continue with photometry!')
        print('Attempting Vizier query as backup!')
        Q, survey_name = twomass_vizier_query(coords, frame_long_str, frame_lat_str, band, width, mag_low_cutoff,
              mag_high_cutoff=mag_high_cutoff, chosen_frame=chosen_frame)

    return Q, survey_name


def vhs_vizier_query(coords, frame_long_str, frame_lat_str, band, width, mag_low_cutoff,
              mag_high_cutoff=PHOTOMETRY_MAG_UPPER_LIMIT, chosen_frame='fk5', errbits=None):
    """Vizier VHS query function"""

    survey_name = 'VHS'
    catNum = 'II/367'
    print('\nQuerying Vizier %s around %s, %s, boxwidth %.2f arcmin, mag lim of %s - %s'
          % (catNum, frame_long_str, frame_lat_str, width, mag_low_cutoff, mag_high_cutoff))

    column_filters = {
        "%sap3" % band: f">{mag_low_cutoff:f}",
    }

    # errbits constraint
    errbit_flg = '<128'
    if errbits is not None:
        errbit_flg = errbits

    column_filters[f"{band}perrbits"] = errbit_flg

    try:
        v = Vizier(columns=['RAJ2000', 'DEJ2000', '%sap3' % band, 'e_%sap3' % band, 'Mclass'],
                   column_filters=column_filters,
                   row_limit=-1)
        Q = v.query_region(SkyCoord(coords, unit=(u.deg, u.deg)), width=str(width) + 'm'
                           , catalog=catNum, cache=False, frame=chosen_frame)
    except (RemoteServiceError, ConnectionError, Timeout) as e:
        print(f'Error: {e}')
        print(
            'Error in Vizier query. Perhaps your image is not in the southern hemisphere sky?  '
            'H band is also not well covered!')
    return Q, survey_name


def vhs_query(coords, frame_long_str, frame_lat_str, band, width, mag_low_cutoff,
              mag_high_cutoff=PHOTOMETRY_MAG_UPPER_LIMIT, chosen_frame='fk5', errbits=None):
    """VHS query function"""

    if chosen_frame != 'fk5':
        raise Exception('Frame other than fk5 detected! *WARNING* Local query currently doesnt '
                        'support galactic coords!')

    ra = coords.ra.deg
    dec = coords.dec.deg

    column_filters = {
        f'{band}AperMag3': f">{mag_low_cutoff:f}",
    }

    # errbits constraint
    errbit_flg = '<128'
    if errbits is not None:
        errbit_flg = errbits

    column_filters[f"{band}ppErrBits"] = errbit_flg

    survey_name = 'VHS'
    print('\nLocal VHS Query around %s, %s, boxwidth %.2f arcmin, mag lim of %s - %s'
          % (frame_long_str, frame_lat_str, width, mag_low_cutoff, mag_high_cutoff))
    try:
        Q = local_query_box(ra_center=ra,
                            dec_center=dec,
                            dbname='prime_vhs_local',
                            tablename='vhs_sources',
                            width=str(width) + 'm',
                            columns=["ra", "dec", f'{band}AperMag3', f'{band}AperMag3Err', 'mergedClass'],
                            column_filters=column_filters
                            )
    except (psycopg2.ProgrammingError, psycopg2.OperationalError) as e:
        print(f'Local VHS query unsuccessful!: Error: {e}')
        print('Attempting Vizier query as backup!')
        Q, survey_name = vhs_vizier_query(coords, frame_long_str, frame_lat_str, band, width, mag_low_cutoff,
              mag_high_cutoff=mag_high_cutoff, chosen_frame=chosen_frame)

    return Q, survey_name


def viking_vizier_query(coords, frame_long_str, frame_lat_str, band, width, mag_low_cutoff,
                    mag_high_cutoff=PHOTOMETRY_MAG_UPPER_LIMIT, chosen_frame='fk5', errbits=None):
    """Vizier VIKING query function"""

    survey_name = 'VIKING'
    catNum = 'II/382/viking4'
    print('\nQuerying Vizier %s around %s, %s, boxwidth %.2f arcmin, mag lim of %s - %s'
          % (catNum, frame_long_str, frame_lat_str, width, mag_low_cutoff, mag_high_cutoff))

    column_filters = {
        "%sap3" % band: f">{mag_low_cutoff:f}",
    }

    # errbits constraint
    errbit_flg = '<128'
    if errbits is not None:
        errbit_flg = errbits
        column_filters[f"Hclass"] = '-1'

    column_filters[f"{band}perrbits"] = errbit_flg

    try:
        v = Vizier(columns=['RAJ2000', 'DEJ2000', '%sap3' % band, 'e_%sap3' % band, 'Mclass'],
                   column_filters=column_filters,
                   row_limit=-1)
        Q = v.query_region(SkyCoord(coords, unit=(u.deg, u.deg)), width=str(width) + 'm'
                           , catalog=catNum, cache=False, frame=chosen_frame)
    except (RemoteServiceError, ConnectionError, Timeout) as e:
        print(f'Error: {e}')
        print(
            'Error in Vizier query. Perhaps your image is not in the southern hemisphere sky?  '
            'H band is also not well covered!'
            ' If you are in S.H., VIKING is only in a relatively smaller strip!')
    return Q, survey_name


def viking_query(coords, frame_long_str, frame_lat_str, band, width, mag_low_cutoff,
              mag_high_cutoff=PHOTOMETRY_MAG_UPPER_LIMIT, chosen_frame='fk5', errbits=None):
    """VIKING query function"""

    if chosen_frame != 'fk5':
        raise Exception('Frame other than fk5 detected! *WARNING* Local query currently doesnt '
                        'support galactic coords!')

    ra = coords.ra.deg
    dec = coords.dec.deg

    column_filters = {
        f'{band}AperMag3': f">{mag_low_cutoff:f}",
    }

    # errbits constraint
    errbit_flg = '<128'
    if errbits is not None:
        errbit_flg = errbits
        column_filters[f"hClass"] = '-1'

    column_filters[f"{band}ppErrBits"] = errbit_flg

    survey_name = 'VIKING'
    print('\nLocal VIKING Query around %s, %s, boxwidth %.2f arcmin, mag lim of %s - %s'
          % (frame_long_str, frame_lat_str, width, mag_low_cutoff, mag_high_cutoff))
    try:
        Q = local_query_box(ra_center=ra,
                            dec_center=dec,
                            dbname='prime_vhs_local',
                            tablename='viking_sources',
                            width=str(width) + 'm',
                            columns=["ra", "dec", f'{band}AperMag3', f'{band}AperMag3Err', 'mergedClass'],
                            column_filters=column_filters
                            )
    except (psycopg2.ProgrammingError, psycopg2.OperationalError) as e:
        print(f'Local VIKING query unsuccessful!: Error: {e}')
        print('Attempting Vizier query as backup!')
        Q, survey_name = viking_vizier_query(coords, frame_long_str, frame_lat_str, band, width, mag_low_cutoff,
              mag_high_cutoff=mag_high_cutoff, chosen_frame=chosen_frame)

    return Q, survey_name


def vvv_query(coords, frame_long_str, frame_lat_str, band, width, mag_low_cutoff,
                    mag_high_cutoff=PHOTOMETRY_MAG_UPPER_LIMIT, chosen_frame='fk5', errbits=None):
    """Vizier VVV query function"""

    survey_name = 'VVV'
    catNum = 'II/348/vvv2'
    print('\nQuerying Vizier %s around %s, %s, boxwidth %.2f arcmin, mag lim of %s - %s'
          % (catNum, frame_long_str, frame_lat_str, width, mag_low_cutoff, mag_high_cutoff))

    column_filters = {
        "%sap3" % band: f">{mag_low_cutoff:f}",
    }

    # errbits constraint
    errbit_flg = '<128'
    if errbits is not None:
        errbit_flg = errbits

    column_filters[f"{band}perrbits"] = errbit_flg

    try:
        v = Vizier(columns=['RAJ2000', 'DEJ2000', '%sap3' % band, 'e_%sap3' % band, 'Mclass'],
                   column_filters=column_filters,
                   row_limit=-1)
        Q = v.query_region(SkyCoord(coords, unit=(u.deg, u.deg)), width=str(width) + 'm'
                           , catalog=catNum, cache=False, frame=chosen_frame)
    except (RemoteServiceError, ConnectionError, Timeout) as e:
        print(f'Error: {e}')
        print(
            'Error in Vizier query. Perhaps your image is not in the southern hemisphere sky?  '
            'H band is also not well covered!')
    return Q, survey_name


def las_vizier_query(coords, frame_long_str, frame_lat_str, band, width, mag_low_cutoff,
                    mag_high_cutoff=PHOTOMETRY_MAG_UPPER_LIMIT, chosen_frame='fk5', errbits=None):
    """Vizier UKIDSS LAS query function"""

    survey_name = 'LAS'
    catNum = 'II/319/las9'
    print('\nQuerying Vizier %s around %s, %s, boxwidth %.2f arcmin, mag lim of %s - %s'
          % (catNum, frame_long_str, frame_lat_str, width, mag_low_cutoff, mag_high_cutoff))

    mag_col_name = f'{band}mag'
    errbit_col_name = f"{band}flags"
    if band == 'J':
        mag_col_name = f'{band}mag1'
        errbit_col_name = f"{band}flags1"

    column_filters = {
        "%s" % mag_col_name: f">{mag_low_cutoff:f}",
    }

    # errbits constraint
    errbit_flg = '<16'
    if errbits is not None:
        errbit_flg = errbits

    column_filters[errbit_col_name] = errbit_flg

    try:
        v = Vizier(columns=['RAJ2000', 'DEJ2000', '%s' % mag_col_name, 'e_%s' % mag_col_name],
                   column_filters=column_filters,
                   row_limit=-1)
        Q = v.query_region(SkyCoord(coords, unit=(u.deg, u.deg)), width=str(width) + 'm',
                           catalog=catNum, cache=False, frame=chosen_frame)
    except (RemoteServiceError, ConnectionError, Timeout) as e:
        print(f'Error: {e}')
        print(
            'Error in Vizier query. Perhaps your image is not in the southern hemisphere sky?'
            '\n perhaps check UKIDSS coverage maps?')

    return Q, survey_name


def las_query(coords, frame_long_str, frame_lat_str, band, width, mag_low_cutoff,
              mag_high_cutoff=PHOTOMETRY_MAG_UPPER_LIMIT, chosen_frame='fk5', errbits=None):
    """UKIDSS LAS query function"""

    if chosen_frame != 'fk5':
        raise Exception('Frame other than fk5 detected! *WARNING* Local query currently doesnt '
                        'support galactic coords!')

    ra = coords.ra.deg
    dec = coords.dec.deg

    column_filters = {
        f'{band}AperMag3': f">{mag_low_cutoff:f}",
    }

    # errbits constraint
    errbit_flg = '<128'
    if errbits is not None:
        errbit_flg = errbits

    column_filters[f"{band}ppErrBits"] = errbit_flg

    survey_name = 'LAS'
    print('\nLocal LAS Query around %s, %s, boxwidth %.2f arcmin, mag lim of %s - %s'
          % (frame_long_str, frame_lat_str, width, mag_low_cutoff, mag_high_cutoff))
    try:
        Q = local_query_box(ra_center=ra,
                            dec_center=dec,
                            dbname='prime_vhs_local',
                            tablename='las_sources',
                            width=str(width) + 'm',
                            columns=["ra", "dec", f'{band}AperMag3', f'{band}AperMag3Err', 'mergedClass'],
                            column_filters=column_filters
                            )
    except (psycopg2.ProgrammingError, psycopg2.OperationalError) as e:
        print(f'Local LAS query unsuccessful!: Error: {e}')
        print('Attempting Vizier query as backup!')
        Q, survey_name = las_vizier_query(coords, frame_long_str, frame_lat_str, band, width, mag_low_cutoff,
              mag_high_cutoff=mag_high_cutoff, chosen_frame=chosen_frame)

    return Q, survey_name


def uhs_query(coords, frame_long_str, frame_lat_str, band, width, mag_low_cutoff,
              mag_high_cutoff=PHOTOMETRY_MAG_UPPER_LIMIT, chosen_frame='fk5', errbits=None):
    """UKIDSS UHS query function"""

    if chosen_frame != 'fk5':
        raise Exception('Frame other than fk5 detected! *WARNING* Local query currently doesnt '
                        'support galactic coords!')

    ra = coords.ra.deg
    dec = coords.dec.deg

    column_filters = {
        f'{band}AperMag3': f">{mag_low_cutoff:f}",
    }

    # errbits constraint
    errbit_flg = '<128'
    if errbits is not None:
        errbit_flg = errbits

    column_filters[f"{band}ppErrBits"] = errbit_flg

    survey_name = 'UHS'
    print('\nLocal UHS Query around %s, %s, boxwidth %.2f arcmin, mag lim of %s - %s'
          % (frame_long_str, frame_lat_str, width, mag_low_cutoff, mag_high_cutoff))
    try:
        Q = local_query_box(ra_center=ra,
                            dec_center=dec,
                            dbname='prime_vhs_local',
                            tablename='uhs_sources',
                            width=str(width) + 'm',
                            columns=["ra", "dec", f'{band}AperMag3', f'{band}AperMag3Err', 'mergedClass'],
                            column_filters=column_filters
                            )
    except (psycopg2.ProgrammingError, psycopg2.OperationalError) as e:
        print(f'Local UHS query unsuccessful!: Error: {e}')
        raise Exception('UKIDSS UHS is not available with vizier!')

    return Q, survey_name


def gps_query(coords, frame_long_str, frame_lat_str, band, width, mag_low_cutoff,
              mag_high_cutoff=PHOTOMETRY_MAG_UPPER_LIMIT, chosen_frame='fk5', errbits=None):
    """UKIDSS GPS query function"""

    if chosen_frame != 'fk5':
        raise Exception('Frame other than fk5 detected! *WARNING* Local query currently doesnt '
                        'support galactic coords!')

    ra = coords.ra.deg
    dec = coords.dec.deg

    column_filters = {
        f'{band}AperMag3': f">{mag_low_cutoff:f}",
    }

    # errbits constraint
    errbit_flg = '<128'
    if errbits is not None:
        errbit_flg = errbits

    column_filters[f"{band}ppErrBits"] = errbit_flg

    survey_name = 'GPS'
    print('\nLocal GPS Query around %s, %s, boxwidth %.2f arcmin, mag lim of %s - %s'
          % (frame_long_str, frame_lat_str, width, mag_low_cutoff, mag_high_cutoff))
    try:
        Q = local_query_box(ra_center=ra,
                            dec_center=dec,
                            dbname='prime_vhs_local',
                            tablename='gps_sources',
                            width=str(width) + 'm',
                            columns=["ra", "dec", f'{band}AperMag3', f'{band}AperMag3Err', 'mergedClass'],
                            column_filters=column_filters
                            )
    except (psycopg2.ProgrammingError, psycopg2.OperationalError) as e:
        print(f'Local GPS query unsuccessful!: Error: {e}')
        raise Exception('UKIDSS GPS is not available with vizier!')

    return Q, survey_name


def gcs_query(coords, frame_long_str, frame_lat_str, band, width, mag_low_cutoff,
              mag_high_cutoff=PHOTOMETRY_MAG_UPPER_LIMIT, chosen_frame='fk5', errbits=None):
    """UKIDSS GCS query function"""

    if chosen_frame != 'fk5':
        raise Exception('Frame other than fk5 detected! *WARNING* Local query currently doesnt '
                        'support galactic coords!')

    ra = coords.ra.deg
    dec = coords.dec.deg

    column_filters = {
        f'{band}AperMag3': f">{mag_low_cutoff:f}",
    }

    # errbits constraint
    errbit_flg = '<128'
    if errbits is not None:
        errbit_flg = errbits

    column_filters[f"{band}ppErrBits"] = errbit_flg

    survey_name = 'GCS'
    print('\nLocal GCS Query around %s, %s, boxwidth %.2f arcmin, mag lim of %s - %s'
          % (frame_long_str, frame_lat_str, width, mag_low_cutoff, mag_high_cutoff))
    try:
        Q = local_query_box(ra_center=ra,
                            dec_center=dec,
                            dbname='prime_vhs_local',
                            tablename='gcs_sources',
                            width=str(width) + 'm',
                            columns=["ra", "dec", f'{band}AperMag3', f'{band}AperMag3Err', 'mergedClass'],
                            column_filters=column_filters
                            )
    except (psycopg2.ProgrammingError, psycopg2.OperationalError) as e:
        print(f'Local GCS query unsuccessful!: Error: {e}')
        raise Exception('UKIDSS GCS is not available with vizier!')

    return Q, survey_name


def panstarrs_query(coords, frame_long_str, frame_lat_str, band, width, mag_low_cutoff,
                    mag_high_cutoff=PHOTOMETRY_MAG_UPPER_LIMIT, chosen_frame='fk5', errbits=None):
    """Vizier Pan-STARRS query function"""

    survey_name = 'PanSTARRS'
    catNum = 'II/389/ps1_dr2'
    print('\nQuerying Vizier %s around %s, %s, boxwidth %.2f arcmin, mag lim of %s - %s'
          % (catNum, frame_long_str, frame_lat_str, width, mag_low_cutoff, mag_high_cutoff))

    column_filters = {
        "%smag" % band: f">{mag_low_cutoff:f}",
    }

    # errbits constraint
    errbit_flg = "!=1 && !=2"
    if errbits is not None:
        # errbit_flg = errbits
        column_filters['Nd'] = '>6'

    column_filters[f"{band.lower()}Flags"] = errbit_flg

    try:
        v = Vizier(columns=['RAJ2000', 'DEJ2000', '%smag' % band.lower(), 'e_%smag' % band.lower()],
                   column_filters=column_filters,
                   row_limit=-1)
        Q = v.query_region(SkyCoord(coords, unit=(u.deg, u.deg)), width=str(width) + 'm',
                           catalog=catNum, cache=False, frame=chosen_frame)
    except (RemoteServiceError, ConnectionError, Timeout) as e:
        print(f'Error: {e}')
        print(
            'Error in Vizier query. Perhaps your image is not in the southern hemisphere sky?')
    return Q, survey_name


def skymapper_query(coords, frame_long_str, frame_lat_str, band, width, mag_low_cutoff,
                    mag_high_cutoff=PHOTOMETRY_MAG_UPPER_LIMIT, chosen_frame='fk5', errbits=None):
    """Vizier Skymapper query function"""

    survey_name = 'Skymapper'
    catNum = 'II/379/smssdr4'
    print('\nQuerying Vizier %s around %s, %s, boxwidth %.2f arcmin, mag lim of %s - %s'
          % (catNum, frame_long_str, frame_lat_str, width, mag_low_cutoff, mag_high_cutoff))

    column_filters = {
        "%sPSF" % band.lower(): f"{mag_low_cutoff:f}..{mag_high_cutoff:f}",
    }

    # errbits constraint
    errbit_flg = "<4"
    if errbits is not None:
        # errbit_flg = errbits
        column_filters['Nd'] = '>6'

    column_filters[f"{band.lower()}Flag"] = errbit_flg

    try:
        v = Vizier(columns=['RAICRS', 'DEICRS', '%sPSF' % band.lower(), 'e_%sPSF' % band.lower()],
                   column_filters=column_filters,
                   row_limit=-1)
        Q = v.query_region(SkyCoord(coords, unit=(u.deg, u.deg)), width=str(width) + 'm'
                           , catalog=catNum, cache=False, frame=chosen_frame)
    except (RemoteServiceError, ConnectionError, Timeout) as e:
        print(f'Error: {e}')
        print(
            'Error in Vizier query. Perhaps your image is not in the southern hemisphere sky?'
            ' Perhaps check skymapper coverage maps?')
    return Q, survey_name


def sdss_query(coords, frame_long_str, frame_lat_str, band, width, mag_low_cutoff,
                    mag_high_cutoff=PHOTOMETRY_MAG_UPPER_LIMIT, chosen_frame='fk5', errbits=None):
    """Vizier SDSS query function"""

    survey_name = 'SDSS'
    catNum = 'V/154/sdss16'
    print('\nQuerying Vizier %s around %s, %s, boxwidth %.2f arcmin, mag lim of %s - %s'
          % (catNum, frame_long_str, frame_lat_str, width, mag_low_cutoff, mag_high_cutoff))

    column_filters = {
        "%spmag" % band.lower(): f"{mag_low_cutoff:f}..{mag_high_cutoff:f}",
        "%sFlag" % band.lower(): "<4",
        "clean": "=1",
    }

    # errbits constraint
    if errbits is not None:
        errbit_flg = errbits
        column_filters['Nd'] = '>6'
        column_filters[f"{band.lower()}Flags"] = errbit_flg

    try:
        v = Vizier(columns=['RA_ICRS', 'DE_ICRS', '%spmag' % band.lower(), 'e_%spmag' % band.lower()],
                   column_filters=column_filters,
                   row_limit=-1)
        Q = v.query_region(SkyCoord(coords, unit=(u.deg, u.deg)), width=str(width) + 'm'
                           , catalog=catNum, cache=False, frame=chosen_frame)
    except (RemoteServiceError, ConnectionError, Timeout) as e:
        print(f'Error: {e}')
        print(
            'Error in Vizier query. Perhaps your image is not in the southern hemisphere sky?'
            ' Perhaps check SDSS coverage maps?')
    return Q, survey_name


def des_query(coords, frame_long_str, frame_lat_str, band, width, mag_low_cutoff,
                    mag_high_cutoff=PHOTOMETRY_MAG_UPPER_LIMIT, chosen_frame='fk5', errbits=None):
    """Vizier DES query function"""

    if band == 'Z':
        band = band.lower()

    survey_name = 'DES'
    catNum = 'II/371/des_dr2'
    print('\nQuerying Vizier %s around %s, %s, boxwidth %.2f arcmin, mag lim of %s - %s'
          % (catNum, frame_long_str, frame_lat_str, width, mag_low_cutoff, mag_high_cutoff))

    column_filters = {
        "%smag" % band: f"{mag_low_cutoff:f}..{mag_high_cutoff:f}",
        "%sFlag" % band.lower(): "<4",
        "clean": "=1",
    }

    # errbits constraint
    # errbit_flg = "<4"
    # if errbits is not None:
    #     errbit_flg = errbits
    #     column_filters['Nd'] = '>6'
    #
    # column_filters[f"{band.lower()}Flags"] = errbit_flg

    try:
        v = Vizier(columns=['RA_ICRS', 'DE_ICRS', '%smag' % band, 'e_%smag' % band],
                   column_filters=column_filters,
                   row_limit=-1)
        Q = v.query_region(SkyCoord(coords, unit=(u.deg, u.deg)), width=str(width) + 'm',
                           catalog=catNum, cache=False, frame=chosen_frame)
    except (RemoteServiceError, ConnectionError, Timeout) as e:
        print(f'Error: {e}')
        print(
            'Error in Vizier query. Perhaps your image is not in the southern hemisphere sky?'
            ' Perhaps check DES coverage maps?')
    return Q, survey_name


PHOTOMETRY_QUERY_FUNCTIONS = {
    'Z': [viking_query, vvv_query, vhs_query, gcs_query, panstarrs_query, sdss_query, des_query],
    'Y': [viking_query, vvv_query, vhs_query, las_query, gcs_query, panstarrs_query, des_query],
    'J': [viking_query, vvv_query, vhs_query, las_query, uhs_query, gcs_query, twomass_query],
    'H': [viking_query, vvv_query, vhs_query, las_query, uhs_query, gcs_query, twomass_query]
}


# GAIA COMPLETION QUERY


def gaia_crsmtch_check(coords, width, chosen_frame, w, data, crop, Q):
    crop = int(crop)
    max_x = data.shape[0]
    max_y = data.shape[1]

    # gaia query
    mag_low_cutoff = 3
    catNum = 'I/350/gaiaedr3'
    # mag_lims = f">{mag_low_cutoff:f}"

    try:
        print(f' Querying {catNum} and crossmatching to determine catalog completion..')
        G, _ = gaia_query(coords=coords, width=width, chosen_frame=chosen_frame)

        # crsmtch check
        gaia_colnames = G[0].colnames
        G_RA = gaia_colnames[0]
        G_DEC = gaia_colnames[1]

        query_colnames = Q[0].colnames
        Q_RA = query_colnames[0]
        Q_DEC = query_colnames[1]

        G_imCoords = w.all_world2pix(G[0][G_RA], G[0][G_DEC], 1)
        Q_imCoords = w.all_world2pix(Q[0][Q_RA], Q[0][Q_DEC], 1)

        good_G_stars = G[0][
            np.where((G_imCoords[0] > crop) & (G_imCoords[0] < (max_x - crop)) & (G_imCoords[1] > crop) & (
                    G_imCoords[1] < (max_y - crop)))]
        good_Q_stars = Q[0][
            np.where((Q_imCoords[0] > crop) & (Q_imCoords[0] < (max_x - crop)) & (Q_imCoords[1] > crop) & (
                    Q_imCoords[1] < (max_y - crop)))]

        GaiaCatCoords = SkyCoord(ra=good_G_stars[G_RA], dec=good_G_stars[G_DEC], frame='icrs', unit='degree')
        QueryCatCoords = SkyCoord(ra=good_Q_stars[Q_RA], dec=good_Q_stars[Q_DEC], frame='icrs', unit='degree')

        print(' Gaia cropped source total = ', len(good_G_stars))
        print(f' Chosen survey cropped source total = ', len(good_Q_stars))

        gaia_crsmtch_thresh = 1.0
        idx_gaia, idx_query, d2d, d3d = QueryCatCoords.search_around_sky(GaiaCatCoords,
                                                                         gaia_crsmtch_thresh * u.arcsec)

        df = pd.DataFrame({
            'idx_gaia': idx_gaia,
            'idx_query': idx_query,
            'd2d': d2d.to(u.arcsec).value  # example in arcsec
        })

        # Sort by separation and drop duplicates of gaia index, keeping the closest
        gaia_matches_closest = df.sort_values('d2d').drop_duplicates('idx_gaia', keep='first')

        print(f' Crossmatched Gaia source num = {len(gaia_matches_closest)}')
        gaia_completion = len(gaia_matches_closest) / len(good_G_stars)
        print(' Completion = %.2f' % gaia_completion)
    except AttributeError:
        print(' Gaia sources not found!  Skipping completion check!')
        gaia_completion = 1
    except Exception as e:
        print(f'Error in Vizier GAIA query & Survey Crossmatch: {e}')
        print(' Assuming bad completion!')
        gaia_completion = 0

    return gaia_completion

# GENERAL PHOTOMETRY QUERY FUNCTION


def query(
        raImage, decImage, band, w, data, crop=defaults['crop'], acc_comp_lvl=defaults['gaia_acc_comp_lvl'], survey=None,
        given_catalog_path=None, mag_lower_lim=PHOTOMETRY_MAG_LOWER_LIMIT, mag_upper_lim=PHOTOMETRY_MAG_UPPER_LIMIT, bulge=False,
        magtype=defaults['magtype']
):

    mag_low_cutoff = mag_lower_lim
    mag_high_cutoff = mag_upper_lim

    # query box width
    width = PHOTOMETRY_QUERY_WIDTH

    if given_catalog_path:
        Q = ascii.read(given_catalog_path)
        Q = Q[Q[f'{band}MAG_{magtype}'] > mag_low_cutoff]
        chosen_survey = 'PRIME'

        return Q, chosen_survey, mag_low_cutoff

    else:
        if bulge:
            # if bulge field, change to galactic coords for query
            print('Galactic bulge field detected!  Adjusting query parameters accordingly...')
            coords = SkyCoord(ra=raImage * u.degree, dec=decImage * u.degree, frame='fk5')
            coords = coords.galactic  # galactic conversion for bulge fields
            chosen_frame = 'galactic'
            frame_long = coords.l.deg
            frame_long_str = 'l = %.4f' % frame_long
            frame_lat = coords.b.deg
            frame_lat_str = 'b = %.4f' % frame_lat
            print('Converting coords to galactic: %s, %s' % (frame_long_str, frame_lat_str))

        else:
            coords = SkyCoord(ra=raImage * u.degree, dec=decImage * u.degree, frame='fk5')
            chosen_frame = 'fk5'
            frame_long = raImage
            frame_long_str = 'RA: %.4f' % frame_long
            frame_lat = decImage
            frame_lat_str = 'DEC: %.4f' % frame_lat

    # specified survey query

    if survey:
        Q = None
        chosen_survey = survey

        fctns_for_band = [fctn for fctn in PHOTOMETRY_QUERY_FUNCTIONS[band]]
        chosen_fctn = [fctn for fctn in fctns_for_band if chosen_survey.lower() in fctn.__name__.lower()]

        Q, chosen_survey = chosen_fctn[0](
            coords=coords, frame_long_str=frame_long_str, frame_lat_str=frame_lat_str, band=band,
            width=width,mag_low_cutoff=mag_low_cutoff, mag_high_cutoff=mag_high_cutoff, chosen_frame=chosen_frame
        )

        if not Q:
            raise ReferenceError('Chosen catalog failed to produce a result!  No coverage?')
        else:
            if len(Q[0]) > 0:
                print('Queried source total = ', len(Q[0]))
            else:
                raise ValueError('No sources found in survey catalog!  No coverage?')

        return Q, chosen_survey, mag_low_cutoff

    # auto query

    fctns_for_band = [fctn for fctn in PHOTOMETRY_QUERY_FUNCTIONS[band]]
    last_idx = fctns_for_band[-1]

    Q = None
    chosen_survey = None

    for fctn in fctns_for_band:
        Q, chosen_survey = fctn(
            coords=coords, frame_long_str=frame_long_str, frame_lat_str=frame_lat_str, band=band,
            width=width,mag_low_cutoff=mag_low_cutoff, mag_high_cutoff=mag_high_cutoff, chosen_frame=chosen_frame
        )

        if Q:
            if len(Q[0]) > 0:
                print('Queried source total = ', len(Q[0]))
                if fctn != last_idx:
                    gaia_comp = gaia_crsmtch_check(coords, width, chosen_frame, w, data, crop, Q)
                    if gaia_comp >= acc_comp_lvl:
                        print(f' Completion acceptable (>{acc_comp_lvl})! Moving on w/ catalog!\n')
                        break
                    else:
                        print(f' Gaia completion < {acc_comp_lvl}, defaulting to next catalog..\n')
            else:
                print(' No sources found in survey in this area!')

    if Q is None:
        raise ReferenceError('All catalogs failed to produce a result!')

    return Q, chosen_survey, mag_low_cutoff


# LOCAL SHIFT COMPLEX QUERY FUNCTION

"""def complex_query(raImage, decImage, band, boxsize, maglow=12, maghigh=14, bulge=False):
    # new automatic survey picking

    # current catalogs
    if bulge:
        # if bulge field, change to galactic coords for query
        print('Galactic bulge field detected!  Adjusting query parameters accordingly...')
        coords = SkyCoord(ra=raImage * u.degree, dec=decImage * u.degree, frame='fk5')
        coords = coords.galactic  # galactic conversion for bulge fields
        chosen_frame = 'galactic'
        frame_long = coords.l.deg
        frame_long_str = 'l = %.4f' % frame_long
        frame_lat = coords.b.deg
        frame_lat_str = 'b = %.4f' % frame_lat
        print('Converting coords to galactic: %s, %s' % (frame_long_str, frame_lat_str))

        mag_high_cutoff = maghigh
        mag_low_cutoff = maglow

    else:
        coords = SkyCoord(ra=raImage * u.degree, dec=decImage * u.degree, frame='fk5')
        chosen_frame = 'fk5'
        frame_long = raImage
        frame_long_str = 'RA: %.4f' % frame_long
        frame_lat = decImage
        frame_lat_str = 'DEC: %.4f' % frame_lat

        print('Non-bulge field, implementing wide mag range...')
        mag_high_cutoff = 16
        mag_low_cutoff = maglow

    width = boxsize

    not_null = '!=null'

    errbitoptions_2mass = ['!= null','~0??', '~?0?']     # 2mass errbits column, 1st is for J and 2nd is for H
    if band == 'J':
        errbit_2mass = errbitoptions_2mass[1]
    elif band == 'H':
        errbit_2mass = errbitoptions_2mass[2]
    else:
        errbit_2mass = errbitoptions_2mass[0]

    if bulge:
        acc_source_num = 200    # total number of sources allowed in full query before trying smaller box
    else:
        acc_source_num = 400

    # auto query

    # contains changes in bounds to iterate through if too many sources:
    # format: [boxsize multiplier, mag lim scalar change, errbits column constraint]
    bounds_change_list = [
        [1.0, 0, 0, '<=16', '!= null'],
        [1.0, 0.5, 0, '<=16', '!= null'],
        [1.0, 0.5, 0, '<16', errbit_2mass],
        [1.0, 0.5, 0.5, '<16', errbit_2mass],
        [1.0, 0.5, 0.75, '<16', errbit_2mass],
        [1.0, 0.5, 1.0, '<16', errbit_2mass],
        [1.0, 0.5, 1.25, '<16', errbit_2mass],
        [0.85, 0.5, 1.25, '<16', errbit_2mass],
    ]

    fctns_for_band = [fctn for fctn in PHOTOMETRY_QUERY_FUNCTIONS[band]]
    last_idx = fctns_for_band[-1]
    Q = None
    chosen_survey = None
    success_flag = False

    for bounds in bounds_change_list:
        boxscale, low_lim_scalar, high_lim_scalar, errbits_constraint, errbits_2M = bounds

        effective_boxsize = width * boxscale
        eff_mag_low_cutoff = mag_low_cutoff + low_lim_scalar
        eff_mag_high_cutoff = mag_high_cutoff - high_lim_scalar
        errbits = [errbits_constraint, errbits_2M]

        for fctn in fctns_for_band:

            print(
                '\nQuerying %s around %s, %s, boxwidth %.2f arcmin, '
                'mag lim of %s - %s, err constraints: %s, %s'
                % (fctn.__name__, frame_long_str, frame_lat_str,
                   effective_boxsize, eff_mag_low_cutoff, eff_mag_high_cutoff,
                   errbits_constraint, errbits_2M)
            )

            try:
                Q, chosen_survey = fctn(
                    coords=coords,
                    frame_long_str=frame_long_str,
                    frame_lat_str=frame_lat_str,
                    band=band,
                    width=effective_boxsize,  # pass tightened boxsize
                    mag_low_cutoff=eff_mag_low_cutoff,
                    mag_high_cutoff=eff_mag_high_cutoff,
                    chosen_frame=chosen_frame,
                    errbits=errbits,  # pass tightened errbits
                )
            except Exception as e:
                print(f'Error in query via {fctn.__name__}: {e}')
                continue

            if Q is None or len(Q[0]) == 0:
                print(f'No sources found via {fctn.__name__}, trying next survey...')
                continue

            print('Queried source total = ', len(Q[0]))

            # AB → VEGA correction applies regardless of bounds iteration
            if chosen_survey in ['DES_Y', 'DES_Z', 'Skymapper']:
                print('Adjusting AB mags to VEGA...')
                offset = AB_OFFSET_DICT.get(band, 0.0)
                eff_mag_high_cutoff -= offset  # adjust the local variable, not the outer one

            if len(Q[0]) <= acc_source_num:
                success_flag = True
                break  # good count — exit survey loop
            else:
                print(f'Too many sources (>{acc_source_num}), tightening bounds...')
                Q = None
                break  # exit survey loop, try next bounds tier

        if success_flag:
            break  # exit bounds loop entirely

    if Q is None:
        raise ReferenceError('All catalogs and bounds combinations failed to produce a result!')

    # success_flag = False
    # last_idx = catalogs[-1][1]
    # for chosen_survey, catNum in catalogs:
    #     for k in keycheck:
    #         if catNum in k:
    #             print('%s catalog found!' % k)
    #             vhs_table = result[''.join(k)]
    #             cols = vhs_table.colnames
    #             vhs_band_col = vhs_table[cols[2]]
    #
    #             if not np.all(vhs_band_col.mask):
    #                 print('Survey has coverage in %s band!' % band)
    #                 print('Survey = %s' % k)
    #
    #                 # AB surveys
    #                 if chosen_survey in ['DES_Y', 'DES_Z', 'Skymapper']:
    #                     print('Adjusting AB mags to VEGA...')
    #                     offset = AB_OFFSET_DICT.get(band, 0.0)
    #                     mag_high_cutoff -= offset
    #
    #                 no_sources_flag = False
    #
    #                 for bounds in bounds_change_list:
    #                     boxscale, low_lim_scalar, high_lim_scalar, errbits_constraint, errbits_2M = bounds
    #                     errbits = [errbits_constraint, errbits_2M]
    #
    #                     effective_boxsize = boxsize * boxscale
    #                     eff_mag_low_cutoff = mag_low_cutoff + low_lim_scalar
    #                     eff_mag_high_cutoff = mag_high_cutoff - high_lim_scalar
    #
    #                     print('\nQuerying Vizier %s around %s, %s, boxwidth %.2f arcmin, mag lim of %s - %s, '
    #                           'err constraints: %s, %s'
    #                           % (catNum, frame_long_str, frame_lat_str, effective_boxsize,
    #                              eff_mag_low_cutoff, eff_mag_high_cutoff, errbits_constraint, errbits_2M))
    #                     try:
    #                         if catNum == last_idx:
    #                             print('\nVizier catalogs exhausted, switching to local 2MASS query...')
    #                             Q = local_query_box(ra_center=raImage,
    #                                                 dec_center=decImage,
    #                                                 width=str(effective_boxsize) + 'm',
    #                                                 dbname='prime_2mass_local',
    #                                                 tablename='twomass_local',
    #                                                 columns=["ra", "dec", f"{band.lower()}mag", f"e_{band.lower()}mag"],
    #                                                 column_filters={
    #                                                     f"{band.lower()}mag": f"BETWEEN {eff_mag_low_cutoff:f} AND {eff_mag_high_cutoff:f}",
    #                                                     "cc_flg": errbits_2M,
    #                                                 }
    #                                                 )
    #                         else:
    #                             v = Vizier(columns=[cols[0], cols[1], cols[2]],
    #                                        column_filters={
    #                                            cols[2]: f"{eff_mag_low_cutoff:f}..{eff_mag_high_cutoff:f}",
    #                                            f"{band.lower()}Flag": "<4",
    #                                            f"{band}perrbits": errbits_constraint,
    #                                            f"{band}1perrb": errbits_constraint,
    #                                            f"{band}flags": errbits_constraint,
    #                                            "Cflg": errbits_2M,
    #                                            "Hclass": "== -1",
    #                                            "Class": "== 0",
    #                                            "Nd": ">6"
    #                                        }, row_limit=-1)
    #
    #                             Q = v.query_region(SkyCoord(coords, unit=(u.deg, u.deg)),
    #                                                width=str(effective_boxsize) + 'm',
    #                                                catalog=catNum, cache=False,
    #                                                frame=chosen_frame)
    #
    #                         if Q and len(Q[0]) > 0:
    #                             print('Queried source total = ', len(Q[0]))
    #                             if len(Q[0]) <= acc_source_num:
    #                                 success_flag = True  # Mark success
    #                                 break
    #                             else:
    #                                 print("Too many sources (>%i), trying different bounds..." % acc_source_num)
    #                                 no_sources_flag = True
    #                                 continue
    #                         else:
    #                             print(f"No sources found in {catNum}, trying fallback if available...")
    #                             break
    #
    #                     except Exception as e:
    #                         print('Error in Vizier query.')
    #                         print(f"Error details: {e}")
    #                         continue
    #
    #                 if success_flag:
    #                     break  # Break out of keycheck loop as well
    #
    #     if success_flag:
    #         break  # Break out of catalogs loop

    return Q, coords, catNum, cols[2], eff_mag_low_cutoff, eff_mag_high_cutoff, effective_boxsize, errbits"""

# ASTROMETRY QUERY FUNCTION (FOR SCAMP & SUCH)


def local_scamp_query(coords, width, chosen_frame):
    G, _ = gaia_query(coords=coords, width=width, chosen_frame=chosen_frame)
    gaia_data = G[0]
    # Mag error dummy column
    magerr = Column(data=np.ones(len(gaia_data[gaia_data.colnames[0]]))*2, name='mag_err', dtype=np.float32)
    gaia_data['mag_err'] = magerr

    # flags dummy column
    n = len(gaia_data[gaia_data.colnames[0]])
    flags_col = Column(np.zeros(n, dtype=np.int32), name='FLAGS', dtype=np.int32)
    gaia_data['FLAGS'] = flags_col

    # ra / dec error conversion, mas to degrees
    gaia_data['ra_error'] = gaia_data['ra_error'] / (60 * 60 * 1000)
    gaia_data['dec_error'] = gaia_data['dec_error'] / (60 * 60 * 1000)

    gaia_ldac_hdul = convert_table_to_ldac(gaia_data)

    # if 'FLAGS' not in gaia_data.colnames:
    #     n = len(gaia_data[gaia_data.colnames[0]])
    #     flags_col = Column(np.zeros(n, dtype=np.int32), name='FLAGS', dtype=np.int32)
    #     gaia_tbl = Table(gaia_data).copy()
    #     gaia_tbl.add_column(flags_col)
    # else:
    #     gaia_tbl = Table(gaia_data)
    #
    # doner = fits.open(donor_cat_path)
    # imhead = doner[1]
    #
    # hdu_imhead = fits.BinTableHDU(imhead.data, header=imhead.header)
    # hdu_imhead.header['EXTNAME'] = 'LDAC_IMHEAD'
    #
    # hdu_objects = fits.BinTableHDU(gaia_tbl)
    # hdu_objects.header['EXTNAME'] = 'LDAC_OBJECTS'
    #
    # primary = fits.PrimaryHDU()
    #
    # # Recreate LDAC format
    # hdul = fits.HDUList([primary,
    #                      hdu_imhead,
    #                      hdu_objects]
    # )

    ra = coords.ra.deg
    dec = coords.dec.deg

    gaia_cat_name = f'GAIA_CAT_{ra}_{dec}.ldac'

    gaia_ldac_hdul.writeto(gaia_cat_name, overwrite=True)

    with fits.open(f'GAIA_CAT_{ra}_{dec}.ldac', mode='update') as hdul:
        hdul.flush()

    return gaia_cat_name, gaia_data.colnames
