from astroquery.vsa import Vsa
from astropy.coordinates import Angle, SkyCoord
import astropy.units as u
import pandas as pd
import io
import re
import requests
import warnings
from astropy.table import Table
from astropy.table import join
from astropy.utils.metadata import MergeConflictWarning
from astropy.coordinates import Angle
from astropy.table import vstack
from datetime import datetime as dt
import numpy as np
from astropy.io.votable import parse_single_table
import psycopg2
import math
from astroquery.utils import TableList
import os
import argparse
import time
from psycopg2.extras import execute_values

from vista_cat_test import connection_details


# PURE SQL QUERY METHOD
warnings.simplefilter('ignore', MergeConflictWarning)

def vista_sql_query(ra, dec, rad, database, split_flag=False, vvv_col_set=None, init_flag=False):

    MAX_RETRIES = 5
    RETRY_DELAY = 30

    if init_flag:

        columns = ['sourceID']

        survey_dict = {
            'VHSDR7': ['110', columns],
            'VIKINGDR5': ['140', columns],
            'VVVDR5': ['120', columns]
        }

        survey_id, columns = survey_dict[database]

    else:

        base_columns = [
            'sourceID', 'ra', 'dec', 'l', 'b', 'priOrSec', 'mergedClassStat', 'mergedClass', 'pStar', 'pGalaxy', 'pNoise',
            'pSaturated', 'yPetroMag', 'yPetroMagErr', 'yPsfMag', 'yPsfMagErr', 'ySerMag2D', 'ySerMag2DErr', 'yAperMag3',
            'yAperMag3Err', 'yAperMag4', 'yAperMag4Err', 'yAperMag6', 'yAperMag6Err', 'yAperMagNoAperCorr3',
            'yAperMagNoAperCorr4', 'yAperMagNoAperCorr6', 'yEll', 'yPA', 'yErrBits', 'yClass', 'yClassStat', 'yppErrBits',
            'jPetroMag', 'jPetroMagErr', 'jPsfMag', 'jPsfMagErr', 'jSerMag2D', 'jSerMag2DErr', 'jAperMag3',
            'jAperMag3Err', 'jAperMag4', 'jAperMag4Err', 'jAperMag6', 'jAperMag6Err', 'jAperMagNoAperCorr3',
            'jAperMagNoAperCorr4', 'jAperMagNoAperCorr6', 'jEll', 'jPA', 'jErrBits', 'jClass', 'jClassStat', 'jppErrBits',
            'hPetroMag', 'hPetroMagErr', 'hPsfMag', 'hPsfMagErr', 'hSerMag2D', 'hSerMag2DErr', 'hAperMag3',
            'hAperMag3Err', 'hAperMag4', 'hAperMag4Err', 'hAperMag6', 'hAperMag6Err', 'hAperMagNoAperCorr3',
            'hAperMagNoAperCorr4', 'hAperMagNoAperCorr6', 'hEll', 'hPA', 'hErrBits', 'hClass', 'hClassStat', 'hppErrBits',
            'ksPetroMag', 'ksPetroMagErr', 'ksPsfMag', 'ksPsfMagErr', 'ksSerMag2D', 'ksSerMag2DErr', 'ksAperMag3',
            'ksAperMag3Err', 'ksAperMag4', 'ksAperMag4Err', 'ksAperMag6', 'ksAperMag6Err', 'ksAperMagNoAperCorr3',
            'ksAperMagNoAperCorr4', 'ksAperMagNoAperCorr6', 'ksEll', 'ksPA', 'ksErrBits', 'ksClass', 'ksClassStat',
            'ksppErrBits'
        ]

        vmc_columns = [
            'sourceID', 'ra', 'dec', 'l', 'b', 'priOrSec', 'mergedClassStat', 'mergedClass', 'pStar', 'pGalaxy', 'pNoise',
            'pSaturated', 'yPetroMag', 'yPetroMagErr', 'yPsfMag', 'yPsfMagErr', 'ySerMag2D', 'ySerMag2DErr', 'yAperMag3',
            'yAperMag3Err', 'yAperMag4', 'yAperMag4Err', 'yAperMag6', 'yAperMag6Err', 'yAperMagNoAperCorr3',
            'yAperMagNoAperCorr4', 'yAperMagNoAperCorr6', 'yEll', 'yPA', 'yErrBits', 'yClass', 'yClassStat', 'yppErrBits',
            'jPetroMag', 'jPetroMagErr', 'jPsfMag', 'jPsfMagErr', 'jSerMag2D', 'jSerMag2DErr', 'jAperMag3',
            'jAperMag3Err', 'jAperMag4', 'jAperMag4Err', 'jAperMag6', 'jAperMag6Err', 'jAperMagNoAperCorr3',
            'jAperMagNoAperCorr4', 'jAperMagNoAperCorr6', 'jEll', 'jPA', 'jErrBits', 'jClass', 'jClassStat', 'jppErrBits',
            'ksPetroMag', 'ksPetroMagErr', 'ksPsfMag', 'ksPsfMagErr', 'ksSerMag2D', 'ksSerMag2DErr', 'ksAperMag3',
            'ksAperMag3Err', 'ksAperMag4', 'ksAperMag4Err', 'ksAperMag6', 'ksAperMag6Err', 'ksAperMagNoAperCorr3',
            'ksAperMagNoAperCorr4', 'ksAperMagNoAperCorr6', 'ksEll', 'ksPA', 'ksErrBits', 'ksClass', 'ksClassStat',
            'ksppErrBits'
        ]

        vvv_columns = [
            'sourceID', 'ra', 'dec', 'l', 'b', 'priOrSec', 'mergedClassStat', 'mergedClass', 'pStar', 'pGalaxy', 'pNoise',
            'pSaturated',
            'zAperMag1', 'zAperMag1Err', 'zAperMag3', 'zAperMag3Err', 'zAperMag4', 'zAperMag4Err', 'zEll', 'zPA',
            'zErrBits', 'zClass', 'zClassStat', 'zppErrBits',
            'z_1AperMag1', 'z_1AperMag1Err', 'z_1AperMag3', 'z_1AperMag3Err', 'z_1AperMag4', 'z_1AperMag4Err', 'z_1Ell',
            'z_1PA',
            'z_1ErrBits', 'z_1Class', 'z_1ClassStat', 'z_1ppErrBits',
            'z_2AperMag1', 'z_2AperMag1Err', 'z_2AperMag3', 'z_2AperMag3Err', 'z_2AperMag4', 'z_2AperMag4Err', 'z_2Ell',
            'z_2PA',
            'z_2ErrBits', 'z_2Class', 'z_2ClassStat', 'z_2ppErrBits',
            'yAperMag1', 'yAperMag1Err', 'yAperMag3', 'yAperMag3Err', 'yAperMag4', 'yAperMag4Err', 'yEll', 'yPA',
            'yErrBits', 'yClass', 'yClassStat', 'yppErrBits',
            'y_1AperMag1', 'y_1AperMag1Err', 'y_1AperMag3', 'y_1AperMag3Err', 'y_1AperMag4', 'y_1AperMag4Err', 'y_1Ell',
            'y_1PA',
            'y_1ErrBits', 'y_1Class', 'y_1ClassStat', 'y_1ppErrBits',
            'y_2AperMag1', 'y_2AperMag1Err', 'y_2AperMag3', 'y_2AperMag3Err', 'y_2AperMag4', 'y_2AperMag4Err', 'y_2Ell',
            'y_2PA',
            'y_2ErrBits', 'y_2Class', 'y_2ClassStat', 'y_2ppErrBits',
            'jAperMag1', 'jAperMag1Err', 'jAperMag3', 'jAperMag3Err', 'jAperMag4', 'jAperMag4Err', 'jEll', 'jPA',
            'jErrBits', 'jClass', 'jClassStat', 'jppErrBits',
            'j_1AperMag1', 'j_1AperMag1Err', 'j_1AperMag3', 'j_1AperMag3Err', 'j_1AperMag4', 'j_1AperMag4Err', 'j_1Ell',
            'j_1PA',
            'j_1ErrBits', 'j_1Class', 'j_1ClassStat', 'j_1ppErrBits',
            'j_2AperMag1', 'j_2AperMag1Err', 'j_2AperMag3', 'j_2AperMag3Err', 'j_2AperMag4', 'j_2AperMag4Err', 'j_2Ell',
            'j_2PA',
            'j_2ErrBits', 'j_2Class', 'j_2ClassStat', 'j_2ppErrBits',
            'hAperMag1', 'hAperMag1Err', 'hAperMag3', 'hAperMag3Err', 'hAperMag4', 'hAperMag4Err', 'hEll', 'hPA',
            'hErrBits', 'hClass', 'hClassStat', 'hppErrBits',
            'h_1AperMag1', 'h_1AperMag1Err', 'h_1AperMag3', 'h_1AperMag3Err', 'h_1AperMag4', 'h_1AperMag4Err', 'h_1Ell',
            'h_1PA',
            'h_1ErrBits', 'h_1Class', 'h_1ClassStat', 'h_1ppErrBits',
            'h_2AperMag1', 'h_2AperMag1Err', 'h_2AperMag3', 'h_2AperMag3Err', 'h_2AperMag4', 'h_2AperMag4Err', 'h_2Ell',
            'h_2PA',
            'h_2ErrBits', 'h_2Class', 'h_2ClassStat', 'h_2ppErrBits',
            'ksAperMag1', 'ksAperMag1Err', 'ksAperMag3', 'ksAperMag3Err', 'ksAperMag4', 'ksAperMag4Err', 'ksEll', 'ksPA',
            'ksErrBits', 'ksClass', 'ksClassStat', 'ksppErrBits',
            'ks_1AperMag1', 'ks_1AperMag1Err', 'ks_1AperMag3', 'ks_1AperMag3Err', 'ks_1AperMag4', 'ks_1AperMag4Err',
            'ks_1Ell', 'ks_1PA',
            'ks_1ErrBits', 'ks_1Class', 'ks_1ClassStat', 'ks_1ppErrBits',
            'ks_2AperMag1', 'ks_2AperMag1Err', 'ks_2AperMag3', 'ks_2AperMag3Err', 'ks_2AperMag4', 'ks_2AperMag4Err',
            'ks_2Ell', 'ks_2PA',
            'ks_2ErrBits', 'ks_2Class', 'ks_2ClassStat', 'ks_2ppErrBits'
        ]

        vvv_columns_1 = [
            'sourceID', 'ra', 'dec', 'l', 'b', 'priOrSec', 'mergedClassStat', 'mergedClass', 'pStar', 'pGalaxy', 'pNoise',
            'pSaturated',
            'zAperMag1', 'zAperMag1Err', 'zAperMag3', 'zAperMag3Err', 'zAperMag4', 'zAperMag4Err', 'zEll', 'zPA',
            'zErrBits', 'zClass', 'zClassStat', 'zppErrBits',
            'z_1AperMag1', 'z_1AperMag1Err', 'z_1AperMag3', 'z_1AperMag3Err', 'z_1AperMag4', 'z_1AperMag4Err', 'z_1Ell',
            'z_1PA',
            'z_1ErrBits', 'z_1Class', 'z_1ClassStat', 'z_1ppErrBits',
            'z_2AperMag1', 'z_2AperMag1Err', 'z_2AperMag3', 'z_2AperMag3Err', 'z_2AperMag4', 'z_2AperMag4Err', 'z_2Ell',
            'z_2PA',
            'z_2ErrBits', 'z_2Class', 'z_2ClassStat', 'z_2ppErrBits',
            'yAperMag1', 'yAperMag1Err', 'yAperMag3', 'yAperMag3Err', 'yAperMag4', 'yAperMag4Err', 'yEll', 'yPA',
            'yErrBits', 'yClass', 'yClassStat', 'yppErrBits',
            'y_1AperMag1', 'y_1AperMag1Err', 'y_1AperMag3', 'y_1AperMag3Err', 'y_1AperMag4', 'y_1AperMag4Err', 'y_1Ell',
            'y_1PA',
            'y_1ErrBits', 'y_1Class', 'y_1ClassStat', 'y_1ppErrBits',
            'y_2AperMag1', 'y_2AperMag1Err', 'y_2AperMag3', 'y_2AperMag3Err', 'y_2AperMag4', 'y_2AperMag4Err', 'y_2Ell',
            'y_2PA',
            'y_2ErrBits', 'y_2Class', 'y_2ClassStat', 'y_2ppErrBits',
        ]

        vvv_columns_2 = [
            'sourceID',
            'jAperMag1', 'jAperMag1Err', 'jAperMag3', 'jAperMag3Err', 'jAperMag4', 'jAperMag4Err', 'jEll', 'jPA',
            'jErrBits', 'jClass', 'jClassStat', 'jppErrBits',
            'j_1AperMag1', 'j_1AperMag1Err', 'j_1AperMag3', 'j_1AperMag3Err', 'j_1AperMag4', 'j_1AperMag4Err', 'j_1Ell',
            'j_1PA',
            'j_1ErrBits', 'j_1Class', 'j_1ClassStat', 'j_1ppErrBits',
            'j_2AperMag1', 'j_2AperMag1Err', 'j_2AperMag3', 'j_2AperMag3Err', 'j_2AperMag4', 'j_2AperMag4Err', 'j_2Ell',
            'j_2PA',
            'j_2ErrBits', 'j_2Class', 'j_2ClassStat', 'j_2ppErrBits',
            'hAperMag1', 'hAperMag1Err', 'hAperMag3', 'hAperMag3Err', 'hAperMag4', 'hAperMag4Err', 'hEll', 'hPA',
            'hErrBits', 'hClass', 'hClassStat', 'hppErrBits',
            'h_1AperMag1', 'h_1AperMag1Err', 'h_1AperMag3', 'h_1AperMag3Err', 'h_1AperMag4', 'h_1AperMag4Err', 'h_1Ell',
            'h_1PA',
            'h_1ErrBits', 'h_1Class', 'h_1ClassStat', 'h_1ppErrBits',
            'h_2AperMag1', 'h_2AperMag1Err', 'h_2AperMag3', 'h_2AperMag3Err', 'h_2AperMag4', 'h_2AperMag4Err', 'h_2Ell',
            'h_2PA',
            'h_2ErrBits', 'h_2Class', 'h_2ClassStat', 'h_2ppErrBits',
            'ksAperMag1', 'ksAperMag1Err', 'ksAperMag3', 'ksAperMag3Err', 'ksAperMag4', 'ksAperMag4Err', 'ksEll', 'ksPA',
            'ksErrBits', 'ksClass', 'ksClassStat', 'ksppErrBits',
            'ks_1AperMag1', 'ks_1AperMag1Err', 'ks_1AperMag3', 'ks_1AperMag3Err', 'ks_1AperMag4', 'ks_1AperMag4Err',
            'ks_1Ell', 'ks_1PA',
            'ks_1ErrBits', 'ks_1Class', 'ks_1ClassStat', 'ks_1ppErrBits',
            'ks_2AperMag1', 'ks_2AperMag1Err', 'ks_2AperMag3', 'ks_2AperMag3Err', 'ks_2AperMag4', 'ks_2AperMag4Err',
            'ks_2Ell', 'ks_2PA',
            'ks_2ErrBits', 'ks_2Class', 'ks_2ClassStat', 'ks_2ppErrBits'
        ]

        vvv_columns_3 = [
            'sourceID',
            'jAperMag1', 'jAperMag1Err', 'jAperMag3', 'jAperMag3Err', 'jAperMag4', 'jAperMag4Err', 'jEll', 'jPA',
            'jErrBits', 'jClass', 'jClassStat', 'jppErrBits',
            'j_1AperMag1', 'j_1AperMag1Err', 'j_1AperMag3', 'j_1AperMag3Err', 'j_1AperMag4', 'j_1AperMag4Err', 'j_1Ell',
            'j_1PA',
            'j_1ErrBits', 'j_1Class', 'j_1ClassStat', 'j_1ppErrBits',
            'j_2AperMag1', 'j_2AperMag1Err', 'j_2AperMag3', 'j_2AperMag3Err', 'j_2AperMag4', 'j_2AperMag4Err', 'j_2Ell',
            'j_2PA',
            'j_2ErrBits', 'j_2Class', 'j_2ClassStat', 'j_2ppErrBits',
            'hAperMag1', 'hAperMag1Err', 'hAperMag3', 'hAperMag3Err', 'hAperMag4', 'hAperMag4Err', 'hEll', 'hPA',
            'hErrBits', 'hClass', 'hClassStat', 'hppErrBits',
            'h_1AperMag1', 'h_1AperMag1Err', 'h_1AperMag3', 'h_1AperMag3Err', 'h_1AperMag4', 'h_1AperMag4Err', 'h_1Ell',
            'h_1PA',
            'h_1ErrBits', 'h_1Class', 'h_1ClassStat', 'h_1ppErrBits',
        ]

        vvv_columns_4 = [
            'sourceID',
            'h_2AperMag1', 'h_2AperMag1Err', 'h_2AperMag3', 'h_2AperMag3Err', 'h_2AperMag4', 'h_2AperMag4Err', 'h_2Ell',
            'h_2PA',
            'h_2ErrBits', 'h_2Class', 'h_2ClassStat', 'h_2ppErrBits',
            'ksAperMag1', 'ksAperMag1Err', 'ksAperMag3', 'ksAperMag3Err', 'ksAperMag4', 'ksAperMag4Err', 'ksEll',
            'ksPA',
            'ksErrBits', 'ksClass', 'ksClassStat', 'ksppErrBits',
            'ks_1AperMag1', 'ks_1AperMag1Err', 'ks_1AperMag3', 'ks_1AperMag3Err', 'ks_1AperMag4', 'ks_1AperMag4Err',
            'ks_1Ell', 'ks_1PA',
            'ks_1ErrBits', 'ks_1Class', 'ks_1ClassStat', 'ks_1ppErrBits',
            'ks_2AperMag1', 'ks_2AperMag1Err', 'ks_2AperMag3', 'ks_2AperMag3Err', 'ks_2AperMag4', 'ks_2AperMag4Err',
            'ks_2Ell', 'ks_2PA',
            'ks_2ErrBits', 'ks_2Class', 'ks_2ClassStat', 'ks_2ppErrBits'
        ]

        if database == 'VVVDR5' and split_flag:

            survey_dict = {
                'VHSDR7': ['110', base_columns],
                'VIKINGDR5': ['140', base_columns],
                'VVVDR5': ['120', (vvv_columns_1, vvv_columns_2, vvv_columns_3, vvv_columns_4)],
            }

            try:
                survey_id, columns = survey_dict[database]
                columns = columns[vvv_col_set-1]
            except KeyError:
                raise KeyError('vvv_col_set not specified!  When using VVVDR5, specify a set of columns: 1, 2, 3 or 4')

        else:
            survey_dict = {
                'VHSDR7': ['110', base_columns],
                'VIKINGDR5': ['140', base_columns],
                'VVVDR5': ['120', vvv_columns],
                'VMCDR7': ['130', vmc_columns]
            }

            try:
                survey_id, columns = survey_dict[database]
            except KeyError:
                raise KeyError('Database not supported! Use supported db.')

    joined_cols = ','.join(columns)

    BASE_URL = "http://vsa.roe.ac.uk:8080/vdfs/WSASQL"

    params = {
        "programmeID": survey_id,
        "database": database,
        "archive": "VSA",
        "formaction": "region",
        "from": "source",
        "ra": ra,
        "dec": dec,
        "sys": "J",
        "name": "",
        "radius": rad,
        "xSize": "",
        "ySize": "",
        "boxAlignment": "RADec",
        "emailAddress": "",
        "format": "VOT",
        "compress": "",
        "rows": 30,
        "select": joined_cols,
        "where": "(priOrSec=0 OR priOrSec=frameSetID)",
    }

    response = requests.get(BASE_URL, params=params)
    response.raise_for_status()

    # print(response.text[:2000])

    for attempt in range(1, MAX_RETRIES + 1):
        response = requests.get(BASE_URL, params=params)
        response.raise_for_status()
        match = re.search(r'href="(http://vsa\.roe\.ac\.uk/tmp/tmp_sql/results[^"]+\.xml)"', response.text)
        if not match:
            print(f'  Attempt {attempt}/{MAX_RETRIES}: download link not found, retrying in {RETRY_DELAY}s...')
            time.sleep(RETRY_DELAY)
            continue

        csv_url = match.group(1)
        print(f" Found CSV URL: {csv_url}")
        csv_response = requests.get(csv_url)
        csv_response.raise_for_status()

        if not csv_response.content or len(csv_response.content) == 0:
            print(len(csv_response.content))
            # confirm empty content
            print(f'  Server returned 0 rows with no valid file, returning empty table.')
            return Table()

        # attempt to parse, only fall back if genuinely fails
        try:
            table = parse_single_table(io.BytesIO(csv_response.content)).to_table()
            return table
        except ValueError as e:
            if "Query returned 0 result rows" in response.text:
                print(f'  Server returned 0 rows with unparseable file, returning empty table.')
                return Table()
            print(f'  Attempt {attempt}/{MAX_RETRIES}: parse error ({e}), retrying in {RETRY_DELAY}s...')
            time.sleep(RETRY_DELAY)
            continue

    raise ValueError(f"Could not find VOTable download link after {MAX_RETRIES} attempts")


def ukirt_sql_query(ra, dec, rad, database, col_split_set=None):

    MAX_RETRIES = 5
    RETRY_DELAY = 30

    z_mag_cols = [
        'zHallMag', 'zHallMagErr', 'zPetroMag', 'zPetroMagErr', 'zPsfMag', 'zPsfMagErr', 'zSerMag2D', 'zSerMag2DErr', 'zAperMag3',
        'zAperMag3Err', 'zAperMag4', 'zAperMag4Err', 'zAperMag6', 'zAperMag6Err', 'zEll', 'zPA', 'zErrBits', 'zClass', 'zClassStat', 'zppErrBits',
    ]

    y_mag_cols = [
        'yHallMag', 'yHallMagErr', 'yPetroMag', 'yPetroMagErr', 'yPsfMag', 'yPsfMagErr', 'ySerMag2D', 'ySerMag2DErr', 'yAperMag3',
        'yAperMag3Err', 'yAperMag4', 'yAperMag4Err', 'yAperMag6', 'yAperMag6Err', 'yEll', 'yPA', 'yErrBits', 'yClass', 'yClassStat', 'yppErrBits',
    ]

    k_epoch_mag_cols = [
        'k_1HallMag', 'k_1HallMagErr', 'k_1PetroMag', 'k_1PetroMagErr', 'k_1PsfMag', 'k_1PsfMagErr', 'k_1SerMag2D', 'k_1SerMag2DErr', 'k_1AperMag3',
        'k_1AperMag3Err', 'k_1AperMag4', 'k_1AperMag4Err', 'k_1AperMag6', 'k_1AperMag6Err', 'k_1Ell', 'k_1PA', 'k_1ErrBits',
        'k_1Class', 'k_1ClassStat', 'k_1ppErrBits',

        'k_2HallMag', 'k_2HallMagErr', 'k_2PetroMag', 'k_2PetroMagErr', 'k_2PsfMag', 'k_2PsfMagErr', 'k_2SerMag2D',
        'k_2SerMag2DErr', 'k_2AperMag3',
        'k_2AperMag3Err', 'k_2AperMag4', 'k_2AperMag4Err', 'k_2AperMag6', 'k_2AperMag6Err', 'k_2Ell', 'k_2PA',
        'k_2ErrBits',
        'k_2Class', 'k_2ClassStat', 'k_2ppErrBits',
    ]

    j_epoch_mag_cols = [
        'j_1HallMag', 'j_1HallMagErr', 'j_1PetroMag', 'j_1PetroMagErr', 'j_1PsfMag', 'j_1PsfMagErr', 'j_1SerMag2D',
        'j_1SerMag2DErr', 'j_1AperMag3',
        'j_1AperMag3Err', 'j_1AperMag4', 'j_1AperMag4Err', 'j_1AperMag6', 'j_1AperMag6Err', 'j_1Ell', 'j_1PA',
        'j_1ErrBits', 'j_1Class',
        'j_1ClassStat', 'j_1ppErrBits',

        'j_2HallMag', 'j_2HallMagErr', 'j_2PetroMag', 'j_2PetroMagErr', 'j_2PsfMag', 'j_2PsfMagErr', 'j_2SerMag2D',
        'j_2SerMag2DErr', 'j_2AperMag3',
        'j_2AperMag3Err', 'j_2AperMag4', 'j_2AperMag4Err', 'j_2AperMag6', 'j_2AperMag6Err', 'j_2Ell', 'j_2PA',
        'j_2ErrBits', 'j_2Class',
        'j_2ClassStat', 'j_2ppErrBits',

    ]

    h2_mag_cols = [
        'h2HallMag', 'h2HallMagErr', 'h2PetroMag', 'h2PetroMagErr', 'h2PsfMag', 'h2PsfMagErr', 'h2SerMag2D', 'h2SerMag2DErr', 'h2AperMag3',
        'h2AperMag3Err', 'h2AperMag4', 'h2AperMag4Err', 'h2AperMag6', 'h2AperMag6Err', 'h2Ell', 'h2PA', 'h2ErrBits', 'h2Class', 'h2ClassStat', 'h2ppErrBits',
    ]

    columns = [
        'sourceID', 'ra', 'dec', 'l', 'b', 'priOrSec', 'mergedClassStat', 'mergedClass', 'pStar', 'pGalaxy', 'pNoise',
        'pSaturated',

        'jHallMag', 'jHallMagErr', 'jPetroMag', 'jPetroMagErr', 'jPsfMag', 'jPsfMagErr', 'jSerMag2D', 'jSerMag2DErr', 'jAperMag3',
        'jAperMag3Err', 'jAperMag4', 'jAperMag4Err', 'jAperMag6', 'jAperMag6Err', 'jEll', 'jPA', 'jErrBits', 'jClass', 'jClassStat', 'jppErrBits',

        'hHallMag', 'hHallMagErr', 'hPetroMag', 'hPetroMagErr', 'hPsfMag', 'hPsfMagErr', 'hSerMag2D', 'hSerMag2DErr', 'hAperMag3',
        'hAperMag3Err', 'hAperMag4', 'hAperMag4Err', 'hAperMag6', 'hAperMag6Err', 'hEll', 'hPA', 'hErrBits', 'hClass', 'hClassStat', 'hppErrBits',

        'kHallMag', 'kHallMagErr', 'kPetroMag', 'kPetroMagErr', 'kPsfMag', 'kPsfMagErr', 'kSerMag2D', 'kSerMag2DErr', 'kAperMag3',
        'kAperMag3Err', 'kAperMag4', 'kAperMag4Err', 'kAperMag6', 'kAperMag6Err', 'kEll', 'kPA', 'kErrBits', 'kClass', 'kClassStat', 'kppErrBits'
    ]

    chosen_database = database

    if database == 'LASDR11':
        chosen_database = 'UKIDSSDR11PLUS'
        y_insert_idx = columns.index('pSaturated') + 1
        columns[y_insert_idx:y_insert_idx] = y_mag_cols
        j_epoch_insert_idx = columns.index('jppErrBits') + 1
        columns[j_epoch_insert_idx:j_epoch_insert_idx] = j_epoch_mag_cols

    elif database == 'GPSDR11':
        chosen_database = 'UKIDSSDR11PLUS'
        k_insert_idx = columns.index('kppErrBits') + 1
        columns[k_insert_idx:k_insert_idx] = k_epoch_mag_cols
        h_insert_idx = columns.index('k_2ppErrBits') + 1
        columns[h_insert_idx:h_insert_idx] = h2_mag_cols
        columns = [col for col in columns if 'Mag6' not in col]

    elif database == 'GCSDR11':
        chosen_database = 'UKIDSSDR11PLUS'
        z_insert_idx = columns.index('pSaturated') + 1
        columns[z_insert_idx:z_insert_idx] = z_mag_cols
        y_insert_idx = columns.index('zppErrBits') + 1
        columns[y_insert_idx:y_insert_idx] = y_mag_cols
        k_insert_idx = columns.index('kppErrBits') + 1
        columns[k_insert_idx:k_insert_idx] = k_epoch_mag_cols

    if col_split_set == 1:
        columns = columns[:int(len(columns) / 2)]
    elif col_split_set == 2:
        columns = [columns[0]] + columns[int(len(columns) / 2):]

    joined_cols = ','.join(columns)

    survey_dict = {
        'LASDR11': ['101', columns],
        'UHSDR3': ['107', columns],
        'GPSDR11': ['102', columns],
        'GCSDR11': ['103', columns],
    }

    try:
        survey_id, columns = survey_dict[database]
    except KeyError:
        raise KeyError('Database not supported! Use supported db.')

    BASE_URL = "http://wsa.roe.ac.uk:8080/wsa/WSASQL"

    params = {
        "database": chosen_database,
        "programmeID": survey_id,
        "from": "source",
        "formaction": "region",
        "ra": ra,
        "dec": dec,
        "sys": "J",
        "name": "",
        "radius": rad,
        "xSize": "",
        "ySize": "",
        "boxAlignment": "RADec",
        "emailAddress": "",
        "format": "VOT",
        "compress": "",
        "rows": 30,
        "select": joined_cols,
        "where": "(priOrSec=0 OR priOrSec=frameSetID)",
    }

    response = requests.get(BASE_URL, params=params)
    response.raise_for_status()

    # print(response.text[:2000])

    for attempt in range(1, MAX_RETRIES + 1):
        response = requests.get(BASE_URL, params=params)
        response.raise_for_status()
        match = re.search(r'href="(http://wsa\.roe\.ac\.uk/tmp/tmp_sql/results[^"]+\.xml)"', response.text)
        if not match:
            print(f'  Attempt {attempt}/{MAX_RETRIES}: download link not found, retrying in {RETRY_DELAY}s...')
            # time.sleep(RETRY_DELAY)
            # continue

        csv_url = match.group(1)
        print(f" Found CSV URL: {csv_url}")
        csv_response = requests.get(csv_url)
        csv_response.raise_for_status()

        if not csv_response.content or len(csv_response.content) == 0:
            print(len(csv_response.content))
            # confirm empty content
            print(f'  Server returned 0 rows with no valid file, returning empty table.')
            return Table()

        # attempt to parse, only fall back if genuinely fails
        try:
            table = parse_single_table(io.BytesIO(csv_response.content)).to_table()
            return table
        except ValueError as e:
            if "Query returned 0 result rows" in response.text:
                print(f'  Server returned 0 rows with unparseable file, returning empty table.')
                return Table()
            print(f'  Attempt {attempt}/{MAX_RETRIES}: parse error ({e}), retrying in {RETRY_DELAY}s...')
            time.sleep(RETRY_DELAY)
            continue

    raise ValueError(f"Could not find VOTable download link after {MAX_RETRIES} attempts")


def adjust_sql_table_and_export(table, storage_dir, database, id):
    # Match astroquery dtypes based on column name patterns
    for col in table.colnames:
        # int16 columns
        if col in ('mergedClass', 'yClass', 'jClass', 'j_1Class', 'j_2Class', 'hClass', 'ksClass', 'kClass', 'zClass',
                   'k_1Class', 'k_2Class', 'h2Class'):
            table[col] = table[col].astype(np.int16)
        # int32 columns
        elif col.endswith('ErrBits') or col.endswith('ppErrBits'):
            table[col] = table[col].astype(np.int32)
        # int64 columns (keep as-is)
        elif col in ('sourceID', 'priOrSec'):
            table[col] = table[col].astype(np.int64)
        # float64 columns (keep as-is)
        elif col in ('ra', 'dec', 'l', 'b', 'distance'):
            table[col] = table[col].astype(np.float64)
        # everything else → float32
        else:
            table[col] = table[col].astype(np.float32)

    surv_basename = database.split('DR')[0]
    file_name = os.path.join(storage_dir, f'{surv_basename}_{id}.fits.gz')
    table.write(file_name, format='fits', overwrite=True)
    print(f' Written: {file_name}')


def query_vista(storage_dir, ra, dec, rad, id, database):

    ROW_CAP = 2700000

    # rad in arcmin
    start_time = dt.now()

    if database == 'VVVDR5':
        # init_table = vista_sql_query(ra=ra, dec=dec, rad=rad, database=database, init_flag=True)
        # print(f' Init Source #: {len(init_table)}')
        #
        # if len(init_table) > ROW_CAP:

        table1 = vista_sql_query(ra=ra, dec=dec, rad=rad, database=database, split_flag=True, vvv_col_set=1)
        print(f' Table 1 Source #: {len(table1)}')

        if len(table1) == 0:
            print(' No sources found!')
            table = table1

        elif len(table1) > ROW_CAP:
            print(f' Row cap ({ROW_CAP}) reached! Splitting 2nd set of cols into 2...')

            table3 = vista_sql_query(ra=ra, dec=dec, rad=rad, database=database, split_flag=True, vvv_col_set=3)
            print(f' Table 3 Source #: {len(table3)}')
            table4 = vista_sql_query(ra=ra, dec=dec, rad=rad, database=database, split_flag=True, vvv_col_set=4)
            print(f' Table 4 Source #: {len(table4)}')

            print(' Joining 2nd col set...')
            del table3['distance']
            del table4['distance']

            table2 = join(table3, table4, keys='sourceID')

            print(' Joining tables...')
            del table1['distance']
            table = join(table1, table2, keys='sourceID')

        else:
            table2 = vista_sql_query(ra=ra, dec=dec, rad=rad, database=database, split_flag=True, vvv_col_set=2)
            print(f' Table 2 Source #: {len(table2)}')

            print(' Joining tables...')
            del table1['distance']
            table = join(table1, table2, keys='sourceID')

        # else:
        #     print(f'Under row cap ({ROW_CAP}), continuing with all cols...')
        #     table = vista_sql_query(ra=ra, dec=dec, rad=rad, database=database)
        #     print(f' Source #: {len(table)}')

    else:
        table = vista_sql_query(ra=ra, dec=dec, rad=rad, database=database)
        print(f' Source #: {len(table)}')

    adjust_sql_table_and_export(table=table, storage_dir=storage_dir, database=database, id=id)
    end_time = dt.now()
    print(f' Full SQL query & fits write time: {(end_time - start_time).total_seconds()}\n')


def query_ukirt(storage_dir, ra, dec, rad, id, database):

    # ROW_CAP = 3720000

    # rad in arcmin
    start_time = dt.now()

    if database == 'GPSDR11':
        table1 = ukirt_sql_query(ra=ra, dec=dec, rad=rad, database=database, col_split_set=1)
        print(f' Tbl 1 Source #: {len(table1)}')

        if len(table1) == 0:
            print(' No sources found!')
            table = table1

        else:
            table2 = ukirt_sql_query(ra=ra, dec=dec, rad=rad, database=database, col_split_set=2)
            print(f' Tbl 2 Source #: {len(table2)}')

            print(' Joining tables...')
            del table1['distance']
            table = join(table1, table2, keys='sourceID')
            table[:10].pprint_all()

    else:
        table = ukirt_sql_query(ra=ra, dec=dec, rad=rad, database=database)
        print(f' Source #: {len(table)}')

    adjust_sql_table_and_export(table=table, storage_dir=storage_dir, database=database, id=id)

    end_time = dt.now()
    print(f' Full SQL query & fits write time: {(end_time - start_time).total_seconds()}\n')
#%%
# query_ukirt('/mnt/photometry/local_catalog_db_test/ukirt_test/', '05:25:37.77', '20:52:09.5', 90, 'test_query',
#             database='GCSDR11' )
#%% VISTA DB LOADER


def load_to_db(red_data, dbname, tablename):
    start_time = dt.now()
    conn = psycopg2.connect(
        dbname=dbname,
        user=connection_details['user'],
        password=connection_details['password'],
        host=connection_details['host'],
        port=connection_details['port']
    )

    print(f' Row count: {len(red_data)}')
    # Get column names from the SQL table
    cur = conn.cursor()
    cur.execute("""
        SELECT column_name 
        FROM information_schema.columns 
        WHERE table_name = %s
        ORDER BY ordinal_position
    """, (tablename,))

    sql_columns = [row[0] for row in cur.fetchall()]
    sql_columns_str = ', '.join(sql_columns)

    int16_cols = [
        "k_1Class",
        "k_2Class",
        "h2Class",
    ]

    print(' Converting table to dataframe')
    df = red_data.to_pandas(use_nullable_int=True)

    for col in int16_cols:
        if col in df.columns:
            df[col] = (
                pd.to_numeric(df[col], errors="coerce")
                .astype("Int16")  # nullable integer dtype
            )

    df = df.where(pd.notnull(df), None)

    # for col in df.columns:
    #     if df[col].dtype == np.float64 or df[col].dtype == np.float32:
    #         df[col] = df[col].astype(object)
    #     elif df[col].dtype == np.int64 or df[col].dtype == np.int32 or df[col].dtype == np.int16:
    #         df[col] = df[col].astype(object)

    conv_time = dt.now()
    print(' Conversion time:', (conv_time - start_time).total_seconds())

    staging_table = f"staging_{tablename}"

    cur.execute(f"""
        CREATE TEMP TABLE IF NOT EXISTS {staging_table}
        (LIKE {tablename} INCLUDING DEFAULTS)
    """)

    cur.execute(f"TRUNCATE {staging_table}")

    buffer = io.StringIO()

    df.to_csv(
        buffer,
        index=False,
        header=False,
        sep='\t',
        na_rep='\\N'
    )
    buffer.seek(0)

    cur.copy_from(
        buffer,
        staging_table,
        sep='\t',
        null='\\N',
        columns=sql_columns
    )

    insert_query = f"""
        INSERT INTO {tablename} ({sql_columns_str})
        SELECT {sql_columns_str}
        FROM {staging_table}
        ON CONFLICT (sourceID) DO NOTHING
    """

    cur.execute(insert_query)
    conn.commit()
    conn.close()

    end_time = dt.now()
    print('Full catalog to database loading time:', (end_time - start_time).total_seconds())


# sql_query = Table.read('/mnt/fits_data/PRIME_local_catalogs/vista/vmc/fields_0_999/VMC_field971.fits.gz')
# if len(sql_query) == 0:
#     print('No data found in file!  Must not be covered by survey!  Skipping loading process.')
# else:
#     load_to_db(red_data=sql_query, dbname='prime_2mass_local', tablename='vmc_sources_test')

#%% BOX QUERY
from photometrus.query import ang_convert

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
        # use radial query near poles
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

    # conn = psycopg2.connect(service=dbname)

    conn = psycopg2.connect(
        dbname="prime_2mass_local",
        user=connection_details['user'],
        password=connection_details['password'],
        host=connection_details['host'],
        port=connection_details['port']
    )

    # Get column names from the SQL table
    cur = conn.cursor()
    cur.execute("""
        SELECT column_name 
        FROM information_schema.columns 
        WHERE table_name = %s
        ORDER BY ordinal_position
    """, (tablename,))

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

# Q = local_query_box(74.91, -69.81, 'prime_2mass_local', 'vmc_sources_test', width='40m')
#
# Q[0][:10].pprint_all()
# print(len(Q[0]))

#%% NO SOURCE FILE REMOVER

# folder_path = '/mnt/fits_data/PRIME_local_catalogs/vista/vhs/fields_21000_21999/'
# size_threshold = 200 * 1024     # 2kb threshold
#
# for filename in os.listdir(folder_path):
#     file_path = os.path.join(folder_path, filename)
#
#     if os.path.isfile(file_path):
#         if os.path.getsize(file_path) < size_threshold:
#             print(f"Removing {filename} ({os.path.getsize(file_path)} bytes)")
#             os.remove(file_path)


#%% VHS MASS DOWNLOADING
# storage_dir = '/mnt/fits_data/PRIME_local_catalogs/vista/vhs/'


def mass_vista_download(storage_directory=None,
                grid_file='/home/alex/PycharmProjects/prime-photometry/photometrus/configs/obsable_all_sky_grid.csv',
                field_range=None,
                database=None
):

    ukirt_surveys = [
        'LASDR11', 'UHSDR3', 'GPSDR11', 'GCSDR11'
    ]

    if not storage_directory:
        raise Exception('Specify storage directory!')

    df = pd.read_csv(grid_file, comment="#")
    Field_arr = df['ObjectName']
    RA_arr = df['RA']
    DEC_arr = df['DEC']

    surv_basename = database.split('DR')[0]

    def get_field_number(id):
        # extracts field num
        return int(re.search(r'\d+', id).group())

    def field_exists(id, storage_directory):
        # determine if field w/ id already exists in subdirs
        num = get_field_number(id)
        subdir_start = (num // 1000) * 1000
        subdir_end = subdir_start + 999
        subdir = os.path.join(storage_directory, f'fields_{subdir_start}_{subdir_end}')
        return os.path.exists(os.path.join(subdir, f'{surv_basename}_{id}.fits.gz'))

    unprocessed = [(id, ra, dec) for id, ra, dec in zip(Field_arr, RA_arr, DEC_arr)
        if not field_exists(id, storage_directory)
    ]

    total = len(Field_arr)

    # filter by range if specified
    if field_range is not None:
        range_start, range_end = field_range
        print(f'Field range specified: {range_start} - {range_end}')
        unprocessed = [(id, ra, dec) for id, ra, dec in unprocessed
                       if range_start <= get_field_number(id) <= range_end]
        total = abs(range_start - range_end) + 1

    remaining = len(unprocessed)
    print(f"Fields remaining: {remaining}/{total} ({total - remaining} already processed)")

    for id, pointing_ra, pointing_dec in unprocessed:
        # determine subdirectory based on field id
        num = get_field_number(id)
        subdir_start = (num // 1000) * 1000
        subdir_end = subdir_start + 999
        subdir = os.path.join(storage_directory, f'fields_{subdir_start}_{subdir_end}')
        os.makedirs(subdir, exist_ok=True)

        print(f' Querying {surv_basename} for grid #: {id}, RA = {pointing_ra}, DEC = {pointing_dec}')
        if database in ukirt_surveys:
            query_ukirt(storage_dir=subdir, ra=pointing_ra, dec=pointing_dec, rad=55, id=id, database=database)
        else:
            query_vista(storage_dir=subdir, ra=pointing_ra, dec=pointing_dec, rad=55, id=id, database=database)


# mass_vhs_download(storage_directory=storage_dir)

#%%

def mass_vista_db_loader(storage_directory='/mnt/fits_data/PRIME_local_catalogs/vista/vmc/',
                       dbname=None, tablename=None, field_range=None
):

    if not tablename:
        raise Exception(f'Tablename not specified!, database = {dbname}')
    else:
        print(f'DBNAME = {dbname}, TBLNAME = {tablename}\n')

    def get_field_number(id):
        return int(re.search(r'\d+', id).group())

    # determine which subdirectories to scan based on field_range
    if field_range is not None:
        range_start, range_end = field_range

        start_block = (range_start // 1000) * 1000
        end_block = (range_end // 1000) * 1000

        subdirs = [
            os.path.join(storage_directory, f'fields_{s}_{s + 999}')
            for s in range(start_block, end_block + 1, 1000)
        ]
    else:
        # subdirectory scan
        subdirs = sorted([
            os.path.join(storage_directory, d)
            for d in os.listdir(storage_directory)
            if os.path.isdir(os.path.join(storage_directory, d)) and d.startswith('fields_')
        ])

    # get all fits.gz files
    fits_files = []
    for subdir in subdirs:
        if not os.path.exists(subdir):
            print(f'Subdirectory not found, skipping: {subdir}')
            continue
        for fname in sorted(os.listdir(subdir)):
            if not fname.endswith('.fits.gz'):
                continue
            # check field num
            num = get_field_number(fname)
            if field_range is not None and not (range_start <= num <= range_end):
                continue
            fits_files.append(os.path.join(subdir, fname))

    total = len(fits_files)
    print(f'Found {total} fits files to load')

    sorted_fits_files = (sorted(fits_files, key=lambda x: int(re.split(r'[._]', x)[-3][5:])))
    for i, fpath in enumerate(sorted_fits_files, 1):
        fname = os.path.basename(fpath)
        print(f'\n[{i}/{total}] Loading {fname}')
        print(fpath)
        try:
            data = Table.read(fpath)
            if len(data) == 0:
                print('No data found in file!  Must not be covered by survey!  Skipping loading process.')
            else:
                load_to_db(red_data=data, dbname=dbname, tablename=tablename)
            # print(f'load_to_db({fpath}, {dbname}, {tablename})')
        except Exception as e:
            print(f'  ERROR loading {fname}: {e}')
            continue

    print(f'\nDatabase loading complete: {total} files processed')


# mass_vhs_db_loader(field_range=(995,2995))

#%%

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-db_load', action='store_true', help='optional arg to load already downloaded files'
                                                              'into db')
    parser.add_argument('-download', action='store_true', help='optional arg to download files via sql querying')
    parser.add_argument('--range', type=int, nargs=2, metavar=('START', 'END'))
    parser.add_argument('--database', type=str, metavar='DB', help='Used during download, specify which '
                                                                   'VISTA / UKIRT db you want to query')
    parser.add_argument('--storage_directory', type=str, metavar='SD', help='Specify'
                                                                            ' storage directory where downloaded .fits.gz'
                                                                            ' files are stored.')
    parser.add_argument('--dbname', type=str, metavar='DBNAME', help='Used during db loading, specify name of '
                                                                     'psql database youre loading into, '
                                                                     'default = prime_vhs_local',
                        default='prime_vhs_local')
    parser.add_argument('--tablename', type=str, metavar='TBLNAME', help='Used during db loading, specify name of '
                                                                     'psql database table youre loading into, ex. vhs_sources')
    args = parser.parse_args()

    field_range = tuple(args.range) if args.range else None
    if args.download:
        mass_vista_download(field_range=field_range, database=args.database, storage_directory=args.storage_directory)
    elif args.db_load:
        mass_vista_db_loader(field_range=field_range, dbname=args.dbname, tablename=args.tablename,
                           storage_directory=args.storage_directory)
    else:
        raise Exception('specify -download or -db_load for functionality!')


if __name__ == "__main__":
    main()
