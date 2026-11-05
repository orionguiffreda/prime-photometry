"""
Utility functions for photometry
"""

import random
import os
import shutil
import numpy as np
from contextlib import contextmanager, redirect_stdout, redirect_stderr
import time
from astropy.io import fits
from filelock import FileLock, Timeout

from photometrus.settings import AB_OFFSET_DICT

from photometrus.utils.defaults import PROCESSING_DEFAULTS as defaults


# Context Managers
@contextmanager
def open_fits_robust(imageName, mode='update', timeout=30):
    """
    Provides robustness against file update errors from multiprocess

    Parameters
    ----------
    imageName : str
        Filename of input image
    mode : str
        Mode for fits.open (default 'update')
    timeout : int
        Max seconds to wait for the lock before raising (default 30)
    """

    lock_path = imageName + ".lock"
    lock = FileLock(lock_path, timeout=timeout)

    try:
        with lock:
            hdul = fits.open(imageName, mode=mode)
            try:
                yield hdul
            finally:
                hdul.close()
    except Timeout:
        raise Timeout(f"Could not acquire lock for '{imageName}' within {timeout}s")


@contextmanager
def log_output(enable, filename):
    """
    Option to log output to log file (cmd line prints go to log instead)

    Parameters
    ----------
    enable: bool
        To enable logging of cmd line outputs instead of just writing to terminal
    filename: str
        Filename for log file
    """

    if enable:
        with open(filename, "w") as f, redirect_stdout(f), redirect_stderr(f):
            yield
    else:
        yield


# UTILITY FUNCTIONS

# parallel running - filename switching for independent psf fit
def end_name_gen(ext='ecsv', magtype=defaults['magtype']):
    """
    Makes output file extension for plots & such depending on magtype

    Parameters
    ----------
    ext: str
        File extension
    magtype: str
        Optional, specify the magtype for pruning (default = AUTO)
    """
    if magtype:
        end_name = f'_{magtype}.{ext}'
    else:
        end_name = f'.{ext}'
    return end_name


# Read LDAC tables
def get_table_from_ldac(filename, frame=1):
    """
    Load an astropy table from a fits_ldac by frame (Since the ldac format has column
    info for odd tables, giving it twce as many tables as a regular fits BinTableHDU,
    match the frame of a table to its corresponding frame in the ldac file).

    Parameters
    ----------
    filename: str
        Name of the file to open
    frame: int
        Number of the frame in a regular fits file (default = 1)
    """

    from astropy.table import Table
    if frame > 0:
        frame = frame * 2
    tbl = Table.read(filename, hdu=frame)
    return tbl


def get_band(filepath):
    """
    Automatically retrieve band of image from header

    Parameters
    ----------
    filepath: str
        full filepath to fits image
    """

    band = fits.getval(filepath, 'FILTER2')
    if band == 'Open':
        band = 'Z'
    return band


def grb_rad_convert(rad):
    """
    Automatically convert input grb threshold to arcsec

    Parameters
    ----------
    rad: str
        Input grb threshold value, defaults to arcsec, but can specify in arcmin / deg using format: "rad=30_arcmin"
    """

    radsplit = rad.split('_')
    if len(radsplit) == 1 or radsplit[1] == 'arcsec':
        arcconvert = float(radsplit[0])
    elif radsplit[1] == 'arcmin':
        rad_f = float(radsplit[0])
        arcconvert = rad_f * 60
    elif radsplit[1] == 'deg':
        rad_f = float(radsplit[0])
        arcconvert = rad_f * 3600
    else:
        print('Only arcsec, arcmin, and deg are supported! Default = arcsec')
        raise Exception('Use supported units.')
    return float(abs(arcconvert))


def ab_convert(mag, band, survey=None, revert=False):
    """
    Converts Vega survey mags to AB using AB_OFFSET_DICT in settings.py

    Parameters
    ----------
    mag: array_like
        Input mag, usually column of values in photometry processing
    band: str
        Filter of survey
    survey: str
        Name of survey
    revert: bool
        Revert AB mags to Vega
    """

    mag = np.asarray(mag)
    # jy zero points
    # from 2MASS
    offset_dict = AB_OFFSET_DICT

    if survey == '2MASS':
        print(' 2MASS to AB')
        zp_dict = {'zp_J': 1594, 'zp_H': 1024}
        pick_zp = 'zp_' + band

        zp = zp_dict[pick_zp]

        flx = zp * 10 ** (-mag / 2.5)
        ab_mag = -2.5 * np.log10(flx / 3631)
        if revert:
            # print('Temporarily reverting PRIME mags to Vega to compare to survey (for GRB or plotting)')
            v_flx = 3631 * 10 ** (-mag / 2.5)
            vega_mag = -2.5 * np.log10(v_flx / zp)
            return vega_mag

    elif 'UKIDSS' in survey:
        print(' UKIDSS to AB')
        # for ukirt conversion: https://adsabs.harvard.edu/full/2006MNRAS.367..454H

        offset = offset_dict['UKIDSS'][band]

        ab_mag = mag+offset
        if revert:
            # print('Temporarily reverting PRIME mags to Vega to compare to survey (for GRB or plotting)')
            vega_mag = mag - offset
            return vega_mag

    elif survey in (
            "DES_Z", "DES_Y", "Skymapper", "SDSS",
            "PanSTARRS", "PanSTARRS_Z", "PRIME"
    ):
        print(' Survey %s is already reported in AB mag, no offset required.' % survey)
        ab_mag = mag

    else:
        # for vista (AB-Vega offsets): https://www.aanda.org/articles/aa/full_html/2015/03/aa24973-14/T3.html
        print(' VISTA to AB')

        offset = offset_dict['VISTA'][band]

        ab_mag = mag+offset
        if revert:
            # print('Temporarily reverting PRIME mags to Vega to compare to survey (for GRB or plotting)')
            vega_mag = mag - offset
            return vega_mag

    return ab_mag


def removal(directory, end_names=None):
    """
    Removes files from directory with specific suffixes

    Parameters
    ----------
    directory : str
        Directory to conduct removal processes on
    end_names: str
        Optional, list of suffixes to remove, using format: end_names=".lock,.psf", default = ".cat,.psf"
    """

    if end_names:
        end_names = end_names.split(',')
    else:
        end_names = ['.cat', '.psf']
    start_names = []    # ['GRB_', 'PSF.']
    preserved_end_names = []    # ['.ecsv', '.fits']
    rem_files = [f for f in os.listdir(directory) if f.endswith(tuple(end_names)) or f.startswith(tuple(start_names))
                 and not f.endswith(tuple(preserved_end_names))]
    for file in rem_files:
        path = os.path.join(directory + file)
        try:
            os.remove(path)
            # print(f"Removed file: {path}")
        except Exception as e:
            print(f"Error removing file: {path} - {e}")


def copy_given_cat_to_dir(given_catalog, name, directory, survey):
    given_cat_ecsvname = '%s.%s.ecsv' % (name, survey)
    new_cat_path = os.path.join(directory, given_cat_ecsvname)

    shutil.copy(given_catalog, new_cat_path)
