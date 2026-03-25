"""
Utility functions for photometry
"""

import random
import numpy as np
from contextlib import contextmanager, redirect_stdout, redirect_stderr
import time
from astropy.io import fits

from photometrus.settings import AB_OFFSET_DICT

from photometrus.utils.defaults import PROCESSING_DEFAULTS as defaults


# Context Managers
@contextmanager
def open_fits_robust(imageName, mode='update', retries=3, delay=1, jitter=1):
    """
    Provides robustness against file update errors from multiprocess

    Parameters
    ----------
    imageName: str
        Filename of input image
    mode: str
        Mode for fits to open file in
    retries: int
        Max retries to open fits file (default = 3)
    delay: int
        Time (s) to delay before attempting to reopen the file (default = 1)
    jitter: int
        Random time betw. 0 and input added to delay (default = 1)
    """

    hdul = None
    for attempt in range(retries):
        try:
            hdul = fits.open(imageName, mode=mode)
            break
        except (FileNotFoundError, OSError) as e:
            if attempt < retries - 1:
                sleep_time = delay + random.uniform(0, jitter)
                print(f"FITS open failed (process collision?), trying again w/ delay: {round(sleep_time,2)}s")
                time.sleep(sleep_time)
            else:
                raise
    try:
        yield hdul
    finally:
        if hdul is not None:
            for attempt in range(retries):
                try:
                    hdul.close()
                    break
                except (FileNotFoundError, OSError) as e:
                    if attempt < retries - 1:
                        sleep_time = delay + random.uniform(0, jitter)
                        print(f"FITS close failed (process collision?), trying again w/ delay: {round(sleep_time,2)}s")
                        time.sleep(sleep_time)
                    else:
                        print(f"*WARNING*: FITS close/flush failed after {retries}")


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
