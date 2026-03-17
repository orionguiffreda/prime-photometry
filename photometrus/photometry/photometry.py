"""
Calibrates photometry for stacked image
"""

import os
import sys
import re
import random
from multiprocessing import Process, Barrier, Event, Queue
from queue import Empty
from contextlib import contextmanager, redirect_stdout, redirect_stderr

import astropy.nddata.utils
import numpy as np
import numpy.ma as ma
import argparse
import pandas as pd
import astropy.units as u
import time
from astroquery.vizier import Vizier
from astroquery.ipac.ned import Ned
from astropy.coordinates import Angle, SkyCoord
from astropy.wcs import WCS
from astropy.wcs import utils
from astropy.stats import sigma_clip, sigma_clipped_stats
from astropy.io import fits
from astropy.io import ascii
from astropy.nddata import Cutout2D
from astropy.io.fits import ImageHDU
from astropy.utils.data import Conf
from bs4 import BeautifulSoup
from regions import CircleSkyRegion
import base64
from astropy.table import Column, MaskedColumn
from astropy.table import Table
import matplotlib.pyplot as plt
import statsmodels.api as sm
import subprocess
from scipy.stats import skew
from scipy import odr
import warnings
from datetime import datetime as dt

from photometrus.settings import (gen_config_file_name, bulge_checker, PHOTOMETRY_MAG_LOWER_LIMIT, PHOTOMETRY_MAG_UPPER_LIMIT,
                                  PHOTOMETRY_QUERY_WIDTH, PHOTOMETRY_QUERY_CATALOGS, PHOTOMETRY_LIM_MAGS,
                                  AB_OFFSET_DICT, PRIME_FILTERS_DICT, get_weight_thresh, set_vizier_mirror, CHIP_ZPS,
                                  MAGTYPES, local_query_box)

from photometrus.utils.defaults import PROCESSING_DEFAULTS as defaults

# %%
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings(action="ignore", module="scipy", message="^One or more")
warnings.filterwarnings(action="ignore", module="numpy", message="Warning: 'partition' will ignore the 'mask' of the MaskedColumn.")

magtype = defaults['magtype']
parallel = defaults['parallel']

# Context Managers
@contextmanager
def open_fits_robust(imageName, mode='update', retries=3, delay=1, jitter=1):
    # provides robustness against file update errors from multiprocess
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


# option to log output to log file (prints go to log instead)
@contextmanager
def log_output(enable, filename):
    if enable:
        with open(filename, "w") as f, redirect_stdout(f), redirect_stderr(f):
            yield
    else:
        yield


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
        Number of the frame in a regular fits file
    """
    from astropy.table import Table
    if frame > 0:
        frame = frame * 2
    tbl = Table.read(filename, hdu=frame)
    return tbl


#%% vega to AB mag conversion
def ab_convert(mag, band, survey=None, revert=False):
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


# efficiency calculation
def gen_efficiency(imageName, PSFSources, band, specific_magtype=None):
    """
    Calculates median efficiency for stacked image and adds to image header

    Parameters
    ----------
    imageName: str
        Filename of stacked image
    PSFSources: str
        Full PRIME source catalog for image
    band: str
        Filter of observation
    """

    coadd_hdr = fits.getheader(imageName)

    if specific_magtype:
        conv_facs = [kword for kword in coadd_hdr if f'FAC_{MAGTYPES[specific_magtype]}' in kword]
    else:
        conv_facs = [kword for kword in coadd_hdr if f'FAC_{MAGTYPES[magtype]}' in kword]

    efficiencies = []

    if conv_facs:
        for conv_fac_name in conv_facs:

            mag_type_name = conv_fac_name.replace('FAC_', '')

            conv_fac = coadd_hdr[conv_fac_name] * (u.microjansky / u.adu)

            # all variables
            diam_tel = 1.8 * u.m
            diam_obstruct = 0.8 * u.m
            exp_time = 2.86 * u.s
            gain = 1.8 * (u.electron / u.adu)
            c = 3*10**8 * (u.m / u.s) # m/s
            h = 6.626*10**-34 * (u.J * u.s)  # J*s
            band_lambda = PRIME_FILTERS_DICT[band]['lambda'] * u.nm    # m
            band_fwhm = PRIME_FILTERS_DICT[band]['fwhm'] * u.nm     # m
            band_freq = (c / band_lambda).to(u.Hz)   # hz

            # telescope area
            a_tel = np.pi * (diam_tel / 2) ** 2 - np.pi * (diam_obstruct / 2) ** 2 # area of telescope, m^2

            eta = (gain * h * band_lambda) / (conv_fac * exp_time * band_fwhm * a_tel)
            eta = eta.decompose().value
            # print('throughput = ',eta)

            # fluxes
            if 'PSF' in conv_fac_name:
                flx_type = 'POINTSOURCE'
            else:
                flx_type = MAGTYPES[magtype]

            try:
                mag_flxs = PSFSources[f'{MAGTYPES[magtype]}_FLUX_DENSITY']          # flux calculated from mag col
                meas_flxs = PSFSources[f'FLUX_{flx_type}']     # measured flux
            except KeyError:
                mag_flxs = meas_flxs = 0

            if not isinstance(mag_flxs, int):
                si_mag_flxs = mag_flxs.to(u.W / (u.m**2 * u.Hz))

                # bandwidth
                # J_bandwidth = c / (J_lambda - 0.5 * J_fwhm) - c / (J_lambda + 0.5 * J_fwhm)     # hz
                bandwidth = ((c * band_fwhm) / band_lambda ** 2).to(u.Hz)   # hz

                # counts
                J_n_photons = (si_mag_flxs * a_tel * bandwidth) / (h * band_freq)    # num of photons
                # PSFSources['Photon Rate'] = J_n_photons.to(1 / u.s)

                meas_flxs = meas_flxs * (u.adu / u.ct)
                J_n_elec = meas_flxs * gain / exp_time
                # PSFSources['Electron Rate'] = J_n_elec.to(u.electron / u.s)

                efficiency = J_n_elec / J_n_photons
                # PSFSources['Efficiency'] = efficiency.decompose().value

                med_efficiency = round(np.nanmedian(efficiency.decompose().value), 5)
                efficiencies.append(med_efficiency)

                with open_fits_robust(imageName) as hdul:
                    hdr = hdul[0].header
                    hdr.set(f'EFF_{mag_type_name}', med_efficiency, 'e- rate / exp. photon rate',
                                 after=conv_fac_name)

                print(f'Med {mag_type_name} efficiency :', med_efficiency)

        for eff in efficiencies:
            if eff > 0.7:
                print(f'\n*WARNING* Med efficiency larger than expected!: {eff}, investigation '
                                 f'required!  Was the catalog AB conversion applied multiple times?\n')
        return efficiencies
    else:
        print('No applicable conversion factors found!  Cannot calculate efficiency!')


# convert ab mag to microjansky
def ab_to_microjy(m_ab):
    """Convert AB magnitude(s) to microjanskys
    """
    m_ab = m_ab.value
    # f_microjy = (3631 * u.Jy * 10**(-0.4 * m_ab)).to(u.microjansky)
    f_jy = (10 ** (23 - (m_ab + 48.6) / 2.5)) * u.jansky
    f_microjy = f_jy.to(u.microjansky)
    return f_microjy

#%%
# import img and get wcs


def img(directory, imageName, crop):
    os.chdir(directory)
    f = fits.open(imageName)
    data = f[0].data  # This is the image array
    header = f[0].header

    # strong the image WCS into an object
    w = WCS(header)

    # Get the RA and Dec of the center of the image
    [raImage, decImage] = w.all_pix2world(data.shape[0] / 2, data.shape[1] / 2, 1)

    # get ra and dec of right, left, top and bottom of img for box size
    [raREdge, decREdge] = w.all_pix2world(data.shape[0] - int(crop), data.shape[1] / 2, 1)
    [raLEdge, decLEdge] = w.all_pix2world(int(crop), data.shape[1] / 2, 1)
    [raTop, decTop] = w.all_pix2world(data.shape[0] / 2, data.shape[1] - int(crop), 1)
    [raBot, decBot] = w.all_pix2world(data.shape[0] / 2, int(crop), 1)

    height = (decTop - decBot) * 60
    width = (raLEdge - raREdge) * 60

    # detect if GB field
    case = header['OBJTYPE']
    bulge = bulge_checker(case)

    # get weight map detection threshold
    try:
        chip = header['CHIP']
    except KeyError:
        chip = 4
    det_cut = get_weight_thresh(chip)

    return data, header, w, raImage, decImage, bulge, det_cut, chip


# %%
# gaia query for catalog completion check
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
        v = Vizier(columns=['RA_ICRS', 'DE_ICRS', 'RPmag'],
                   column_filters={"Dup": "<1", "Nd": ">6"},
                   row_limit=-1)
        G = v.query_region(SkyCoord(coords, unit=(u.deg, u.deg)), width=str(width) + 'm'
                           , catalog=catNum, cache=False, frame=chosen_frame)

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

        # Path = '/mnt/photometry/AT2025wgq/field9614-2025-09-10/J_rerun/stack/crsgaia.reg'
        # newtext = open(Path, 'w+')
        # newtext.write('fk5')
        # for i, j in zip(good_G_stars[idx_gaia][G_RA], good_G_stars[idx_gaia][G_DEC]):
        #     newtext.write('\npoint(%f,%f) # point=circle 5' % (i, j))
        #
        # Path = '/mnt/photometry/AT2025wgq/field9614-2025-09-10/J_rerun/stack/vhs.reg'
        # newtext = open(Path, 'w+')
        # newtext.write('fk5')
        # for i, j in zip(good_Q_stars[Q_RA], good_Q_stars[Q_DEC]):
        #     newtext.write('\npoint(%f,%f) # point=circle 5' % (i, j))

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


# Use astroquery to get catalog search
def query(raImage, decImage, band, w, data, crop, acc_comp_lvl=0.4,
          survey=None, given_catalog_path=None, mag_lower_lim=None, mag_upper_lim=None, bulge=False, no_check=False):
    # query box width
    width = PHOTOMETRY_QUERY_WIDTH

    # mag lower cutoff
    if mag_lower_lim:
        mag_low_cutoff = mag_lower_lim
    else:
        mag_low_cutoff = PHOTOMETRY_MAG_LOWER_LIMIT

    # mag upper cutoff
    if mag_upper_lim:
        mag_high_cutoff = mag_upper_lim
    else:
        mag_high_cutoff = PHOTOMETRY_MAG_UPPER_LIMIT

    if given_catalog_path:
        Q = ascii.read(given_catalog_path)
        Q = Q[Q[f'{band}MAG_{magtype}'] > mag_low_cutoff]
        chosen_survey = 'PRIME'
    else:

        if bulge:           # if bulge field, change to galactic coords for query
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

        # new automatic survey picking
        if not survey:

            # current catalogs
            catalog_dict = PHOTOMETRY_QUERY_CATALOGS

            catalogs = []
            if band == 'J' or band == 'H':
                for k, v in catalog_dict.items():
                    if v[0] == 'J':
                        catalogs.append((k, v[1]))
            elif band == 'Z':
                for k, v in catalog_dict.items():
                    if v[0] == 'Z':
                        catalogs.append((k, v[1]))
            elif band == 'Y':
                for k, v in catalog_dict.items():
                    if v[0] == 'Y':
                        catalogs.append((k, v[1]))
            else:
                print('Only J, H, Y, and Z band are supported!')

            # coords = SkyCoord(ra=[raImage], dec=[decImage], unit=(u.deg, u.deg))

            checkwidth = 28     # smaller radius of initial query, to better avoid cases of being on coverage edge
            # checkwidth covers only cropped part of chip, reducing chance of catalog only being in cropped away area

            # current columns
            v = Vizier(columns=['RAJ2000', 'DEJ2000', 'RAICRS', 'DEICRS', 'RA_ICRS', 'DE_ICRS', '%sap3' % band,
                                'e_%sap3' % band, '%smag' % band, '%smag3' % band, '%smag1' % band, 'e_%smag1' % band,
                                'e_%smag' % band, 'e_%smag3' % band, '%smag' % band.lower(),'e_%smag' % band.lower(),
                                '%sPSF' % band.lower(), 'e_%sPSF' % band.lower(),
                                '%spmag' % band.lower(), 'e_%spmag' % band.lower(),
                                'Mclass'])
            try:
                result = v.query_region(coords, width=str(checkwidth) + 'm', catalog=[f[1] for f in catalogs])
                test = result[0]
            except IndexError:
                print('Sadly, no current surveys available in current area in %s band' % band)
                raise Exception('No surveys available.')

            keys = result.format_table_list()

            print('%s, %s, Box Width: %s arcmin... '
                  '\n%s band initial query resulting in: \n%s' % (
                  frame_long_str, frame_lat_str, checkwidth, band, keys))

            keycheck = result.keys()

            last_idx = catalogs[-1][1]
            for chosen_survey, catNum in catalogs:
                for k in keycheck:
                    if catNum in k:
                        print('%s catalog found!' % k)
                        vhs_table = result[''.join(k)]
                        cols = vhs_table.colnames
                        vhs_band_col = vhs_table[cols[2]]

                        if not np.all(vhs_band_col.mask):
                            print('Survey has coverage in %s band!' % band)
                            print('Survey = %s' % k)

                            if chosen_survey in ['DES_Y', 'DES_Z', 'Skymapper']:
                                print('\nMag upper limit active due to survey choice: %s' % mag_high_cutoff)
                                mag_lims = f"{mag_low_cutoff:f}..{mag_high_cutoff:f}"
                            else:
                                mag_lims = f">{mag_low_cutoff:f}"

                            try:
                                if len(cols) == 4 or len(cols) < 4:
                                    used_cols = ['%s' % cols[0], '%s' % cols[1], '%s' % cols[2], '%s' % cols[3]]
                                elif len(cols) > 4:
                                    used_cols = ['%s' % cols[0], '%s' % cols[1], '%s' % cols[2], '%s' % cols[3], '%s' % cols[4]]

                                if catNum == last_idx:
                                    print('\nVizier catalogs exhausted, switching to local 2MASS query...')
                                    print('Local 2MASS Query around %s, %s, boxwidth %.2f arcmin, mag lim of %s'
                                          % (frame_long_str, frame_lat_str, width, mag_lims))
                                    Q = local_query_box(ra_center=raImage,
                                                        dec_center=decImage,
                                                        width=str(width) + 'm',
                                                        columns=["ra","dec",f"{band.lower()}mag",f"e_{band.lower()}mag"],
                                                        column_filters={
                                                                    f"{band.lower()}mag": f">{mag_low_cutoff:f}",
                                                                }
                                                        )
                                else:
                                    print('\nQuerying Vizier %s around %s, %s, boxwidth %.2f arcmin, mag lim of %s'
                                          % (catNum, frame_long_str, frame_lat_str, width, mag_lims))
                                    v = Vizier(columns=used_cols,
                                               column_filters={"%s" % cols[2]: mag_lims,
                                                               "%sFlag" % band.lower(): "<4",
                                                               "%sflags1" % band.lower(): "<16",
                                                               "%sperrbits" % band: '<128',
                                                               }, row_limit=-1)
                                    Q = v.query_region(SkyCoord(coords, unit=(u.deg, u.deg)),
                                                       width=str(width) + 'm',
                                                       catalog=k, cache=False, frame=chosen_frame)
                                if Q and len(Q[0]) > 0:
                                    print('Queried source total = ', len(Q[0]))
                                    if no_check:
                                        break
                                    elif not no_check and catNum != last_idx:
                                        gaia_comp = gaia_crsmtch_check(coords, width, chosen_frame, w, data, crop, Q)
                                        if gaia_comp >= acc_comp_lvl:
                                            print(f' Completion acceptable (>{acc_comp_lvl})! Moving on w/ catalog!\n')
                                            break
                                        else:
                                            print(f' Gaia completion w/ {catNum} < {acc_comp_lvl}, defaulting to next catalog..\n')
                                else:
                                    # in case survey provides no sources for some reason, try next available
                                    print(f"No sources found in {catNum}, trying fallback if available...")
                            except Exception as e:
                                print('Error in Vizier query.')
                                print(f"Error details: {e}")
                                continue
                        else:
                            print('Query unsuccessful, defaulting to next fallback catalog...')
                    else:
                        continue
                else:
                    continue
                break

        # Old functionality
        else:
            chosen_survey = survey
            if survey == '2MASS':
                catNum = 'II/246'  # changing to 2mass
                print('\nQuerying Vizier %s around %s, %s, boxwidth %.2f arcmin'
                      % (catNum, frame_long_str, frame_lat_str, width))
                try:
                    # You can set the bands for the individual columns (magnitude range, number of detections) inside the Vizier query
                    v = Vizier(columns=['RAJ2000', 'DEJ2000', '%smag' % band, 'e_%smag' % band],
                               column_filters={"%smag" % band: f">{mag_low_cutoff:f}", "Nd": ">6"},
                               row_limit=-1)
                    Q = v.query_region(SkyCoord(coords, unit=(u.deg, u.deg)), width=str(width) + 'm'
                                       , catalog=catNum, cache=False, frame=chosen_frame)
                    # query vizier around (ra, dec) with a radius of boxsize
                    # print(Q[0])
                    print('Queried source total = ', len(Q[0]))
                except:
                    print('I cannnot reach the Vizier database. Is the internet working?')
            elif survey == 'VHS':
                catNum = 'II/367'  # changing to vista
                print('\nQuerying Vizier %s around %s, %s, boxwidth %.2f arcmin'
                      % (catNum, frame_long_str, frame_lat_str, width))
                try:
                    # You can set the bands for the individual columns (magnitude range, number of detections) inside the Vizier query
                    v = Vizier(columns=['RAJ2000', 'DEJ2000', '%sap3' % band, 'e_%sap3' % band],
                               column_filters={"%sap3" % band: f">{mag_low_cutoff:f}"}, row_limit=-1)
                    Q = v.query_region(SkyCoord(coords, unit=(u.deg, u.deg)), width=str(width) + 'm'
                                       , catalog=catNum, cache=False, frame=chosen_frame)
                    # query vizier around (ra, dec) with a radius of boxsize
                    # print(Q[0])
                    print('Queried source total = ', len(Q[0]))
                except:
                    print(
                        'Error in Vizier query. Perhaps your image is not in the southern hemisphere sky?  H band is also not well covered!')
            elif survey == 'VIKING':
                catNum = 'II/343/viking2'
                print('\nQuerying Vizier %s around %s, %s, boxwidth %.2f arcmin'
                      % (catNum, frame_long_str, frame_lat_str, width))
                try:
                    v = Vizier(columns=['RAJ2000', 'DEJ2000', '%sap3' % band, 'e_%sap3' % band],
                               column_filters={"%sap3" % band: f">{mag_low_cutoff:f}"}, row_limit=-1)
                    Q = v.query_region(SkyCoord(coords, unit=(u.deg, u.deg)), width=str(width) + 'm'
                                       , catalog=catNum, cache=False, frame=chosen_frame)
                    print('Queried source total = ', len(Q[0]))
                except:
                    print(
                        'Error in Vizier query. Perhaps your image is not in the southern hemisphere sky?  H band is also not well covered!'
                        ' If you are in S.H., VIKING is only in a relatively smaller strip!')
            elif survey == 'Skymapper':
                catNum = 'II/379/smssdr4'
                print('\nQuerying Vizier %s around %s, %s, boxwidth %.2f arcmin, mag lim of %s - %s'
                      % (catNum, frame_long_str, frame_lat_str, width, mag_low_cutoff, mag_high_cutoff))
                try:
                    v = Vizier(columns=['RAICRS', 'DEICRS', '%sPSF' % band.lower(), 'e_%sPSF' % band.lower()],
                               column_filters={"%sPSF" % band.lower(): f"{mag_low_cutoff:f}..{mag_high_cutoff:f}",
                                               "%sFlag" % band.lower(): "<4"
                                               }, row_limit=-1)
                    Q = v.query_region(SkyCoord(coords, unit=(u.deg, u.deg)), width=str(width) + 'm'
                                       , catalog=catNum, cache=False, frame=chosen_frame)
                    print('Queried source total = ', len(Q[0]))
                except:
                    print(
                        'Error in Vizier query. Perhaps your image is not in the southern hemisphere sky?'
                        '\n perhaps check skymapper coverage maps?')
            elif survey == 'SDSS':
                catNum = 'V/154/sdss16'
                print('\nQuerying Vizier %s around %s, %s, boxwidth %.2f arcmin, mag lim of %s - %s'
                      % (catNum, frame_long_str, frame_lat_str, width, mag_low_cutoff, mag_high_cutoff))
                try:
                    v = Vizier(columns=['RA_ICRS', 'DE_ICRS', '%spmag' % band.lower(), 'e_%spmag' % band.lower()],
                               column_filters={"%spmag" % band.lower(): f"{mag_low_cutoff:f}..{mag_high_cutoff:f}",
                                               "%sFlag" % band.lower(): "<4", "clean": "=1"}, row_limit=-1)
                    Q = v.query_region(SkyCoord(coords, unit=(u.deg, u.deg)), width=str(width) + 'm'
                                       , catalog=catNum, cache=False, frame=chosen_frame)
                    print('Queried source total = ', len(Q[0]))
                except:
                    print(
                        'Error in Vizier query. Perhaps your image is not in the southern hemisphere sky?'
                        '\n perhaps check SDSS coverage maps?')
            elif survey == 'UKIDSS':
                catNum = 'II/319/las9'
                print('\nQuerying Vizier %s around %s, %s, boxwidth %.2f arcmin'
                      % (catNum, frame_long_str, frame_lat_str, width))
                try:
                    v = Vizier(columns=['RAJ2000', 'DEJ2000', '%smag' % band, 'e_%smag' % band],
                               column_filters={"%smag" % band: f">{mag_low_cutoff:f}"}, row_limit=-1)
                    Q = v.query_region(SkyCoord(ra=raImage, dec=decImage, unit=(u.deg, u.deg)), width=str(width) + 'm',
                                       catalog=catNum, cache=False, frame=chosen_frame)
                    print('Queried source total = ', len(Q[0]))
                except:
                    print(
                        'Error in Vizier query. Perhaps your image is not in the southern hemisphere sky?'
                        '\n perhaps check UKIDSS coverage maps?')
            elif survey == 'PanSTARRS':
                catNum = 'II/389/ps1_dr2'
                print('\nQuerying Vizier %s around %s, %s, boxwidth %.2f arcmin'
                      % (catNum, frame_long_str, frame_lat_str, width))
                try:
                    v = Vizier(columns=['RAJ2000', 'DEJ2000', '%smag' % band.lower(), 'e_%smag' % band.lower()],
                               column_filters={"%smag" % band: f">{mag_low_cutoff:f}"}, row_limit=-1)
                    Q = v.query_region(SkyCoord(ra=raImage, dec=decImage, unit=(u.deg, u.deg)), width=str(width) + 'm',
                                       catalog=catNum, cache=False, frame=chosen_frame)
                    print('Queried source total = ', len(Q[0]))
                except:
                    print(
                        'Error in Vizier query. Perhaps your image is not in the southern hemisphere sky?'
                        '\n perhaps check UKIDSS coverage maps?')
            elif survey == 'DES':
                if band == 'Y':
                    catNum = 'II/371/des_dr2'
                    print('\nQuerying Vizier %s around %s, %s, boxwidth %.2f arcmin, mag lim of %s - %s'
                          % (catNum, frame_long_str, frame_lat_str, width, mag_low_cutoff, mag_high_cutoff))
                    try:
                        v = Vizier(columns=['RA_ICRS', 'DE_ICRS', '%smag' % band, 'e_%smag' % band],
                                   column_filters={"%smag" % band: f"{mag_low_cutoff:f}..{mag_high_cutoff:f}",
                                                   "%sFlag" % band.lower(): "<4"}, row_limit=-1)
                        Q = v.query_region(SkyCoord(coords, unit=(u.deg, u.deg)), width=str(width) + 'm',
                                           catalog=catNum, cache=False, frame=chosen_frame)
                        print('Queried source total = ', len(Q[0]))
                    except:
                        print(
                            'Error in Vizier query. Perhaps your image is not in the southern hemisphere sky?'
                            '\n perhaps check DES coverage maps?')
                elif band == 'Z':
                    catNum = 'II/371/des_dr2'
                    print('\nQuerying Vizier %s around %s, %s, boxwidth %.2f arcmin, mag lim of %s - %s'
                          % (catNum, frame_long_str, frame_lat_str, width, mag_low_cutoff, mag_high_cutoff))
                    try:
                        v = Vizier(columns=['RA_ICRS', 'DE_ICRS', '%smag' % band.lower(), 'e_%smag' % band.lower()],
                                   column_filters={"%smag" % band: f"{mag_low_cutoff:f}..{mag_high_cutoff:f}",
                                                   "%sFlag" % band.lower(): "<4"}, row_limit=-1)
                        Q = v.query_region(SkyCoord(coords, unit=(u.deg, u.deg)), width=str(width) + 'm',
                                           catalog=catNum, cache=False, frame=chosen_frame)
                        print('Queried source total = ', len(Q[0]))
                    except:
                        print(
                            'Error in Vizier query. Perhaps your image is not in the southern hemisphere sky?'
                            '\n perhaps check DES coverage maps?')
                else:
                    raise Exception('DES only supports Z and Y band!')
            else:
                print('No supported survey found, currently use either 2MASS, VHS, VIKING, Skymapper, SDSS, UKIDSS, or DES')
                raise Exception('No surveys found.')

    return Q, chosen_survey, mag_low_cutoff


# %%
# run sextractor on swarped img to find sources


def sex1(imageName, det_cut, grb_flag=False):
    print('Running sextractor on img to initially find sources...')
    aper_str = ''

    if grb_flag:
        configFile = gen_config_file_name('sex2.config')
        paramName = gen_config_file_name('tempsource.param')
        catalogName = imageName + '.cat'
    else:
        configFile = gen_config_file_name('sex2.config')
        paramName = gen_config_file_name('photomAUTO.param')
        catalogName = imageName + '.photom.cat'

        # single fixed aperture (for now)
        aper_size = 2.5 / 0.498     # 2.5" in pix
        aper_str = f'-PHOT_APERTURES {aper_size}'
        with fits.open(imageName, mode='update') as hdul:
            hdr = hdul[0].header
            hdr.set('APERS', f'{aper_size}',
                    'fixed aperture size (pix)')
            hdul.close()

    weightName = 'weight'+imageName[5:]
    if os.path.isfile(weightName):
        # imghdr = fits.getheader(imageName)
        # if 'BUNIT' in imghdr:
        #     scale_fac = imghdr['CONV_FAC']
        # else:
        scale_fac = 1
        weightdata = fits.getdata(weightName)
        weightdata = weightdata / scale_fac**2
        weight_med = np.nanmedian(weightdata)
        weight_std = np.nanstd(weightdata)
        detect_cutoff = weight_med - (weight_std * det_cut)
        with fits.open(weightName, mode='update') as hdu:
            whdr = hdu[0].header
            whdr.set('MEDIAN', weight_med, 'Median of weight image', after='EQUINOX')
            whdr.set('DET_CUT', detect_cutoff, 'Pix value cutoff for source detection', after='MEDIAN')
            hdu.close()
        try:
            print('Including weight map!')
            command = (f'sex %s -c %s -CATALOG_NAME %s -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_THRESH %s -WEIGHT_IMAGE %s -PARAMETERS_NAME %s {aper_str}'
                       %
                       (f'"{imageName}[0]"', configFile, catalogName, detect_cutoff, f'"{weightName}[0]"', paramName))
            # print('Executing command: %s' % command)
            subprocess.run(command.split(), check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except subprocess.CalledProcessError as err:
            print('Could not run sextractor with exit error %s' % err)
            raise Exception('Sextractor failed to run, is the stacked image quality adequate?')
    else:
        try:
            command = ('sex %s -c %s -CATALOG_NAME %s -PARAMETERS_NAME %s' %
                       (f'"{imageName}[0]"', configFile, catalogName, paramName))
            # print('Executing command: %s' % command)
            subprocess.run(command.split(), check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except subprocess.CalledProcessError as err:
            print('Could not run sextractor with exit error %s' % err)
            raise Exception('Sextractor failed to run, is the stacked image quality adequate?')
    return catalogName


# %%
# run psfex on sextractor LDAC from previous step


def psfex(catalogName, band, data, crop):
    print('Running PSFex on sextrctr catalogue to generate psf for stars in the img...')
    psfConfigFile = gen_config_file_name('default.psfex')
    psfImageName = 'PSF' + catalogName[5:-4]

    # bright mag prune
    rough_mag_correction = CHIP_ZPS[band]
    hdul = fits.open(catalogName)

    init_data = Table(hdul['LDAC_OBJECTS'].data)
    init_data['MAG_AUTO'] = init_data['MAG_AUTO'] + rough_mag_correction
    pruned_data = init_data[(init_data['MAG_AUTO'] > 14)]
    # pruned_data = init_data

    # additional pruning
    max_x = data.shape[0]
    max_y = data.shape[1]

    if isinstance(pruned_data['FLUX_RADIUS'][0], np.ndarray):
        flux_radius = pruned_data['FLUX_RADIUS'][:, 1]
    else:
        flux_radius = pruned_data['FLUX_RADIUS']

    # pruned_data = pruned_data[
    #     (pruned_data['XWIN_IMAGE'] < (max_x - crop)) & (pruned_data['XWIN_IMAGE'] > crop) &
    #     (pruned_data['YWIN_IMAGE'] < (max_y) - crop) & (pruned_data['YWIN_IMAGE'] > crop) & (flux_radius > 1 / 0.498)]

    # print(len(pruned_data))

    # Recreate LDAC format
    new_data = fits.BinTableHDU(pruned_data, name='LDAC_OBJECTS')
    new_hdul = fits.HDUList([
        hdul[0],  # Primary
        hdul['LDAC_IMHEAD'],  # SExtractor image header
        new_data  # Filtered catalog
    ])

    # prunedcatalogName = catalogName.replace('.cat', '.prune.cat')
    new_hdul.writeto(catalogName, overwrite=True)

    try:
        command = ('psfex %s -c %s -CHECKIMAGE_TYPE SNAPSHOTS -CHECKIMAGE_NAME %s'
                   % (catalogName, psfConfigFile,'PSF.fits'))
        # print('Executing command: %s' % command)
        subprocess.run(command.split(), check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError as err:
        print('Could not run psfex with exit error %s' % err)
    os.rename('PSF_' + catalogName[:-4] + '.fits', psfImageName)


# %%
# feed generated psf model back into sextractor w/ diff param (or could use that param from the start but its slower)


def sex2(imageName, det_cut, catalogName):
    # dynamic aperture photometry adjustment
    # init_cat = get_table_from_ldac(catalogName)
    # if isinstance(init_cat['FLUX_RADIUS'][0], np.ndarray):
    #     acc_sources = init_cat[(init_cat['FLUX_RADIUS'][:, 0] > 1 / 0.498)]
    #     hwhm = np.nanmedian(acc_sources['FLUX_RADIUS'][:, 0])
    #     fwhm = 2 * hwhm
    # else:
    #     acc_sources = init_cat[(init_cat['FLUX_RADIUS'][:, 0] > 1 / 0.498)]
    #     fwhm = np.nanmedian(acc_sources['FLUX_RADIUS'])
    #
    # aper_arr = [round(1 * fwhm,2), round(1.25 * fwhm,2), round(1.5 * fwhm,2), round(1.75 * fwhm,2), round(2 * fwhm,2)]

    # single fixed aperture (for now)
    aper_size = 2.5 / 0.498     # 2.5" in pix
    aper_str = f'-PHOT_APERTURES {aper_size}'
    # print(f'Photometric apertures (pix) adjusted dynamically to seeing: {aper_arr}')
    with fits.open(imageName, mode='update') as hdul:
        hdr = hdul[0].header
        hdr.set('APERS', f'{aper_size}',
                'fixed aperture size (pix)')
        hdul.close()

    print('Feeding psf model back into sextractor for fitting and flux calculation...')
    psfName = imageName + '.psf'
    psfcatalogName = imageName.replace('.fits', '.photom.cat')
    configFile = gen_config_file_name('sex2.config')
    psfparamName = gen_config_file_name('photomPSF.param')
    weightName = 'weight' + imageName[5:]
    if os.path.isfile(weightName):
        # imghdr = fits.getheader(imageName)
        # if 'BUNIT' in imghdr:
        #     scale_fac = imghdr['CONV_FAC']
        # else:
        scale_fac = 1
        weightdata = fits.getdata(weightName)
        weightdata = weightdata / scale_fac**2
        weight_med = np.nanmedian(weightdata)
        weight_std = np.nanstd(weightdata)
        detect_cutoff = weight_med - (weight_std * det_cut)
        try:
            # We are supplying SExtactor with the PSF model with the PSF_NAME option
            command = (f'sex "{imageName}[0]" -c {configFile} -CATALOG_NAME {psfcatalogName} -WEIGHT_TYPE MAP_WEIGHT '
                       f'-WEIGHT_THRESH {detect_cutoff} -WEIGHT_IMAGE "{weightName}[0]" -PSF_NAME {psfName} '
                       f'-PARAMETERS_NAME {psfparamName} {aper_str}')
            # print("Executing command: %s" % command)
            subprocess.run(command.split(), check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except subprocess.CalledProcessError as err:
            print('Could not run sextractor with exit error %s' % err)
            raise Exception('Is there a problem with the sextractor configs or PSF model? Recommend temporarily removing '
                  '"stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL" from this subprocess command to investigate.')
    else:
        try:
            # We are supplying SExtactor with the PSF model with the PSF_NAME option
            command = (f'sex "{imageName}[0]" -c {configFile} -CATALOG_NAME {psfcatalogName} -PSF_NAME {psfName} '
                       f'-PARAMETERS_NAME {psfparamName} {aper_str}')
            # print("Executing command: %s" % command)
            subprocess.run(command.split(), check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except subprocess.CalledProcessError as err:
            print('Could not run sextractor with exit error %s' % err)
            raise Exception('Is there a problem with the sextractor configs or PSF model? Recommend temporarily removing '
                  '"stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL" from this subprocess command to investigate.')
    return psfcatalogName


# %%
# read in tables, crossmatch
def tables(Q, data, w, psfcatalogName, crop, given_catalog_path=None):
    print('Cross-matching sextracted and catalog sources...')
    crop = int(crop)
    max_x = data.shape[0]
    max_y = data.shape[1]
    if given_catalog_path:
        given_cat = Q

        RA = 'ALPHA_J2000'
        DEC = 'DELTA_J2000'

        good_cat_stars = given_cat[(given_cat['FLAGS'] == 0) &
                               (given_cat['XWIN_IMAGE'] < (max_x - crop)) & (given_cat['XWIN_IMAGE'] > crop) &
                               (given_cat['YWIN_IMAGE'] < (max_y) - crop) & (given_cat['YWIN_IMAGE'] > crop)]

        if 'POINTSOURCE' in given_cat.colnames:
            good_cat_stars = good_cat_stars[(good_cat_stars['FLAGS_MODEL'] == 0)]

        try:
            psfsourceTable = get_table_from_ldac(psfcatalogName)
        except FileNotFoundError:
            raise FileNotFoundError(f'{psfcatalogName} not found! Require this file for -grb_only / int cal functionality! Rerun photometry w/ '
                     f'the "-keep" flag.')
        if isinstance(psfsourceTable['FLUX_RADIUS'][0], np.ndarray):
            r50 = psfsourceTable['FLUX_RADIUS'][:, 0]
            r90 = psfsourceTable['FLUX_RADIUS'][:, 1]

            idx_r = psfsourceTable.colnames.index('FLUX_RADIUS')
            psfsourceTable['FLUX_RADIUS'] = r50
            psfsourceTable.add_column(r90, name='FLUX_RADIUS_90', index=idx_r + 1)
            flux_radius = psfsourceTable['FLUX_RADIUS_90']
        else:
            flux_radius = psfsourceTable['FLUX_RADIUS']

        psfsourceTable = psfsourceTable[~np.isnan(flux_radius)]

        PSFSources = psfsourceTable[
            (psfsourceTable['XWIN_IMAGE'] < (max_x - crop)) & (psfsourceTable['XWIN_IMAGE'] > crop)
            & (psfsourceTable['YWIN_IMAGE'] < (max_y) - crop) & (psfsourceTable['YWIN_IMAGE'] > crop) &
            (flux_radius >= 1 / 0.498)]

        cleanPSFSources = psfsourceTable[
            (psfsourceTable['FLAGS'] == 0) & (psfsourceTable['XWIN_IMAGE']
            < (max_x - crop)) & (psfsourceTable['XWIN_IMAGE'] > crop) & (psfsourceTable['YWIN_IMAGE']
            < (max_y) - crop) & (psfsourceTable['YWIN_IMAGE'] > crop) & (flux_radius >= 1 / 0.498)]

        if 'POINTSOURCE' in cleanPSFSources.colnames:
            cleanPSFSources = cleanPSFSources[(cleanPSFSources['FLAGS_MODEL'] == 0)]

        print('Catalogue cropped, source total = ', len(good_cat_stars))

        psfsourceCatCoords = SkyCoord(ra=cleanPSFSources[RA], dec=cleanPSFSources[DEC], frame='icrs', unit='degree')
        massCatCoords = SkyCoord(ra=good_cat_stars[RA], dec=good_cat_stars[DEC], frame='icrs', unit='degree')

        photoDistThresh = 0.6

        idx_psfimage, idx_psfmass, d2d, d3d = massCatCoords.search_around_sky(psfsourceCatCoords,
                                                                              photoDistThresh * u.arcsec)
        print('Found %d good cross-matches' % len(idx_psfmass))

    else:
        colnames = Q[0].colnames
        RA = colnames[0]
        DEC = colnames[1]

        mass_imCoords = w.all_world2pix(Q[0][RA], Q[0][DEC], 1)
        crop_cat_stars = Q[0][np.where(
            (mass_imCoords[0] > crop) & (mass_imCoords[0] < (max_x - crop)) & (mass_imCoords[1] > crop) & (
                        mass_imCoords[1] < (max_y - crop)))]
        print('Approximate catalogue source total in image bounds = ',len(crop_cat_stars))

        good_cat_stars = Q[0].copy()

        try:
            psfsourceTable = get_table_from_ldac(psfcatalogName)
        except FileNotFoundError:
            raise FileNotFoundError(f'{psfcatalogName} not found! Require this file for -grb_only functionality! Rerun photometry w/ '
                     f'the "-keep" flag.')

        if isinstance(psfsourceTable['FLUX_RADIUS'][0], np.ndarray):
            r50 = psfsourceTable['FLUX_RADIUS'][:, 0]
            r90 = psfsourceTable['FLUX_RADIUS'][:, 1]

            idx_r = psfsourceTable.colnames.index('FLUX_RADIUS')
            psfsourceTable['FLUX_RADIUS'] = r50
            psfsourceTable.add_column(r90, name='FLUX_RADIUS_90', index=idx_r + 1)
            flux_radius = psfsourceTable['FLUX_RADIUS_90']
        else:
            flux_radius = psfsourceTable['FLUX_RADIUS']

        if isinstance(psfsourceTable['MAG_APER'][0], np.ndarray):
            mag_aper_col = psfsourceTable['MAG_APER'][:, 0]
            mag_aper_err_col = psfsourceTable['MAGERR_APER'][:, 0]
            flux_aper_col = psfsourceTable['FLUX_APER'][:, 0]
            flux_aper_err_col = psfsourceTable['FLUXERR_APER'][:, 0]

            psfsourceTable['MAG_APER'] = mag_aper_col
            psfsourceTable['MAGERR_APER'] = mag_aper_err_col
            psfsourceTable['FLUX_APER'] = flux_aper_col
            psfsourceTable['FLUXERR_APER'] = flux_aper_err_col
            # idx_apmag = psfsourceTable.colnames.index('MAG_APER')
            # idx_e_apmag = psfsourceTable.colnames.index('MAGERR_APER')
            # idx_apflx = psfsourceTable.colnames.index('FLUX_APER')
            # idx_e_apflx = psfsourceTable.colnames.index('FLUXERR_APER')
            #
            # psfsourceTable.add_column(mag_aper_col, name='MAG_APER', index=idx_apmag)
            # psfsourceTable.add_column(mag_aper_err_col, name='MAGERR_APER', index=idx_e_apmag)
            # psfsourceTable.add_column(flux_aper_col, name='FLUX_APER', index=idx_apflx)
            # psfsourceTable.add_column(flux_aper_err_col, name='FLUXERR_APER', index=idx_e_apflx)

        PSFSources = psfsourceTable[
            (flux_radius >= 1 / 0.498)]

        # PSFSources = psfsourceTable[
        #     (psfsourceTable['XWIN_IMAGE'] < (max_x - crop)) & (psfsourceTable['XWIN_IMAGE'] > crop)
        #     & (psfsourceTable['YWIN_IMAGE'] < (max_y) - crop) & (psfsourceTable['YWIN_IMAGE'] > crop) &
        #     (flux_radius >= 1 / 0.498)]

        cleanPSFSources = psfsourceTable[
            (psfsourceTable['XWIN_IMAGE']< (max_x - crop)) & (psfsourceTable['XWIN_IMAGE'] > crop) &
             (psfsourceTable['YWIN_IMAGE'] < (max_y) - crop) & (psfsourceTable['YWIN_IMAGE'] > crop) &
             (flux_radius >= 1 / 0.498) &
             (psfsourceTable['FLAGS'] == 0)]

        # if magtype == 'PSF':
        #     cleanPSFSources = cleanPSFSources[(cleanPSFSources['FLAGS_MODEL'] == 0)]

        # psfsourceCatCoords = SkyCoord(ra=cleanPSFSources['ALPHA_J2000'], dec=cleanPSFSources['DELTA_J2000'], frame='icrs', unit='degree')

        psfsourceCatCoords = utils.pixel_to_skycoord(cleanPSFSources['X_IMAGE'], cleanPSFSources['Y_IMAGE'], w, origin=1)

        massCatCoords = SkyCoord(ra=good_cat_stars[RA], dec=good_cat_stars[DEC], frame='icrs', unit='degree')

        # Now crossmatch sources
        # Set the cross-match distance threshold to 0.6 arcsec, or just about one pixel
        photoDistThresh = 0.6
        # photoDistThresh = 1.0
        idx_psfimage, idx_psfmass, d2d, d3d = massCatCoords.search_around_sky(psfsourceCatCoords,
                                                                              photoDistThresh * u.arcsec)
        # idx_psfimage are indexes into psfsourceCatCoords for the matched sources, while idx_psfmass are indexes into massCatCoords for the matched sources

        if 'Mclass' in colnames:
            # pruning crossmatches to only include stars if applicable
            print('Found %d good cross-matches before pruning' % len(idx_psfmass))
            star_mask = np.isin(good_cat_stars[colnames[4]], [-1, -2])
            star_matches = star_mask[idx_psfmass]
            idx_psfimage = idx_psfimage[star_matches]
            idx_psfmass = idx_psfmass[star_matches]
            print('Found %d good cross-matches after galaxy pruning (ONLY for zp calc crossmatch)' % len(idx_psfmass))
        else:
            print('Found %d good cross-matches' % len(idx_psfmass))

    return good_cat_stars, cleanPSFSources, PSFSources, idx_psfmass, idx_psfimage, massCatCoords, crop_cat_stars


def queryexport(good_cat_stars, imageName, survey):
    table = good_cat_stars
    chip = imageName[-6]
    num = imageName[-16:-8]
    table.write('%s_C%s_query_%s.ecsv' % (survey, chip, num), overwrite=True)
    print('%s_C%s_query_%s.ecsv written!' % (survey, chip, num))


# parallel running - filename switching for independent psf fit
def end_name_gen(ext='ecsv'):
    if magtype:
        end_name = f'_{magtype}.{ext}'
    else:
        end_name = f'.{ext}'
    return end_name


# parallel running - 2nd hdu generation for independent psf fit

def new_hdu_gen_or_set(hdul):
    hdr = hdul[0].header

    if magtype == 'PSF' and parallel:
        if len(hdul) == 1:
            new_hdu = ImageHDU(data=np.zeros((1, 1)), header=hdr.copy())
            hdul.append(new_hdu)

        target_hdr = hdul[1].header
    else:
        target_hdr = hdr

    return target_hdr

#%% NED galaxy pruning


def gal_match(raImage, decImage):
    coords = SkyCoord(ra=[raImage], dec=[decImage], unit=(u.deg, u.deg))

    print('\nQuerying NED around RA %.4f, Dec %.4f with a radius of 34 arcmin' % (raImage, decImage))
    try:
        gal_query = Ned.query_region(coords, radius=34 * u.arcmin)
        print('NED source total = ', len(gal_query))
    except:
        print('Error in NED query.')


# %%
# derive zero pt / put in swarped header

def zp_write(PSFSources, imageName, zp, e_zp, calmag, calmagerr, weights_noclip, clipped, band, magtype=defaults['magtype']):
    """
    writes ZP info & conversion factors to header and ecsv
    """

    zero_mean = zp
    zero_std = e_zp

    cal_mag_name = f'{band}MAG_{MAGTYPES[magtype]}'
    cal_mag_err_name = f'e_{band}MAG_{MAGTYPES[magtype]}'
    print(f'{MAGTYPES[magtype]} zp = %.4f, zp err = %.6f' % (zero_mean, zero_std))

    tempmagcol = MaskedColumn(calmag, name=cal_mag_name, unit=u.ABmag)
    flux = ab_to_microjy(tempmagcol)
    fluxcol = MaskedColumn(flux, name=f'{MAGTYPES[magtype]}_FLUX_DENSITY', unit=u.microjansky)
    PSFSources.add_column(fluxcol)

    # image / col data conversion to uJy
    if magtype == MAGTYPES['PSF']:
        flx_type = 'POINTSOURCE'
    else:
        flx_type = MAGTYPES[magtype]

    conv_factor_all = fluxcol / PSFSources[f'FLUX_{flx_type}']  # u = uJy / adu
    conv_factor = np.nanmedian(conv_factor_all)

    fluxerr_ujy = PSFSources[f'FLUX_{flx_type}'] * conv_factor
    fluxerrujycol = MaskedColumn(fluxerr_ujy, name=f'E_{MAGTYPES[magtype]}_FLUX_DENSITY', unit=u.microjansky)
    PSFSources.add_column(fluxerrujycol)

    with open_fits_robust(imageName) as hdul:
        hdr = hdul[0].header
        try:
            hdr.set(f'ZP_{MAGTYPES[magtype]}', zero_mean, f'{MAGTYPES[magtype]} Zero Point Offset', after='NINT')
        except KeyError:
            hdr.set(f'ZP_{MAGTYPES[magtype]}', zero_mean, f'{MAGTYPES[magtype]} Zero Point Offset')
        hdr.set(f'eZP_{MAGTYPES[magtype]}', zero_std, f'{MAGTYPES[magtype]} Zero Point Offset Error', after=f'ZP_{MAGTYPES[magtype]}')

        # hdr.set('BUNIT', 'uJy', 'Physical units of the array values IF multiplied by conv_fac', after='EXTEND')
        hdr.set(f'FAC_{MAGTYPES[magtype]}', conv_factor.item(),
                'uJy / ADU PSF Conversion Factor, multiply img by this to get in uJy', before='EQUINOX')
        print(f'Conversion of ADU to uJy calculated for {MAGTYPES[magtype]}, med conversion factor: %.4f' % conv_factor.item())

    magcol = MaskedColumn(calmag, name=cal_mag_name, unit=u.ABmag)
    magerrcol = MaskedColumn(calmagerr, name=cal_mag_err_name, unit=u.ABmag)

    PSFSources.add_column(magcol)
    PSFSources.add_column(magerrcol)

    if magtype == MAGTYPES['APER']:
        aper_size_str = ', 2.5"'
    else:
        aper_size_str = ''

    PSFSources[f'{band}MAG_{MAGTYPES[magtype]}'].description = f'{band} band {MAGTYPES[magtype]} model magnitude{aper_size_str}'
    PSFSources[f'e_{band}MAG_{MAGTYPES[magtype]}'].description = f'Error in {band} band {MAGTYPES[magtype]} model magnitude{aper_size_str}'

    return PSFSources, weights_noclip, clipped


def single_zeropt(good_cat_stars, crop_cat_stars, cleanPSFSources, PSFSources, idx_mass, idx_image, imageName, band, survey, sigma, data,
           crop):
    """
    Calculates and writes zero pt statistics to header for specific magtype
    """

    # getting survey catalog column names
    colnames = good_cat_stars.colnames
    if len(colnames) > 6:
        magcolname = f'{band}MAG_{MAGTYPES[magtype]}'
        magerrcolname = f'e_{band}MAG_{MAGTYPES[magtype]}'
    else:
        magcolname = colnames[2]
        magerrcolname = colnames[3]

    print('Converting Vega surveys to AB to ensure correctness!')
    good_cat_stars[magcolname] = ab_convert(good_cat_stars[magcolname], band=band, survey=survey)
    ab_cat_stars = good_cat_stars

    # calculating zp for specific mag columns
    if magtype == MAGTYPES['PSF']:
        prime_magcolname = 'POINTSOURCE'
    else:
        prime_magcolname = MAGTYPES[magtype]

    prime_mags = [(cleanPSFSources[f'MAG_{prime_magcolname}'][idx_image])]
    primeerr = [(cleanPSFSources[f'MAGERR_{prime_magcolname}'][idx_image])]

    cat_mags = good_cat_stars[magcolname][idx_mass]
    caterr = good_cat_stars[magerrcolname][idx_mass]

    for prime_mag_col, prime_err_col in zip(prime_mags, primeerr):
        if prime_mag_col.ndim == 1:  # standard column zp calc

            mag_col_name = prime_mag_col.name
            mag_err_col_name = prime_err_col.name
            print(f'Running ZP calc for {mag_col_name}, {mag_err_col_name}')

            comberr = np.sqrt(caterr ** 2 + prime_err_col ** 2)

            weights_noclip = 1 / (comberr ** 2)

            offsets = ma.array(cat_mags - prime_mag_col)
            offsets = offsets.data

            # 3 sigma clip
            clipped = sigma_clip(offsets, sigma=sigma)
            offsets = offsets[~clipped.mask]
            weights = np.array(weights_noclip[~clipped.mask])
            print('\nZero point source offsets clipped by %s sigma, total clipped offset # = %s' % (
                sigma, len(offsets)))

            # Compute statistics
            zero_mean = sum(offsets * weights) / sum(weights)
            zero_std = np.sqrt(1 / sum(weights))

            # catalog for all detected sources
            calmag = zero_mean + PSFSources[mag_col_name]
            calmagerr = np.sqrt(PSFSources[mag_err_col_name] ** 2 + zero_std ** 2)

            # writing all zp info to hdr
            PSFSources, weights_noclip, clipped = zp_write(
                                                            PSFSources, imageName, zero_mean, zero_std, calmag, calmagerr,
                                                            weights_noclip, clipped, band, magtype=MAGTYPES[magtype]
                                                            )
        else:
            raise Exception('Multidimensional mag column detected!  Currently, multi MAG_APER is not supported by '
                            'single_zeropoint.')

        # writing other pertinent info to hdr
        print('\nWriting image stats to image header & catalog...')
        with open_fits_robust(imageName) as hdul:
            hdr = hdul[0].header
            hdr.set('N_CRSMCH', len(idx_image), 'Number of Crossmatches', after='NINT')
            hdr.set('N_SRCS', len(PSFSources), 'Total PRIME Sources Number', after='N_CRSMCH')
            hdr.set('SURVEY', survey, 'Chosen Survey for Crossmatch', after='N_SRCS')

        if 'VIGNET' in PSFSources.colnames:
            PSFSources.remove_column('VIGNET')
        PSFSources['FLUX_RADIUS'] = PSFSources['FLUX_RADIUS'] * 0.498
        PSFSources['FLUX_RADIUS'].unit = u.arcsec
        if 'FLUX_RADIUS_90' in PSFSources.colnames:
            PSFSources['FLUX_RADIUS_90'] = PSFSources['FLUX_RADIUS_90'] * 0.498
            PSFSources['FLUX_RADIUS_90'].unit = u.arcsec
            PSFSources['FLUX_RADIUS_90'].description = '90% flux radius'
        print('Total PRIME source # = ', len(PSFSources))

        # efficiency calculation
        gen_efficiency(imageName, PSFSources, band, specific_magtype=MAGTYPES[magtype])

        # col descriptions
        PSFSources['ALPHA_J2000'].description = 'J2000 RA coordinate'
        PSFSources['DELTA_J2000'].description = 'J2000 Dec coordinate'
        PSFSources['FLUX_RADIUS'].description = 'HWHM, 50% flux radius'
        PSFSources['SNR_WIN'].description = 'SNR in a Gaussian window'
        PSFSources['ELONGATION'].description = 'semi-major axis / semi-minor axis'

        # catalog for clean sources
        crop = int(crop)
        max_x = data.shape[0]
        max_y = data.shape[1]

        cleanPSFSources = PSFSources[
            (PSFSources['XWIN_IMAGE'] < (max_x - crop)) & (PSFSources['XWIN_IMAGE'] > crop) &
            (PSFSources['YWIN_IMAGE'] < (max_y) - crop) & (PSFSources['YWIN_IMAGE'] > crop) &
            (PSFSources['FLAGS'] == 0)]

        # writing out full source catalog
        PSFSources.write('%s.%s.%s.ecsv' % (imageName, survey, magtype), overwrite=True)
        print('%s.%s.%s.ecsv written, CSV w/ corrected mags' % (imageName, survey, MAGTYPES[magtype]))

        if abs(len(idx_mass) / len(crop_cat_stars)) < 0.03:
            raise Exception(f'*WARNING* Significant disparity in cropped survey vs. cross-matched source num (<0.03), '
                            f'crsmtch/survey ratio: {round(abs(len(idx_mass) / len(ab_cat_stars)), 3)}'
                            f'\nPhotometry likely NOT reliable, Recommend checking field, '
                            f'likely an astrometric & or stacking issue!')

        return cleanPSFSources, PSFSources, weights_noclip, clipped, ab_cat_stars


def zeropt(good_cat_stars, crop_cat_stars, cleanPSFSources, PSFSources, idx_psfmass, idx_psfimage, imageName, band, survey, sigma, data,
           crop):
    colnames = good_cat_stars.colnames

    if len(colnames) > 6:
        magcolname = f'{band}MAG_{magtype}'
        magerrcolname = f'e_{band}MAG_{magtype}'
    else:
        magcolname = colnames[2]
        magerrcolname = colnames[3]

    print('Converting Vega surveys to AB to ensure correctness!')
    good_cat_stars[magcolname] = ab_convert(good_cat_stars[magcolname], band=band, survey=survey)
    ab_cat_stars = good_cat_stars

    # getting adjusted apertures
    try:
        apers_str = fits.getheader(imageName)['APERS']
        aper_arr = apers_str.split(',')
        aper_arr = [float(ap) for ap in aper_arr]
    except KeyError:
        apers_str = aper_arr = 0

    # calculating zp for all appropriate mag columns

    prime_psf_mags = [(cleanPSFSources['MAG_AUTO'][idx_psfimage])]
    primeerr = [(cleanPSFSources['MAGERR_AUTO'][idx_psfimage])]

    if magtype == MAGTYPES['PSF']:
        prime_psf_mags.append(cleanPSFSources['MAG_POINTSOURCE'][idx_psfimage])
        primeerr.append(cleanPSFSources['MAGERR_POINTSOURCE'][idx_psfimage])
    # if 'MAG_APER' in PSFSources.colnames:
    #     prime_psf_mags.append(cleanPSFSources['MAG_APER'][idx_psfimage])
    #     primeerr.append(cleanPSFSources['MAGERR_APER'][idx_psfimage])

    cat_mags = good_cat_stars[magcolname][idx_psfmass]
    caterr = good_cat_stars[magerrcolname][idx_psfmass]

    psfweights_noclip = []
    psf_clipped = []
    autoweights_noclip = []
    for prime_mag_col, prime_err_col in zip(prime_psf_mags, primeerr):
        if prime_mag_col.ndim == 1:    # standard column zp calc

            mag_col_name = prime_mag_col.name
            mag_err_col_name = prime_err_col.name

            comberr = np.sqrt(caterr ** 2 + prime_err_col ** 2)

            weights_noclip = 1 / (comberr ** 2)

            offsets = ma.array(cat_mags - prime_mag_col)
            offsets = offsets.data

            # 3 sigma clip
            clipped = sigma_clip(offsets, sigma=sigma)
            offsets = offsets[~clipped.mask]
            weights = np.array(weights_noclip[~clipped.mask])
            print('\nZero point source offsets clipped by %s sigma, total clipped offset # = %s' % (sigma, len(offsets)))

            # Compute statistics
            # zero_psfmean = np.average(psfoffsets, weights=psfweights)
            zero_mean = sum(offsets * weights) / sum(weights)
            # zero_psfvar = np.average((psfoffsets - zero_psfmean) ** 2, weights=psfweights)
            # zero_psfstd = np.sqrt(zero_psfvar)
            zero_std = np.sqrt(1 / sum(weights))

            zero_apermeans = []
            zero_aperstds = []
            aperweights_noclip = []
            aper_clipped_all = []

            # catalog for all detected sources
            calmag = zero_mean + PSFSources[mag_col_name]
            calmagerr = np.sqrt(PSFSources[mag_err_col_name] ** 2 + zero_std ** 2)

            # print('Converting mags from Vega to AB for all sources! (if not already in AB)')
            # calmag = ab_convert(calmag, band=band, survey=survey)

            if 'POINTSOURCE' in mag_col_name:

                PSFSources, psfweights_noclip, psf_clipped = zp_write(PSFSources, imageName, zero_mean, zero_std, calmag, calmagerr,
                                                                      weights_noclip, clipped, band, magtype=MAGTYPES['PSF'])

            else:

                PSFSources, autoweights_noclip, auto_clipped = zp_write(PSFSources, imageName, zero_mean, zero_std, calmag, calmagerr,
                                                                      weights_noclip, clipped, band, magtype=MAGTYPES['AUTO'])

            # real unit flux conversion

            # if not 'FLUX_DENSITY' in PSFSources.colnames:
            #
            #     psfflux = magcol.to(u.microjansky)
            #     psffluxcol = MaskedColumn(psfflux, name='PSF_FLUX_DENSITY', unit=u.microjansky)
            #     PSFSources.add_column(psffluxcol)
            #
            #     # image / col data conversion to uJy
            #     conv_factor_all = psffluxcol / PSFSources['FLUX_POINTSOURCE']  # u = uJy / adu
            #     conv_factor = np.nanmedian(conv_factor_all)
            #
            #     fluxerr_ujy = PSFSources['FLUXERR_POINTSOURCE'] * conv_factor
            #     fluxerrujycol = MaskedColumn(fluxerr_ujy, name='E_PSF_FLUX_DENSITY', unit=u.microjansky)
            #     PSFSources.add_column(fluxerrujycol)
            #
            #     with fits.open(imageName, mode='update') as imagehdu:
            #         imagehdr = imagehdu[0].header
            #         # if you replace the pix values, put the if statement back in
            #         # imagehdu[0].data = imagehdu[0].data * conv_factor  # adu * (uJy / adu) = uJy
            #         imagehdr.set('BUNIT', 'uJy', 'Physical units of the array values IF multiplied by conv_fac', after='EXTEND')
            #         imagehdr.set('CONV_FAC', conv_factor, 'uJy / ADU Conversion Factor, multiply img by this to get in uJy', after='BUNIT')
            #         print('Conversion of ADU to uJy calculated, med conversion factor: %.4f' % conv_factor)
            #
            #         # print('BUNIT found already in header, skipping conversion.')
            #         imagehdu.close()

        else:   # vector column MAG_APER zp calc

            valid_apermask = prime_mag_col < 90

            aper_comberr = np.sqrt(caterr[:, np.newaxis] ** 2 + prime_err_col ** 2)

            aperweights_noclip = 1 / (aper_comberr ** 2)

            aperoffsets = ma.array(cat_mags[:, np.newaxis] - prime_mag_col)
            aperoffsets = aperoffsets.data

            # sigma clip & zp calc
            zero_apermeans = []
            zero_aperstds = []
            aper_clipped_all = []

            print('\nFixed aperture ZPs below:')
            for i in range(aperoffsets.shape[1]):
                aperoffsets_i = aperoffsets[:, i]
                aperweights_i = aperweights_noclip[:, i]

                aper_clipped = sigma_clip(aperoffsets_i, sigma=sigma)
                aper_clipped_all.append(aper_clipped)
                mask = ~aper_clipped.mask
                aperoffsets_i = aperoffsets_i[mask]
                aperweights_i = np.array(aperweights_i[mask])

                zero_apermean = sum(aperoffsets_i * aperweights_i) / sum(aperweights_i)
                zero_aperstd = np.sqrt(1 / sum(aperweights_i))

                zero_apermeans.append(zero_apermean)
                zero_aperstds.append(zero_aperstd)

                print(f'    {round(aper_arr[i] * 0.498,2)}" aper zp = %.4f, zp err = %.6f' % (zero_apermean, zero_aperstd))

            zero_apermeans = np.array(zero_apermeans)
            zero_aperstds = np.array(zero_aperstds)

            cal_apermags = zero_apermeans + PSFSources['MAG_APER']
            apermagerrs = np.sqrt(PSFSources['MAGERR_APER'] ** 2 + zero_aperstds ** 2)

            apermags = []
            for i in range(cal_apermags.shape[1]):
                # apmag = ab_convert(cal_apermags[:, i], band=band, survey=survey)
                apmag = cal_apermags[:, i]
                apermags.append(apmag)

            apermags = np.column_stack(apermags)

            apermagcol = MaskedColumn(apermags, name='%sMAG_APER' % band, unit=u.ABmag)
            apermagerrcol = MaskedColumn(apermagerrs, name='e_%sMAG_APER' % band, unit=u.ABmag)

            PSFSources.add_column(apermagcol)
            PSFSources.add_column(apermagerrcol)

            PSFSources['%sMAG_APER' % band].description = (f'{band} Band aperture magnitudes: '
                                                           f'{round(aper_arr[0], 2)}, {round(aper_arr[1], 2)}, {round(aper_arr[2], 2)},'
                                                           f' {round(aper_arr[3], 2)}, {round(aper_arr[4], 2)} diameters')
            PSFSources['e_%sMAG_APER' % band].description = \
                (f'Error in {band} band aperture magnitudes (pix)')

    # zero_psfmean, zero_psfmed, zero_psfstd = sigma_clipped_stats(psfoffsets)
    # print('PSF Mean ZP: %.2f\nPSF Median ZP: %.2f\nPSF STD ZP: %.2f'%(zero_psfmean, zero_psfmed, zero_psfstd))

    # writing zp to header
    print('\nWriting ZP info to image header...')
    with fits.open(imageName, mode='update') as hdul:
        hdr = hdul[0].header
        hdr.set('N_CRSMCH', len(idx_psfimage), 'Number of Crossmatches', after='NINT')
        hdr.set('N_SRCS', len(PSFSources), 'Total PRIME Sources Number', after='N_CRSMCH')
        hdr.set('SURVEY', survey, 'Chosen Survey for Crossmatch', after='N_SRCS')

        if len(zero_apermeans) > 0:
            for idx, (zp, err) in enumerate(zip(zero_apermeans, zero_aperstds)):
                hdr.set(f'e_ZPap{idx}', err, f'{(idx + 1) * 2}" Aperture Zero Point Offset Error', after='e_ZP_PSF')
                hdr.set(f'ZPap{idx}', zp, f'{(idx + 1) * 2}" Aperture Zero Point Offset', after='e_ZP_PSF')
        hdul.close()

    if 'VIGNET' in PSFSources.colnames:
        PSFSources.remove_column('VIGNET')
    PSFSources['FLUX_RADIUS'] = PSFSources['FLUX_RADIUS'] * 0.498
    PSFSources['FLUX_RADIUS'].unit = u.arcsec
    if 'FLUX_RADIUS_90' in PSFSources.colnames:
        PSFSources['FLUX_RADIUS_90'] = PSFSources['FLUX_RADIUS_90'] * 0.498
        PSFSources['FLUX_RADIUS_90'].unit = u.arcsec
        PSFSources['FLUX_RADIUS_90'].description = '90% flux radius'
    print('Total PRIME source # = ', len(PSFSources))

    # efficiency calculation
    gen_efficiency(imageName, PSFSources, band)

    # col descriptions
    PSFSources['ALPHA_J2000'].description = 'J2000 RA coordinate'
    PSFSources['DELTA_J2000'].description = 'J2000 Dec coordinate'
    PSFSources['FLUX_RADIUS'].description = 'HWHM, 50% flux radius'
    PSFSources['SNR_WIN'].description = 'SNR in a Gaussian window'
    PSFSources['ELONGATION'].description = 'semi-major axis / semi-minor axis'

    # PSFSources['%sMAG_PSF' % band] = -2.5 * np.log10(PSFSources['FLUX_DENSITY']) + 23.9
    # PSFSources['e_%sMAG_PSF' % band] = 1.086 * (PSFSources['E_FLUX_DENSITY'] / PSFSources['FLUX_DENSITY'])

    # catalog for clean sources

    crop = int(crop)
    max_x = data.shape[0]
    max_y = data.shape[1]

    cleanPSFSources = PSFSources[
        (PSFSources['XWIN_IMAGE'] < (max_x - crop)) & (PSFSources['XWIN_IMAGE'] > crop) &
        (PSFSources['YWIN_IMAGE'] < (max_y) - crop) & (PSFSources['YWIN_IMAGE'] > crop) &
        (PSFSources['FLAGS'] == 0)]

    if magtype == 'PSF' and parallel:
        PSFSources.write('%s.%s.%s.ecsv' % (imageName, survey, magtype), overwrite=True)
        print('%s.%s.%s.ecsv written, CSV w/ corrected mags' % (imageName, survey, magtype))
    else:
        if magtype == 'PSF':
            PSFSources.write('%s.%s.%s.ecsv' % (imageName, survey, magtype), overwrite=True)
            print('%s.%s.%s.ecsv written, CSV w/ corrected mags' % (imageName, survey, magtype))
        PSFSources.write('%s.%s.ecsv' % (imageName, survey), overwrite=True)
        print('%s.%s.ecsv written, CSV w/ corrected mags' % (imageName, survey))

    # catalog conversion to AB

    # print('Converting Vega surveys to AB to ensure correctness!\n')
    # good_cat_stars[magcolname] = ab_convert(good_cat_stars[magcolname], band=band, survey=survey)
    # ab_cat_stars = good_cat_stars

    if abs(len(idx_psfmass) / len(crop_cat_stars)) < 0.03:
        raise Exception(f'*WARNING* Significant disparity in cropped survey vs. cross-matched source num (<0.03), '
                        f'survey/crsmtch ratio: {round(abs(len(idx_psfmass) / len(ab_cat_stars)), 3)}'
                        f'\nPhotometry likely NOT reliable, Recommend checking field, '
                        f'likely an astrometric & or stacking issue!')

    return (cleanPSFSources, PSFSources, psfweights_noclip, psf_clipped, ab_cat_stars, aperweights_noclip, aper_clipped_all,
            autoweights_noclip, auto_clipped)


#%%


# New Source Search
def newsourcesearch(ra, dec, imageName, survey, band, thresh, massCatCoords, good_cat_stars, directory, chip, grbname=defaults['grb_name'],
                    mag_low_lim=PHOTOMETRY_MAG_LOWER_LIMIT,  **kwargs):

    source_ra = ra
    source_dec = dec
    ab_cat_stars = good_cat_stars
    header = fits.getheader(imageName)
    w = WCS(header)

    if grbname != defaults['grb_name']:
        grb_name = grbname
        newsrcname = grbname
    else:
        grb_name = defaults['grb_name']
        newsrcname = 'New_PRIME'

    if len(imageName) <= 16:
        num = 'img'
    else:
        num = imageName[-16:-8]

    # large region postage stamp cutout fctn (png & fits)
    def grb_cutout(imageName, GRBcoords, photoDistThresh, regprimename=None, regsurvname=None):
        imgdata = fits.getdata(imageName)
        img = fits.open(imageName)
        head = img[0].header
        w = WCS(head)

        savename = '%s_%s_Cutout_%s_%s' % (grb_name, band, survey, num)
        threshname = '%s_query_thresh.reg' % grb_name

        size = 4 * photoDistThresh * u.arcsec
        try:
            cutout = Cutout2D(imgdata, GRBcoords, size, wcs=w, copy=True)
            region = CircleSkyRegion(center=GRBcoords[0], radius=Angle(thresh, unit='arcsec'))
            pix_region = region.to_pixel(cutout.wcs)
        except astropy.nddata.utils.NoOverlapError:
            print(' Area of GRB threshold not found within image, cannot generate cutout!')
            return savename, threshname

        mean, median, sigma_cut = sigma_clipped_stats(cutout.data)
        plt.figure(10, figsize=(8, 8))
        plt.imshow(cutout.data, vmin=median - 3 * sigma_cut, vmax=median + 3 * sigma_cut, origin='lower', cmap='viridis')
        pix_region.plot(color='cyan', ls='--', label='Input GRB threshold')

        if regprimename:
            primeregs = open('%s_srcs.reg' % newsrcname, 'r')
            plt_primeregs = []
            primeallregs = [reg for reg in primeregs if reg != 'fk5\n']
            for reg in primeallregs:
                nums = re.findall(r'[-+]?\d*\.?\d+', reg)
                srcra = float(nums[0])
                srcdec = float(nums[1])
                srcrad = float(nums[2])
                srccoords = SkyCoord(ra=[srcra], dec=[srcdec], frame='icrs', unit='degree')
                srcreg = CircleSkyRegion(center=srccoords[0], radius=Angle(srcrad, unit='arcsec'))
                plt_primeregs.append(srcreg)

            for reg in plt_primeregs:
                pix_reg = reg.to_pixel(cutout.wcs)
                pix_reg.plot(color='red', ls='-', label='PRIME Source')

        if regsurvname:
            survregs = open('%s_%s_srcs.reg' % (grb_name, survey), 'r')
            plt_survregs = []
            survallregs = [reg for reg in survregs if reg != 'fk5\n']
            for reg in survallregs:
                nums = re.findall(r'[-+]?\d*\.?\d+', reg)
                srcra = float(nums[0])
                srcdec = float(nums[1])
                srcrad = float(nums[2])
                srccoords = SkyCoord(ra=[srcra], dec=[srcdec], frame='icrs', unit='degree')
                srcreg = CircleSkyRegion(center=srccoords[0], radius=Angle(srcrad, unit='arcsec'))
                plt_survregs.append(srcreg)

            for reg in plt_survregs:
                pix_reg = reg.to_pixel(cutout.wcs)
                pix_reg.plot(color='magenta', ls='-', label='%s Source' % survey)

        handles, labels = plt.gca().get_legend_handles_labels()
        seen = set()
        filtered_handles = []
        filtered_labels = []
        for h, l in zip(handles, labels):
            if l not in seen:
                filtered_handles.append(h)
                filtered_labels.append(l)
                seen.add(l)
        plt.legend(filtered_handles, filtered_labels, loc='best')
        plt.savefig(savename + '.png', dpi=300)
        plt.close()

        fits.writeto(savename + '.fits', cutout.data, cutout.wcs.to_header(), overwrite=True)
        return savename, threshname

    # new source reg gen
    def source_reg_gen(src_ra=0, src_dec=0, rad=2, src_survey=None, append=False):
        if not src_survey:
            name = '%s_srcs.reg' % newsrcname
            color = 'red'
        else:
            name = '%s_%s_srcs.reg' % (grb_name, survey)
            color = 'yellow'
        if not append:
            newtext = open(name, 'w+')
            newtext.write('fk5')
            if src_ra != 0 and src_dec != 0:
                newtext.write(f'\ncircle({src_ra}, {src_dec}, {rad}") # color={color}')
        else:
            newtext = open(name, 'a')
            newtext.write(f'\ncircle({src_ra}, {src_dec}, {rad}") # color={color}')

        return name

    # html gen for new sources
    def html_gen(data, directory, savename, threshname, band, survey, ra, dec, thresh, primename=None):
        df = data.to_pandas()
        tbl_html = df.to_html(index=False, classes="my-table")

        with open(savename + '.png', "rb") as img_file:
            encoded = base64.b64encode(img_file.read()).decode("utf-8")

        fits_items = [savename + '.fits', threshname, primename]
        fits_items = [os.path.join(directory, f) for f in fits_items if f is not None]

        base = fits_items[0]
        regions = " ".join(f"-regions {f}" for f in fits_items[1:])

        ds9_command = f"ds9 {base} {regions} &"

        final_html = f"""
        <div style="text-align: center; font-family: Arial, sans-serif;">

            <h2>GRB Information</h2>
            <p style="margin-top: 0; margin-bottom: 10px; font-size: 16px; color: #555;">
                RA = {ra}, Dec = {dec}, threshold = {thresh}"
            </p>

            <img src="data:image/png;base64,{encoded}" alt="GRB Cutout Region" width="600" style="margin-bottom: 10px;">

            <p style="margin-top: 0px; margin-bottom: 10px; font-size: 16px; color: #555;">
                To see source regions, open GRB stamp in DS9 through terminal: 
            </p>
            <p style="margin-top: 0px; margin-bottom: 20px; font-size: 14px; color: #030303;">
                {ds9_command}
            </p>

            <div style="display: inline-block; text-align: left;">
                {tbl_html}
            </div>

        </div>
        """

        with open('%s_Source_%s_Data_%s_%s.html' % (newsrcname, band, survey, num), "w") as f:
            f.write(final_html)

    # crossmatch for all detected sources
    print('Large grb radius inputted, using new source search to find'
          ' detected sources brighter than survey lim mag w/ no crossmatch')

    # sexigesimal conversion
    try:
        float(source_ra)
        source_ra = source_ra
        source_dec = source_dec
        deci_sky_coords = SkyCoord(ra=[source_ra], dec=[source_dec], frame='icrs', unit='degree')
    except ValueError:
        coords = source_ra + ' ' + source_dec
        print('Sexagesimal RA = %s & Dec = %s' % (source_ra, source_dec))

        deci_sky_coords = SkyCoord(coords, frame='icrs', unit=(u.hourangle, u.deg))
        deci_coords = deci_sky_coords.to_string()
        deci_coords = deci_coords.split(' ')

        source_ra = deci_coords[0]
        source_dec = deci_coords[1]

    if len(ab_cat_stars.colnames) > 6:
        RA = 'ALPHA_J2000'
        DEC = 'DELTA_J2000'
    else:
        RA = ab_cat_stars.colnames[0]
        DEC = ab_cat_stars.colnames[1]

    massCatCoords = SkyCoord(ra=ab_cat_stars[RA], dec=ab_cat_stars[DEC], frame='icrs', unit='degree')
    print('Catalog cropped #:', len(massCatCoords))

    # initial crossmatch
    mag_ecsvname = '%s.%s.ecsv' % (imageName, survey)
    mag_ecsvtable = ascii.read(mag_ecsvname)
    mag_ecsvSources = mag_ecsvtable[(mag_ecsvtable['FLAGS'] == 0) & (mag_ecsvtable['FLAGS_MODEL'] == 0)]
    # print(len(mag_ecsvSources))
    if RA == 'ALPHA_J2000':
        mag_ecsvsourceCatCoords = SkyCoord(ra=mag_ecsvSources['ALPHA_J2000'], dec=mag_ecsvSources['DELTA_J2000'],
                                           frame='icrs',
                                           unit='degree')
    else:
        mag_ecsvsourceCatCoords = utils.pixel_to_skycoord(mag_ecsvSources['X_IMAGE'], mag_ecsvSources['Y_IMAGE'], w, origin=1)

    photoDistThresh = 1.0   # set higher to combat offset wcs in certain sources, maybe change for denser fields?
    idx_psfimage_noclean, idx_psfmass_noclean, d2d, d3d = massCatCoords.search_around_sky(mag_ecsvsourceCatCoords,
                                                                          photoDistThresh * u.arcsec)

    mask = np.ones(len(mag_ecsvSources), dtype=bool)
    mask[idx_psfimage_noclean] = False
    PSFsources_nomatch = mag_ecsvSources[mask]  # removing previous crossmatched sources
    print('# of sources found after removing survey crossmatches: %i' % len(PSFsources_nomatch))

    sourcecoords = SkyCoord(ra=[source_ra], dec=[source_dec], frame='icrs', unit='degree')
    PSFsources_nomatchCatCoords = SkyCoord(ra=PSFsources_nomatch['ALPHA_J2000'], dec=PSFsources_nomatch['DELTA_J2000'],
                                       frame='icrs',
                                       unit='degree')
    idx_inputcoords, idx_PSFsources_nomatch, d2dd, d3dd = PSFsources_nomatchCatCoords.search_around_sky(sourcecoords,
                                                                                                      thresh * u.arcsec)
    PSFsources_nomatch = PSFsources_nomatch[idx_PSFsources_nomatch]     # implementing error radius
    print('# of sources found after implementing err radius: %i' % len(PSFsources_nomatch))

    # lim mag pruning
    for f in PHOTOMETRY_LIM_MAGS.keys():
        if survey == f:
            lim_mag = PHOTOMETRY_LIM_MAGS[f]
            break
    else:
        # limiting mag est.
        # lim_mag = ab_cat_stars['%sMAG_PSF' % band][(ab_cat_stars['SNR_WIN'] < 6) & (ab_cat_stars['SNR_WIN'] > 5)]
        lim_mag = 19.7
        print(f'Error in finding lim mag for survey catalogs, no matching catalog found? '
              f'Using generous PRIME lim: {lim_mag}')

    print(f' {survey} limiting mag: {lim_mag}')
    PSFsources_new = PSFsources_nomatch[(PSFsources_nomatch[f'{band}MAG_{MAGTYPES[magtype]}'] < lim_mag) &
                                        (PSFsources_nomatch[f'{band}MAG_{MAGTYPES[magtype]}'] > mag_low_lim)]
    print('# of sources found after removing sources dimmer than %.2f & brighter than %.2f in %s mag: %i'
          % (mag_low_lim, lim_mag, MAGTYPES[magtype], len(PSFsources_new)))

    # Table gen

    if len(PSFsources_new) > 0:
        PSFsources_new.write('%s_Sources.%s.%s.%s.ecsv' % (newsrcname, imageName, survey, num), overwrite=True)
        print('New source full catalog written!')

        all_magtypes = set(MAGTYPES.keys())
        prime_mag_cols = sorted([col for col in PSFsources_new.colnames if any(mag in col for mag in all_magtypes)
                                 and band in col and f'{band}MAG' in col and 'e_' not in col])

        mag_auto_ar = []
        mag_auto_err_ar = []
        mag_psf_ar = []
        mag_psf_err_ar = []
        mag_aper_ar = []
        mag_aper_err_ar = []
        ra_ar = []
        dec_ar = []
        rad_ar = []
        snr_ar = []
        source_reg_gen()

        for i in PSFsources_new:
            if f'{band}MAG_AUTO' in prime_mag_cols:
                grb_mag_auto = i[f'{band}MAG_AUTO']
                mag_auto_ar.append(grb_mag_auto)
                grb_mag_auto_err = i[f'e_{band}MAG_AUTO']
                mag_auto_err_ar.append(grb_mag_auto_err)
            if f'{band}MAG_PSF' in prime_mag_cols:
                grb_mag_psf = i[f'{band}MAG_PSF']
                mag_psf_ar.append(grb_mag_psf)
                grb_mag_psf_err = i[f'e_{band}MAG_PSF']
                mag_psf_err_ar.append(grb_mag_psf_err)
            if f'{band}MAG_APER' in prime_mag_cols:
                grb_mag_aper = i[f'{band}MAG_APER']
                mag_aper_ar.append(grb_mag_aper)
                grb_mag_aper_err = i[f'e_{band}MAG_APER']
                mag_aper_err_ar.append(grb_mag_aper_err)
            grb_ra = i['ALPHA_J2000']
            ra_ar.append(grb_ra)
            grb_dec = i['DELTA_J2000']
            dec_ar.append(grb_dec)
            grb_rad = i['FLUX_RADIUS']
            rad_ar.append(grb_rad)
            grb_snr = i['SNR_WIN']
            snr_ar.append(grb_snr)

            regprimename = source_reg_gen(src_ra=grb_ra, src_dec=grb_dec, rad=grb_rad, append=True)

        grbdata = Table()
        grbdata['RA'] = np.round(np.array(ra_ar), decimals=5) * u.deg
        grbdata['DEC'] = np.round(np.array(dec_ar), decimals=5) * u.deg
        if f'{band}MAG_AUTO' in prime_mag_cols:
            grbdata[f'{band}autoMag'] = np.round(np.array(mag_auto_ar), decimals=3) * u.ABmag
            grbdata[f'{band}autoMag_Err'] = np.round(np.array(mag_auto_err_ar), decimals=3) * u.ABmag
        if f'{band}MAG_PSF' in prime_mag_cols:
            grbdata[f'{band}psfMag'] = np.round(np.array(mag_psf_ar), decimals=3) * u.ABmag
            grbdata[f'{band}psfMag_Err'] = np.round(np.array(mag_psf_err_ar), decimals=3) * u.ABmag
        if f'{band}MAG_APER' in prime_mag_cols:
            grbdata[f'{band}aperMag'] = np.round(np.array(mag_aper_ar), decimals=3) * u.ABmag
            grbdata[f'{band}aperMag_Err'] = np.round(np.array(mag_aper_err_ar), decimals=3) * u.ABmag
        grbdata['Radius'] = np.round(np.array(rad_ar), decimals=2) * u.arcsec
        grbdata['SNR'] = np.round(np.array(snr_ar), decimals=2)

        print('New source catalog writen!')
        grbdata.write('%s_Source_%s_Data_%s_%s.ecsv' % (newsrcname, band, survey, num), overwrite=True)

        # generate stamp & html
        print('Generating location cutout and html files!')
        savename, threshname = grb_cutout(imageName=imageName, GRBcoords=deci_sky_coords, photoDistThresh=thresh,
                                          regprimename=regprimename)
        html_gen(grbdata, directory, savename, threshname, band, survey, ra=source_ra,
                 dec=source_dec, thresh=thresh, primename=regprimename)

    else:
        print('No sources remaining after pruning!  Cannot write new source catalog!')


# %% optional GRB-specific photom


def GRB(ra, dec, imageName, survey, band, thresh, massCatCoords, good_cat_stars, directory, chip, coordlist=None,
        grbname=defaults['grb_name'], **kwargs):
    mag_ecsvname = '%s.%s.ecsv' % (imageName, survey)
    mag_ecsvtable = ascii.read(mag_ecsvname)
    mag_ecsvcleanSources = mag_ecsvtable  # [(mag_ecsvtable['FLAGS'] == 0) & (mag_ecsvtable['FLAGS_MODEL'] == 0)]
    # mag_ecsvcleanSources['%sMAG_PSF' % band] = (
    #     ab_convert(mag_ecsvcleanSources['%sMAG_PSF' % band], band=band, survey=survey, revert=True))
    mag_ecsvsourceCatCoords = SkyCoord(ra=mag_ecsvcleanSources['ALPHA_J2000'], dec=mag_ecsvcleanSources['DELTA_J2000'],
                                       frame='icrs',
                                       unit='degree')

    if len(imageName) <= 16:
        num = 'img'
    else:
        num = imageName[-16:-8]

    # postage stamp cutout fctn (png & fits)
    def grb_cutout(imageName, GRBcoords, photoDistThresh, loc=None, append=False, regprimename=None, regsurvname=None):
        imgdata = fits.getdata(imageName)
        img = fits.open(imageName)
        head = img[0].header
        w = WCS(head)

        if loc:
            savename = '%s_%s_C%i_Cutout_%s_%s_loc_%d' % (grbname, band, chip, survey, num, loc)
            threshname = '%s_queries_thresh.reg' % grbname
        else:
            savename = '%s_%s_C%i_Cutout_%s_%s' % (grbname, band, chip, survey, num)
            threshname = '%s_query_thresh.reg' % grbname

        size = 4 * photoDistThresh * u.arcsec
        try:
            cutout = Cutout2D(imgdata, GRBcoords[0], size, wcs=w, copy=True)
            region = CircleSkyRegion(center=GRBcoords[0], radius=Angle(thresh, unit='arcsec'))
            pix_region = region.to_pixel(cutout.wcs)
        except astropy.nddata.utils.NoOverlapError:
            print(' Area of GRB threshold not found within image, cannot generate cutout!')
            return savename, threshname

        mean, median, sigma_cut = sigma_clipped_stats(cutout.data)
        plt.figure(10, figsize=(8, 8))
        plt.imshow(cutout.data, vmin=median - 3 * sigma_cut, vmax=median + 3 * sigma_cut, origin='lower',
                   cmap='viridis')
        pix_region.plot(color='cyan', ls='--', label='Input GRB threshold')

        if regprimename:
            primeregs = open('%s_PRIME_srcs.reg' % grbname, 'r')
            plt_primeregs = []
            primeallregs = [reg for reg in primeregs if reg != 'fk5\n']
            for reg in primeallregs:
                nums = re.findall(r'[-+]?\d*\.?\d+', reg)
                srcra = float(nums[0])
                srcdec = float(nums[1])
                srcrad = float(nums[2])
                srccoords = SkyCoord(ra=[srcra], dec=[srcdec], frame='icrs', unit='degree')
                srcreg = CircleSkyRegion(center=srccoords[0], radius=Angle(srcrad, unit='arcsec'))
                plt_primeregs.append(srcreg)

            for reg in plt_primeregs:
                pix_reg = reg.to_pixel(cutout.wcs)
                pix_reg.plot(color='red', ls='-', label='PRIME Source')

        if regsurvname:
            survregs = open('%s_%s_srcs.reg' % (grbname, survey), 'r')
            plt_survregs = []
            survallregs = [reg for reg in survregs if reg != 'fk5\n']
            for reg in survallregs:
                nums = re.findall(r'[-+]?\d*\.?\d+', reg)
                srcra = float(nums[0])
                srcdec = float(nums[1])
                srcrad = float(nums[2])
                srccoords = SkyCoord(ra=[srcra], dec=[srcdec], frame='icrs', unit='degree')
                srcreg = CircleSkyRegion(center=srccoords[0], radius=Angle(srcrad, unit='arcsec'))
                plt_survregs.append(srcreg)

            for reg in plt_survregs:
                pix_reg = reg.to_pixel(cutout.wcs)
                pix_reg.plot(color='magenta', ls='-', label='%s Source' % survey)

        handles, labels = plt.gca().get_legend_handles_labels()
        seen = set()
        filtered_handles = []
        filtered_labels = []
        for h, l in zip(handles, labels):
            if l not in seen:
                filtered_handles.append(h)
                filtered_labels.append(l)
                seen.add(l)
        plt.legend(filtered_handles, filtered_labels, loc='best')
        plt.savefig(savename + '.png', dpi=300)
        plt.close()

        fits.writeto(savename + '.fits', cutout.data, cutout.wcs.to_header(), overwrite=True)

        # ds9 regions
        if append:
            newtext = open(threshname, 'a')  # input threshold
            newtext.write(f'\ncircle({ra}, {dec}, {photoDistThresh}") # color=cyan width=2 text={{Query Thresh}}')
        else:
            newtext = open(threshname, 'w+')  # input threshold
            newtext.write('fk5')
            newtext.write(f'\ncircle({ra}, {dec}, {photoDistThresh}") # color=cyan width=2 text={{Query Thresh}}')

        print(' Exported cutout & ds9 region of GRB search area!')

        return savename, threshname

    # source ds9 region writing
    def source_reg_gen(src_ra=0, src_dec=0, rad=2, src_survey=None, append=False):
        if not src_survey:
            name = '%s_PRIME_srcs.reg' % grbname
            color = 'red'
        else:
            name = '%s_%s_srcs.reg' % (grbname, survey)
            color = 'yellow'
        if not append:
            newtext = open(name, 'w+')
            newtext.write('fk5')
            if src_ra != 0 and src_dec != 0:
                newtext.write(f'\ncircle({src_ra}, {src_dec}, {rad}") # color={color}')
        else:
            newtext = open(name, 'a')
            newtext.write(f'\ncircle({src_ra}, {src_dec}, {rad}") # color={color}')

        return name

    # mag diff calc betw. survey and prime for existing source crsmtches
    def mag_diff_calc(survey_cat, prime_cat, survey_idx, prime_idx, d2d, band):
        # prime & survey mag cols
        all_magtypes = set(MAGTYPES.keys())
        prime_mag_cols = sorted([col for col in prime_cat.colnames if any(mag in col for mag in all_magtypes)
                                 and band in col and f'{band}MAG' in col and 'e_' not in col])
        prime_err_cols = sorted([col for col in prime_cat.colnames if any(mag in col for mag in all_magtypes)
                                 and band in col and f'e_{band}MAG' in col])

        colnames = survey_cat.colnames

        # mag diff calc for all applicable cols
        mag_diff_ar = []
        flag_ar = []

        for prime_mags, prime_errs in zip(prime_mag_cols, prime_err_cols):
            col_magtype = prime_mags.split('_')[-1]

            if len(colnames) > 6:
                magcolname = f'{band}MAG_{col_magtype}'
                magerrcolname = f'e_{band}MAG_{col_magtype}'
            else:
                magcolname = colnames[2]
                magerrcolname = colnames[3]

            if len(prime_idx) > 1:
                if 0 in prime_idx:
                    if not np.isscalar(prime_cat[f'{band}MAG_{col_magtype}']):
                        mag_diff = float(survey_cat[magcolname][survey_idx][0]) - float(
                            prime_cat[f'{band}MAG_{col_magtype}'][prime_idx][0])

                        comb_err = np.sqrt(survey_cat[magerrcolname][survey_idx][0] ** 2 +
                                           prime_cat[f'e_{band}MAG_{col_magtype}'][prime_idx][0] ** 2)
                    else:
                        mag_diff = float(survey_cat[magcolname][survey_idx][0]) - float(prime_cat[f'{band}MAG_{col_magtype}'])

                        comb_err = np.sqrt(survey_cat[magerrcolname][survey_idx][0] ** 2 +
                                           prime_cat[f'e_{band}MAG_{col_magtype}'] ** 2)

                    sep = d2d[0]
                else:
                    mag_diff = float(survey_cat[magcolname][survey_idx]) - float(
                        prime_cat[f'{band}MAG_{col_magtype}'][prime_idx])
                    # errors in quad
                    comb_err = np.sqrt(survey_cat[magerrcolname][survey_idx] ** 2 +
                                       prime_cat[f'e_{band}MAG_{col_magtype}'][prime_idx] ** 2)
                    sep = d2d
            else:
                # survey mag - prime mag
                mag_diff = float(survey_cat[magcolname][survey_idx].item()) - float(prime_cat[f'{band}MAG_{col_magtype}'].item())

                # errors in quad
                comb_err = np.sqrt(survey_cat[magerrcolname][survey_idx] ** 2 +
                                   prime_cat[f'e_{band}MAG_{col_magtype}'] ** 2)

                sep = d2d

            err_3_sig = 3 * comb_err

            if abs(mag_diff) > err_3_sig:
                flg = 2
                # print(' At least 1 source in threshold has a significant mag difference to the catalog!')
            else:
                flg = 1

            if not np.isscalar(sep):
                sep = sep[0]
            elif len(sep) > 1:
                sep = sep[0]

            sep_dist = sep.to(u.arcsec)
            sep_dist = sep_dist / u.arcsec

            mag_diff_ar.append((mag_diff, col_magtype))
            # mag_diff_ar.append(mag_diff)
            flag_ar.append((np.int16(flg), col_magtype))

        return mag_diff_ar, float(sep_dist), flag_ar


    def mag_diff_ar_grabber(mag_diff_ar, flag_ar, output, col_magtype='AUTO'):
        # grabbing correct mag diff values
        if isinstance(mag_diff_ar[0], tuple):
            mag_diff_val = [t[0] for t in mag_diff_ar if t[1] == col_magtype]
            flag_val = [t[0] for t in flag_ar if t[1] == col_magtype]
            if output == 'mag':
                return mag_diff_val[0]
            elif output == 'flag':
                return flag_val[0]
        else:
            if output == 'mag':
                return 99
            elif output == 'flag':
                return np.int16(0)


    # generation of html file
    def html_gen(data, directory, savename, threshname, band, survey, ra, dec, thresh, survname=None, primename=None):
        # building table
        newtbl = data.copy()
        # newtbl.remove_columns([f'{band}apMag', f'{band}apMag_Err'])

        df = newtbl.to_pandas()

        descriptions = {}
        for col in newtbl.colnames:
            desc = newtbl[col].description
            if desc is None or str(desc).strip() == "":
                desc = f"Description of {col}"
            descriptions[col] = desc

        init_html = df.to_html(index=False)

        soup = BeautifulSoup(init_html, "html.parser")

        for th in soup.find_all("th"):
            col_name = th.text.strip()
            if col_name in descriptions:
                th['title'] = descriptions[col_name]

        tbl_html = str(soup)

        # adding img & command
        with open(savename + '.png', "rb") as img_file:
            encoded = base64.b64encode(img_file.read()).decode("utf-8")

        fits_items = [savename + '.fits', threshname, survname, primename]
        fits_items = [os.path.join(directory, f) for f in fits_items if f is not None]

        base = fits_items[0]
        regions = " ".join(f"-regions {f}" for f in fits_items[1:])

        ds9_command = f"ds9 {base} {regions} &"

        # final html gen
        final_html = f"""
        <div style="text-align: center; font-family: Arial, sans-serif;">

            <h2>GRB Information</h2>
            <p style="margin-top: 0; margin-bottom: 10px; font-size: 16px; color: #555;">
                RA = {ra}, Dec = {dec}, threshold = {thresh}"
            </p>

            <img src="data:image/png;base64,{encoded}" alt="GRB Cutout Region" width="600" style="margin-bottom: 10px;">

            <p style="margin-top: 0px; margin-bottom: 10px; font-size: 16px; color: #555;">
                To see source regions, open GRB stamp in DS9 through terminal: 
            </p>
            <p style="margin-top: 0px; margin-bottom: 20px; font-size: 14px; color: #030303;">
                {ds9_command}
            </p>

            <div style="display: inline-block; text-align: left;">
                {tbl_html}
            </div>

        </div>
        """

        with open('%s_Multisource_%s_C%i_Data_%s_%s.html' % (grbname, band, chip, survey, num), "w") as f:
            f.write(final_html)

    # sexigesimal conversion
    try:
        float(ra)
        ra = ra
        dec = dec
    except ValueError:
        coords = ra + ' ' + dec
        print('Sexagesimal RA = %s & Dec = %s' % (ra, dec))

        deci_coords = SkyCoord(coords, frame='icrs', unit=(u.hourangle, u.deg)).to_string()
        deci_coords = deci_coords.split(' ')

        ra = deci_coords[0]
        dec = deci_coords[1]

    photoDistThresh = thresh
    if coordlist:
        for i in range(len(coordlist)):
            print('Checking GRB location %i: %s' % (i, coordlist[i]))
        coordlist = tuple(eval(i) for i in coordlist)
        try:
            idx_GRBpsfdict = {}
            keys = np.arange(0, len(coordlist), 1)
            for i in range(len(coordlist)):
                GRBcoords = SkyCoord(ra=[coordlist[i][0]], dec=[coordlist[i][1]], frame='icrs', unit='degree')
                idx_GRB, idx_GRBcleanpsf, d2d, d3d = mag_ecsvsourceCatCoords.search_around_sky(GRBcoords,
                                                                                               photoDistThresh * u.arcsec)

                # survey source crsmtch
                idx_survey, idx_surveycleanpsf, d2d_surv, d3d_surv = massCatCoords.search_around_sky(GRBcoords,
                                                                                                     photoDistThresh * u.arcsec)
                if len(idx_surveycleanpsf) > 0:
                    print(' %i %s existing sources found within GRB threshold! Writing to DS9 reg files...'
                          % (len(idx_surveycleanpsf), survey))
                    source_reg_gen(src_survey=survey)
                    for idx in idx_surveycleanpsf:
                        src = massCatCoords[idx]
                        source_reg_gen(src.ra.deg, src.dec.deg, src_survey=survey, append=True)
                else:
                    print(' No existing %s sources found within GRB threshold!' % survey)

                if i == 0:
                    grb_cutout(imageName, GRBcoords, photoDistThresh, loc=i)
                else:
                    grb_cutout(imageName, GRBcoords, photoDistThresh, loc=i, append=True)

                idx_GRBpsfdict[keys[i]] = idx_GRBcleanpsf, coordlist[i]
        except NameError:
            print('No Sources found!')

    else:
        print('Checking GRB location %s, %s' % (ra, dec))
        GRBcoords = SkyCoord(ra=[ra], dec=[dec], frame='icrs', unit='degree')
        # prime source crsmtch
        idx_GRB, idx_GRBcleanpsf, d2d, d3d = mag_ecsvsourceCatCoords.search_around_sky(GRBcoords,
                                                                                       photoDistThresh * u.arcsec)
        # survey source crsmtch
        idx_survey, idx_surveycleanpsf, d2d_surv, d3d_surv = massCatCoords.search_around_sky(GRBcoords,
                                                                                             photoDistThresh * u.arcsec)
        if len(idx_surveycleanpsf) > 0:
            print(' %i %s existing sources found within GRB threshold! Writing to DS9 reg files...'
                  % (len(idx_surveycleanpsf), survey))
            source_reg_gen(src_survey=survey)
            for idx in idx_surveycleanpsf:
                src = massCatCoords[idx]
                regsurvname = source_reg_gen(src.ra.deg, src.dec.deg, src_survey=survey, append=True)
        else:
            print(' No existing %s sources found within GRB threshold!' % survey)
            regsurvname = None

        savename, threshname = grb_cutout(imageName, GRBcoords, photoDistThresh)

    # custom col descriptions
    desc = {
        "mag_aper": f'{band} band 2.5" aperture magnitude',
        "mag_aper_err": f'Error in {band} band 2.5" aperture magnitude',
        "mag_err_crsmtch": f'{survey} - PRIME source mag for survey crossmatched source, 99 if no crossmatch',
        "distance": '2D distance (arcsec) betw. PRIME source & input GRB coords',
        "separation": f'2D distance (arcsec) betw. PRIME & crossmatched {survey} source, -1 = no match',
        "crsmtch_flg": (f'Flag for {survey} crossmatch: '
                        f'0 = no match, 1 = match w/ mag diff within 3 sig, '
                        f'2 = match w/ mag diff outside 3 sig')
    }

    mag_col_num = len([col for col in mag_ecsvcleanSources.colnames if any(mag in col for mag in set(MAGTYPES.keys()))
                             and band in col and f'{band}MAG' in col and 'e_' not in col])

    def get_desc(table, colname):
        """Returns col description if col exists"""
        return table[colname].description if colname in table.colnames else None

    # grb table initiation
    grbdata = Table()

    # grb table descriptions, always present
    ra_desc = mag_ecsvcleanSources['ALPHA_J2000'].description
    dec_desc = mag_ecsvcleanSources['DELTA_J2000'].description
    rad_desc = mag_ecsvcleanSources['FLUX_RADIUS'].description
    snr_desc = mag_ecsvcleanSources['SNR_WIN'].description
    elon_desc = mag_ecsvcleanSources['ELONGATION'].description
    dist_desc = desc['distance']
    sep_desc = desc['separation']
    flag_desc = desc['crsmtch_flg']
    mag_diff_desc = desc['mag_err_crsmtch']

    # grb table descriptions, mags
    psf_mag_desc = get_desc(mag_ecsvcleanSources, '%sMAG_PSF' % band)
    psf_mag_err_desc = get_desc(mag_ecsvcleanSources, 'e_%sMAG_PSF' % band)
    aper_mag_desc = get_desc(mag_ecsvcleanSources, '%sMAG_APER' % band)
    aper_mag_err_desc = get_desc(mag_ecsvcleanSources, 'e_%sMAG_APER' % band)
    auto_mag_desc = get_desc(mag_ecsvcleanSources, '%sMAG_AUTO' % band)
    auto_mag_err_desc = get_desc(mag_ecsvcleanSources, 'e_%sMAG_AUTO' % band)

    if coordlist:
        source_reg_gen()
        for key in idx_GRBpsfdict:
            values = idx_GRBpsfdict[key]
            idx_GRBcleanpsf = values[0]
            ra = values[1][0]
            dec = values[1][1]
            print(' idx size = %d for location %d' % (len(idx_GRBcleanpsf), key))
            if len(idx_GRBcleanpsf) == 1:
                print(' GRB source at inputted coords %s and %s, rad = %s arcsec found!' % (ra, dec, photoDistThresh))
                if psf_mag_desc is not None:
                    grb_mag = mag_ecsvcleanSources[idx_GRBcleanpsf]['%sMAG_PSF' % band][0]
                    grb_magerr = mag_ecsvcleanSources[idx_GRBcleanpsf]['e_%sMAG_PSF' % band][0]
                if aper_mag_desc is not None:
                    grb_mag_aper = mag_ecsvcleanSources[idx_GRBcleanpsf]['%sMAG_APER' % band][0]
                    grb_magerr_aper = mag_ecsvcleanSources[idx_GRBcleanpsf]['e_%sMAG_APER' % band][0]
                if auto_mag_desc is not None:
                    grb_mag_auto = mag_ecsvcleanSources[idx_GRBcleanpsf]['%sMAG_AUTO' % band][0]
                    grb_magerr_auto = mag_ecsvcleanSources[idx_GRBcleanpsf]['e_%sMAG_AUTO' % band][0]

                grb_ra = mag_ecsvcleanSources[idx_GRBcleanpsf]['ALPHA_J2000'][0]
                grb_dec = mag_ecsvcleanSources[idx_GRBcleanpsf]['DELTA_J2000'][0]
                grb_rad = mag_ecsvcleanSources[idx_GRBcleanpsf]['FLUX_RADIUS'][0]
                grb_snr = mag_ecsvcleanSources[idx_GRBcleanpsf]['SNR_WIN'][0]
                grb_elon = mag_ecsvcleanSources[idx_GRBcleanpsf]['ELONGATION'][0]
                grb_dist = d2d[0].to(u.arcsec)
                grb_dist = grb_dist / u.arcsec

                all_magtypes = set(MAGTYPES.keys())
                for mag in sorted(all_magtypes):
                    try:
                        print(f' %s {mag} magnitude of GRB is %.2f +/- %.2f' % (band,
                                                                                mag_ecsvcleanSources[idx_GRBcleanpsf][
                                                                                    f'{band}MAG_{mag}'][0],
                                                                                mag_ecsvcleanSources[idx_GRBcleanpsf][
                                                                                    f'e_{band}MAG_{mag}'][0]))
                    except KeyError:
                        pass

                # survey crsmtch check
                idx_both, idx_bothcleanpsf, d2d_crs, d3d_crs = massCatCoords.search_around_sky(
                    mag_ecsvsourceCatCoords[idx_GRBcleanpsf],
                    grb_rad * u.arcsec)
                if len(idx_bothcleanpsf) > 0:
                    # print(' Detected source crossmatched to existing %s source!' % survey)
                    prime_crs_cat = mag_ecsvcleanSources[idx_GRBcleanpsf]
                    mag_diff_crs_ar, sep, survey_flg_ar = mag_diff_calc(survey_cat=good_cat_stars, prime_cat=prime_crs_cat,
                                                                  survey_idx=idx_bothcleanpsf,
                                                                  prime_idx=[idx_GRBcleanpsf][idx_both],
                                                                  d2d=d2d_crs, band=band)
                else:
                    mag_diff_crs_ar = [99] * mag_col_num
                    survey_flg_ar = [np.int16(0)] * mag_col_num
                    sep = -1

                grbdata['RA'] = Column(np.round(np.array([grb_ra]), 5) * u.deg, description=ra_desc)
                grbdata['DEC'] = Column(np.round(np.array([grb_dec]), decimals=5) * u.deg, description=dec_desc)

                if psf_mag_desc is not None:
                    grbdata['%spsfMag' % band] = Column(np.round(np.array([grb_mag]), 3) * u.ABmag,
                        description=psf_mag_desc
                    )
                    grbdata['%spsfMag_Err' % band] = Column(np.round(np.array([grb_magerr]), 3) * u.ABmag,
                        description=psf_mag_err_desc
                    )
                    grbdata['%spsfME_CM' % band] = Column(
                        np.round(np.array([mag_diff_ar_grabber(mag_diff_crs_ar, survey_flg_ar, 'mag', 'PSF')]),3) * u.ABmag,
                        description=mag_diff_desc
                    )
                    grbdata['psfS_CM'] = Column(mag_diff_ar_grabber(mag_diff_crs_ar, survey_flg_ar, 'flag', 'PSF'),
                                                description=flag_desc)

                if aper_mag_desc is not None:
                    grbdata['%saperMag' % band] = Column(np.round(np.array([grb_mag_aper]), 3) * u.ABmag,
                        description=aper_mag_desc
                    )
                    grbdata['%saperMag_Err' % band] = Column(np.round(np.array([grb_magerr_aper]), 3) * u.ABmag,
                        description=aper_mag_err_desc
                    )
                    grbdata['%saperME_CM' % band] = Column(
                        np.round(np.array([mag_diff_ar_grabber(mag_diff_crs_ar, survey_flg_ar, 'mag', 'APER')]), 3) * u.ABmag,
                        description=mag_diff_desc
                    )
                    grbdata['aperS_CM'] = Column(mag_diff_ar_grabber(mag_diff_crs_ar, survey_flg_ar, 'flag', 'APER'),
                                                 description=flag_desc)

                if auto_mag_desc is not None:
                    grbdata['%sautoMag' % band] = Column(np.round(np.array([grb_mag_auto]), 3) * u.ABmag,
                        description=auto_mag_desc
                    )
                    grbdata['%sautoMag_Err' % band] = Column(np.round(np.array([grb_magerr_auto]), 3) * u.ABmag,
                        description=auto_mag_err_desc
                    )
                    grbdata['%sautoME_CM' % band] = Column(
                        np.round(np.array([mag_diff_ar_grabber(mag_diff_crs_ar, survey_flg_ar, 'mag', 'AUTO')]), 3) * u.ABmag,
                        description=mag_diff_desc
                    )
                    grbdata['autoS_CM'] = Column(mag_diff_ar_grabber(mag_diff_crs_ar, survey_flg_ar, 'flag', 'AUTO'),
                                                 description=flag_desc)

                grbdata['Radius'] = Column(np.round(np.array([grb_rad]), decimals=2) * u.arcsec, description=rad_desc)
                grbdata['SNR'] = Column(np.round(np.array([grb_snr]), decimals=2), description=snr_desc)
                grbdata['Elongation'] = Column(np.round(np.array([grb_elon]), decimals=3), description=elon_desc)
                grbdata['Distance'] = Column(np.round(np.array([grb_dist]), decimals=5) * u.arcsec, description=dist_desc)
                grbdata['Separation'] = Column(np.round(np.array([sep]), decimals=5) * u.arcsec, description=sep_desc)

                grbdata.write('%s_%s_C%i_Data_%s_%s_loc_%d.ecsv' % (grbname, band, chip, survey, num, key),
                              overwrite=True)

                source_reg_gen(grb_ra, grb_dec, rad=grb_rad, append=True)
                print(' Generated GRB data table & source DS9 regions!')
            elif len(idx_GRBcleanpsf) > 1:
                print(' Multiple sources detected in search radius (ra = %s, dec = %s, rad = %s arcsec)'
                      ', refer to .ecsv file for source info!' % (ra, dec, photoDistThresh))
                mag_ar = []
                mag_err_ar = []
                apmag_ar = []
                apmag_err_ar = []
                automag_ar = []
                automag_err_ar = []
                ra_ar = []
                dec_ar = []
                rad_ar = []
                snr_ar = []
                elon_ar = []
                dist_ar = []
                sep_ar = []
                diff_ar = []
                crsmtch_ar = []
                idx_GRBcleanpsflist = idx_GRBcleanpsf.tolist()
                for i in idx_GRBcleanpsflist:
                    if psf_mag_desc is not None:
                        grb_mag = mag_ecsvcleanSources[idx_GRBcleanpsf]['%sMAG_PSF' % band][idx_GRBcleanpsflist.index(i)]
                        mag_ar.append(grb_mag)
                        grb_magerr = mag_ecsvcleanSources[idx_GRBcleanpsf]['e_%sMAG_PSF' % band][
                            idx_GRBcleanpsflist.index(i)]
                        mag_err_ar.append(grb_magerr)
                    if aper_mag_desc is not None:
                        grb_mag_aper = mag_ecsvcleanSources[idx_GRBcleanpsf]['%sMAG_APER' % band][idx_GRBcleanpsflist.index(i)]
                        apmag_ar.append(grb_mag_aper)
                        grb_magerr_aper = mag_ecsvcleanSources[idx_GRBcleanpsf]['e_%sMAG_APER' % band][
                            idx_GRBcleanpsflist.index(i)]
                        apmag_err_ar.append(grb_magerr_aper)
                    if auto_mag_desc is not None:
                        grb_mag_auto = mag_ecsvcleanSources[idx_GRBcleanpsf]['%sMAG_AUTO' % band][
                            idx_GRBcleanpsflist.index(i)]
                        automag_ar.append(grb_mag_auto)
                        grb_magerr_auto = mag_ecsvcleanSources[idx_GRBcleanpsf]['e_%sMAG_AUTO' % band][
                            idx_GRBcleanpsflist.index(i)]
                        automag_err_ar.append(grb_magerr_auto)
                    grb_ra = mag_ecsvcleanSources[idx_GRBcleanpsf]['ALPHA_J2000'][idx_GRBcleanpsflist.index(i)]
                    ra_ar.append(grb_ra)
                    grb_dec = mag_ecsvcleanSources[idx_GRBcleanpsf]['DELTA_J2000'][idx_GRBcleanpsflist.index(i)]
                    dec_ar.append(grb_dec)
                    grb_rad = mag_ecsvcleanSources[idx_GRBcleanpsf]['FLUX_RADIUS'][idx_GRBcleanpsflist.index(i)]
                    rad_ar.append(grb_rad)
                    grb_snr = mag_ecsvcleanSources[idx_GRBcleanpsf]['SNR_WIN'][idx_GRBcleanpsflist.index(i)]
                    snr_ar.append(grb_snr)
                    grb_elon = mag_ecsvcleanSources[idx_GRBcleanpsf]['ELONGATION'][idx_GRBcleanpsflist.index(i)]
                    elon_ar.append(grb_elon)
                    dist = (d2d[idx_GRBcleanpsflist.index(i)]).to(u.arcsec)
                    dist = dist / u.arcsec
                    dist_ar.append(dist)

                    # survey crsmtch check
                    checkcoords = SkyCoord(ra=[mag_ecsvsourceCatCoords[i].ra.deg],
                                           dec=[mag_ecsvsourceCatCoords[i].dec.deg],
                                           frame='icrs', unit='degree')
                    idx_both, idx_bothcleanpsf, d2d_crs, d3d_crs = massCatCoords.search_around_sky(
                        checkcoords, grb_rad * u.arcsec)
                    if len(idx_bothcleanpsf) > 0:
                        # print(' Detected source crossmatched to existing %s source!' % survey)
                        prime_crs_cat = mag_ecsvcleanSources[i]
                        mag_diff_crs_ar, sep, survey_flg_ar = mag_diff_calc(survey_cat=good_cat_stars,
                                                                      prime_cat=prime_crs_cat,
                                                                      survey_idx=idx_bothcleanpsf,
                                                                      prime_idx=[i][idx_both],
                                                                      d2d=d2d_crs, band=band)
                    else:
                        mag_diff_crs_ar = [99] * mag_col_num
                        survey_flg_ar = [np.int16(0)] * mag_col_num
                        sep = -1

                    diff_ar.append(mag_diff_crs_ar)
                    crsmtch_ar.append(survey_flg_ar)
                    sep_ar.append(sep)

                    source_reg_gen(grb_ra, grb_dec, rad=grb_rad, append=True)

                grbdata['RA'] = Column(np.round(np.array(ra_ar), 5) * u.deg, description=ra_desc)
                grbdata['DEC'] = Column(np.round(np.array(dec_ar), decimals=5) * u.deg, description=dec_desc)

                if psf_mag_desc is not None:
                    psf_diff_ar = [mag_diff_ar_grabber(md, sf, 'mag', 'PSF') for md, sf in zip(diff_ar, crsmtch_ar)]
                    psf_flg_ar = [mag_diff_ar_grabber(md, sf, 'flag', 'PSF') for md, sf in zip(diff_ar, crsmtch_ar)]

                    grbdata['%spsfMag' % band] = Column(np.round(np.array(mag_ar), 3) * u.ABmag,
                        description=psf_mag_desc
                    )
                    grbdata['%spsfMag_Err' % band] = Column(np.round(np.array(mag_err_ar), 3) * u.ABmag,
                        description=psf_mag_err_desc
                    )
                    grbdata['%spsfME_CM' % band] = Column(np.round(np.array(psf_diff_ar), 3) * u.ABmag,
                        description=mag_diff_desc
                    )
                    grbdata['psfS_CM'] = Column(np.array(psf_flg_ar), description=flag_desc)

                if aper_mag_desc is not None:
                    aper_diff_ar = [mag_diff_ar_grabber(md, sf, 'mag', 'APER') for md, sf in zip(diff_ar, crsmtch_ar)]
                    aper_flg_ar = [mag_diff_ar_grabber(md, sf, 'flag', 'APER') for md, sf in zip(diff_ar, crsmtch_ar)]

                    grbdata['%saperMag' % band] = Column(np.round(np.array(apmag_ar), 3) * u.ABmag,
                        description=aper_mag_desc
                    )
                    grbdata['%saperMag_Err' % band] = Column(np.round(np.array(apmag_err_ar), 3) * u.ABmag,
                        description=aper_mag_err_desc
                    )
                    grbdata['%saperME_CM' % band] = Column(np.round(np.array(aper_diff_ar), 3) * u.ABmag,
                        description=mag_diff_desc
                    )
                    grbdata['aperS_CM'] = Column(np.array(aper_flg_ar), description=flag_desc)

                if auto_mag_desc is not None:
                    auto_diff_ar = [mag_diff_ar_grabber(md, sf, 'mag', 'AUTO') for md, sf in zip(diff_ar, crsmtch_ar)]
                    auto_flg_ar = [mag_diff_ar_grabber(md, sf, 'flag', 'AUTO') for md, sf in zip(diff_ar, crsmtch_ar)]

                    grbdata['%sautoMag' % band] = Column(np.round(np.array(automag_ar), 3) * u.ABmag,
                        description=auto_mag_desc
                    )
                    grbdata['%sautoMag_Err' % band] = Column(np.round(np.array(automag_err_ar), 3) * u.ABmag,
                        description=auto_mag_err_desc
                    )
                    grbdata['%sautoME_CM' % band] = Column(np.round(np.array(auto_diff_ar), 3) * u.ABmag,
                        description=mag_diff_desc
                    )
                    grbdata['autoS_CM'] = Column(np.array(auto_flg_ar), description=flag_desc)

                grbdata['Radius'] = Column(np.round(np.array(rad_ar), decimals=2) * u.arcsec, description=rad_desc)
                grbdata['SNR'] = Column(np.round(np.array(snr_ar), decimals=2), description=snr_desc)
                grbdata['Elongation'] = Column(np.round(np.array(elon_ar), decimals=3), description=elon_desc)
                grbdata['Distance'] = Column(np.round(np.array(dist_ar), decimals=5) * u.arcsec, description=dist_desc)
                grbdata['Separation'] = Column(np.round(np.array(sep_ar), decimals=5) * u.arcsec, description=sep_desc)

                grbdata.write('%s_Multisource_%s_C%i_Data_%s_%s_loc_%d.ecsv' % (grbname, band, chip, survey, num, key),
                              overwrite=True)
                print(' Generated GRB data table & source DS9 regions!')
            else:
                print(
                    ' GRB source at inputted coords %s and %s not found in PRIME catalog, perhaps increase photoDistThresh or '
                    'alter sextractor params?' % (ra, dec))
    else:
        print(' idx size = %d' % len(idx_GRBcleanpsf))
        if len(idx_GRBcleanpsf) == 1:
            print(' GRB source at inputted coords %s and %s, rad = %s arcsec found!' % (ra, dec, photoDistThresh))

            if psf_mag_desc is not None:
                grb_mag = mag_ecsvcleanSources[idx_GRBcleanpsf]['%sMAG_PSF' % band][0]
                grb_magerr = mag_ecsvcleanSources[idx_GRBcleanpsf]['e_%sMAG_PSF' % band][0]
            if aper_mag_desc is not None:
                grb_mag_aper = mag_ecsvcleanSources[idx_GRBcleanpsf]['%sMAG_APER' % band][0]
                grb_magerr_aper = mag_ecsvcleanSources[idx_GRBcleanpsf]['e_%sMAG_APER' % band][0]
            if auto_mag_desc is not None:
                grb_mag_auto = mag_ecsvcleanSources[idx_GRBcleanpsf]['%sMAG_AUTO' % band][0]
                grb_magerr_auto = mag_ecsvcleanSources[idx_GRBcleanpsf]['e_%sMAG_AUTO' % band][0]

            grb_ra = mag_ecsvcleanSources[idx_GRBcleanpsf]['ALPHA_J2000'][0]
            grb_dec = mag_ecsvcleanSources[idx_GRBcleanpsf]['DELTA_J2000'][0]
            grb_rad = mag_ecsvcleanSources[idx_GRBcleanpsf]['FLUX_RADIUS'][0]
            grb_snr = mag_ecsvcleanSources[idx_GRBcleanpsf]['SNR_WIN'][0]
            grb_elon = mag_ecsvcleanSources[idx_GRBcleanpsf]['ELONGATION'][0]
            grb_dist = d2d[0].to(u.arcsec)
            grb_dist = grb_dist / u.arcsec

            print(
                ' Detected GRB ra = %.6f, dec = %.6f, with 50 percent flux radius (HWHM) = %.3f arcsec and SNR = %.3f' % (
                    grb_ra, grb_dec, grb_rad, grb_snr))

            all_magtypes = set(MAGTYPES.keys())
            for mag in sorted(all_magtypes):
                try:
                    print(f' %s {mag} magnitude of GRB is %.2f +/- %.2f' % (band,
                                                                     mag_ecsvcleanSources[idx_GRBcleanpsf][f'{band}MAG_{mag}'][0],
                                                                     mag_ecsvcleanSources[idx_GRBcleanpsf][f'e_{band}MAG_{mag}'][0]))
                except KeyError:
                    pass

            # survey crsmtch check
            idx_both, idx_bothcleanpsf, d2d_crs, d3d_crs = (
                massCatCoords.search_around_sky(mag_ecsvsourceCatCoords[idx_GRBcleanpsf], grb_rad * u.arcsec))
            if len(idx_bothcleanpsf) > 0:
                # print(' Detected source crossmatched to existing %s source!' % survey)
                prime_crs_cat = mag_ecsvcleanSources[idx_GRBcleanpsf]
                mag_diff_crs_ar, sep, survey_flg_ar = mag_diff_calc(survey_cat=good_cat_stars, prime_cat=prime_crs_cat,
                                                              survey_idx=idx_bothcleanpsf, prime_idx=idx_both,
                                                              d2d=d2d_crs, band=band)
            else:
                mag_diff_crs_ar = [99] * mag_col_num
                survey_flg_ar = [np.int16(0)] * mag_col_num
                sep = -1

            grbdata['RA'] = Column(np.round(np.array([grb_ra]), 5) * u.deg, description=ra_desc)
            grbdata['DEC'] = Column(np.round(np.array([grb_dec]), decimals=5) * u.deg, description=dec_desc)

            if psf_mag_desc is not None:
                grbdata['%spsfMag' % band] = Column(np.round(np.array([grb_mag]), 3) * u.ABmag,
                    description=psf_mag_desc
                )
                grbdata['%spsfMag_Err' % band] = Column(np.round(np.array([grb_magerr]), 3) * u.ABmag,
                    description=psf_mag_err_desc
                )
                grbdata['%spsfME_CM' % band] = Column(
                    np.round(np.array([mag_diff_ar_grabber(mag_diff_crs_ar, survey_flg_ar, 'mag', 'PSF')]),3) * u.ABmag,
                    description=mag_diff_desc
                )
                grbdata['psfS_CM'] = Column(mag_diff_ar_grabber(mag_diff_crs_ar, survey_flg_ar, 'flag', 'PSF'),
                                            description=flag_desc)

            if aper_mag_desc is not None:
                grbdata['%saperMag' % band] = Column(np.round(np.array([grb_mag_aper]), 3) * u.ABmag,
                    description=aper_mag_desc
                )
                grbdata['%saperMag_Err' % band] = Column(np.round(np.array([grb_magerr_aper]), 3) * u.ABmag,
                    description=aper_mag_err_desc
                )
                grbdata['%saperME_CM' % band] = Column(
                    np.round(np.array([mag_diff_ar_grabber(mag_diff_crs_ar, survey_flg_ar, 'mag', 'APER')]), 3) * u.ABmag,
                    description=mag_diff_desc
                )
                grbdata['aperS_CM'] = Column(mag_diff_ar_grabber(mag_diff_crs_ar, survey_flg_ar, 'flag', 'APER'),
                                             description=flag_desc)

            if auto_mag_desc is not None:
                grbdata['%sautoMag' % band] = Column(np.round(np.array([grb_mag_auto]), 3) * u.ABmag,
                    description=auto_mag_desc
                )
                grbdata['%sautoMag_Err' % band] = Column(np.round(np.array([grb_magerr_auto]), 3) * u.ABmag,
                    description=auto_mag_err_desc
                )
                grbdata['%sautoME_CM' % band] = Column(
                    np.round(np.array([mag_diff_ar_grabber(mag_diff_crs_ar, survey_flg_ar, 'mag', 'AUTO')]), 3) * u.ABmag,
                    description=mag_diff_desc
                )
                grbdata['autoS_CM'] = Column(mag_diff_ar_grabber(mag_diff_crs_ar, survey_flg_ar, 'flag', 'AUTO'),
                                             description=flag_desc)

            grbdata['Radius'] = Column(np.round(np.array([grb_rad]), decimals=2) * u.arcsec, description=rad_desc)
            grbdata['SNR'] = Column(np.round(np.array([grb_snr]), decimals=2), description=snr_desc)
            grbdata['Elongation'] = Column(np.round(np.array([grb_elon]), decimals=3), description=elon_desc)
            grbdata['Distance'] = Column(np.round(np.array([grb_dist]), decimals=5) * u.arcsec, description=dist_desc)
            grbdata['Separation'] = Column(np.round(np.array([sep]), decimals=5) * u.arcsec, description=sep_desc)

            grbdata.write('%s_%s_C%i_Data_%s_%s.ecsv' % (grbname, band, chip, survey, num), overwrite=True)

            regprimename = source_reg_gen(grb_ra, grb_dec, rad=grb_rad)

            savename, threshname = grb_cutout(imageName, GRBcoords, photoDistThresh,
                                              regprimename=regprimename, regsurvname=regsurvname)

            html_gen(grbdata, directory, savename, threshname, band, survey, ra, dec, thresh, regsurvname, regprimename)

            print(' Generated GRB data table & source DS9 regions!')
        elif len(idx_GRBcleanpsf) > 1:
            print(' Multiple sources detected in search radius (ra = %s, dec = %s, rad = %s arcsec)'
                  ', refer to .ecsv file for source info!' % (ra, dec, photoDistThresh))
            mag_ar = []
            mag_err_ar = []
            apmag_ar = []
            apmag_err_ar = []
            automag_ar = []
            automag_err_ar = []
            ra_ar = []
            dec_ar = []
            rad_ar = []
            snr_ar = []
            elon_ar = []
            dist_ar = []
            sep_ar = []
            diff_ar = []
            crsmtch_ar = []
            idx_GRBcleanpsflist = idx_GRBcleanpsf.tolist()
            source_reg_gen()
            for i in idx_GRBcleanpsflist:
                if psf_mag_desc is not None:
                    grb_mag = mag_ecsvcleanSources[idx_GRBcleanpsf]['%sMAG_PSF' % band][idx_GRBcleanpsflist.index(i)]
                    mag_ar.append(grb_mag)
                    grb_magerr = mag_ecsvcleanSources[idx_GRBcleanpsf]['e_%sMAG_PSF' % band][idx_GRBcleanpsflist.index(i)]
                    mag_err_ar.append(grb_magerr)
                if aper_mag_desc is not None:
                    grb_mag_aper = mag_ecsvcleanSources[idx_GRBcleanpsf]['%sMAG_APER' % band][idx_GRBcleanpsflist.index(i)]
                    apmag_ar.append(grb_mag_aper)
                    grb_magerr_aper = mag_ecsvcleanSources[idx_GRBcleanpsf]['e_%sMAG_APER' % band][idx_GRBcleanpsflist.index(i)]
                    apmag_err_ar.append(grb_magerr_aper)
                    # grb_mag_aper2 = mag_ecsvcleanSources[idx_GRBcleanpsf]['%sMAG_APER' % band][:, 0][
                    #     idx_GRBcleanpsflist.index(i)]
                if auto_mag_desc is not None:
                    grb_mag_auto = mag_ecsvcleanSources[idx_GRBcleanpsf]['%sMAG_AUTO' % band][idx_GRBcleanpsflist.index(i)]
                    automag_ar.append(grb_mag_auto)
                    grb_magerr_auto = mag_ecsvcleanSources[idx_GRBcleanpsf]['e_%sMAG_AUTO' % band][
                        idx_GRBcleanpsflist.index(i)]
                    automag_err_ar.append(grb_magerr_auto)
                grb_ra = mag_ecsvcleanSources[idx_GRBcleanpsf]['ALPHA_J2000'][idx_GRBcleanpsflist.index(i)]
                ra_ar.append(grb_ra)
                grb_dec = mag_ecsvcleanSources[idx_GRBcleanpsf]['DELTA_J2000'][idx_GRBcleanpsflist.index(i)]
                dec_ar.append(grb_dec)
                grb_rad = mag_ecsvcleanSources[idx_GRBcleanpsf]['FLUX_RADIUS'][idx_GRBcleanpsflist.index(i)]
                rad_ar.append(grb_rad)
                grb_snr = mag_ecsvcleanSources[idx_GRBcleanpsf]['SNR_WIN'][idx_GRBcleanpsflist.index(i)]
                snr_ar.append(grb_snr)
                grb_elon = mag_ecsvcleanSources[idx_GRBcleanpsf]['ELONGATION'][idx_GRBcleanpsflist.index(i)]
                elon_ar.append(grb_elon)
                dist = (d2d[idx_GRBcleanpsflist.index(i)]).to(u.arcsec)
                dist = dist / u.arcsec
                dist_ar.append(dist)

                # survey crsmtch check
                checkcoords = SkyCoord(ra=[mag_ecsvsourceCatCoords[i].ra.deg], dec=[mag_ecsvsourceCatCoords[i].dec.deg],
                                       frame='icrs', unit='degree')
                idx_both, idx_bothcleanpsf, d2d_crs, d3d_crs = massCatCoords.search_around_sky(
                    checkcoords, grb_rad * u.arcsec)

                if len(idx_bothcleanpsf) > 0:
                    # print(' Detected source crossmatched to existing %s source!' % survey)
                    prime_crs_cat = mag_ecsvcleanSources[i]
                    mag_diff_crs_ar, sep, survey_flg_ar = mag_diff_calc(survey_cat=good_cat_stars, prime_cat=prime_crs_cat,
                                                                  survey_idx=idx_bothcleanpsf, prime_idx=idx_both,
                                                                  d2d=d2d_crs, band=band)
                else:
                    mag_diff_crs_ar = [99] * mag_col_num
                    survey_flg_ar = [np.int16(0)] * mag_col_num
                    sep = -1

                diff_ar.append(mag_diff_crs_ar)
                crsmtch_ar.append(survey_flg_ar)
                sep_ar.append(sep)

                regprimename = source_reg_gen(grb_ra, grb_dec, rad=grb_rad, append=True)

            grbdata['RA'] = Column(np.round(np.array(ra_ar), 5) * u.deg, description=ra_desc)
            grbdata['DEC'] = Column(np.round(np.array(dec_ar), decimals=5) * u.deg, description=dec_desc)

            if psf_mag_desc is not None:
                psf_diff_ar = [mag_diff_ar_grabber(md, sf, 'mag', 'PSF') for md, sf in zip(diff_ar, crsmtch_ar)]
                psf_flg_ar = [mag_diff_ar_grabber(md, sf, 'flag', 'PSF') for md, sf in zip(diff_ar, crsmtch_ar)]

                grbdata['%spsfMag' % band] = Column(np.round(np.array(mag_ar), 3) * u.ABmag,
                    description=psf_mag_desc
                )
                grbdata['%spsfMag_Err' % band] = Column(np.round(np.array(mag_err_ar), 3) * u.ABmag,
                    description=psf_mag_err_desc
                )
                grbdata['%spsfME_CM' % band] = Column(np.round(np.array(psf_diff_ar), 3) * u.ABmag,
                    description=mag_diff_desc
                )
                grbdata['psfS_CM'] = Column(np.array(psf_flg_ar), description=flag_desc)

            if aper_mag_desc is not None:
                aper_diff_ar = [mag_diff_ar_grabber(md, sf, 'mag', 'APER') for md, sf in zip(diff_ar, crsmtch_ar)]
                aper_flg_ar = [mag_diff_ar_grabber(md, sf, 'flag', 'APER') for md, sf in zip(diff_ar, crsmtch_ar)]

                grbdata['%saperMag' % band] = Column(np.round(np.array(apmag_ar), 3) * u.ABmag,
                    description=aper_mag_desc
                )
                grbdata['%saperMag_Err' % band] = Column(np.round(np.array(apmag_err_ar), 3) * u.ABmag,
                    description=aper_mag_err_desc
                )
                grbdata['%saperME_CM' % band] = Column(np.round(np.array(aper_diff_ar), 3) * u.ABmag,
                    description=mag_diff_desc
                )
                grbdata['aperS_CM'] = Column(np.array(aper_flg_ar), description=flag_desc)

            if auto_mag_desc is not None:
                auto_diff_ar = [mag_diff_ar_grabber(md, sf, 'mag', 'AUTO') for md, sf in zip(diff_ar, crsmtch_ar)]
                auto_flg_ar = [mag_diff_ar_grabber(md, sf, 'flag', 'AUTO') for md, sf in zip(diff_ar, crsmtch_ar)]

                grbdata['%sautoMag' % band] = Column(np.round(np.array(automag_ar), 3) * u.ABmag,
                    description=auto_mag_desc
                )
                grbdata['%sautoMag_Err' % band] = Column(np.round(np.array(automag_err_ar), 3) * u.ABmag,
                    description=auto_mag_err_desc
                )
                grbdata['%sautoME_CM' % band] = Column(np.round(np.array(auto_diff_ar), 3) * u.ABmag,
                    description=mag_diff_desc
                )
                grbdata['autoS_CM'] = Column(np.array(auto_flg_ar), description=flag_desc)

            grbdata['Radius'] = Column(np.round(np.array(rad_ar), decimals=2) * u.arcsec, description=rad_desc)
            grbdata['SNR'] = Column(np.round(np.array(snr_ar), decimals=2), description=snr_desc)
            grbdata['Elongation'] = Column(np.round(np.array(elon_ar), decimals=3), description=elon_desc)
            grbdata['Distance'] = Column(np.round(np.array(dist_ar), decimals=5) * u.arcsec, description=dist_desc)
            grbdata['Separation'] = Column(np.round(np.array(sep_ar), decimals=5) * u.arcsec, description=sep_desc)

            grbdata.write('%s_Multisource_%s_C%i_Data_%s_%s.ecsv' % (grbname, band, chip, survey, num), overwrite=True)

            savename, threshname = grb_cutout(imageName, GRBcoords, photoDistThresh,
                                              regprimename=regprimename, regsurvname=regsurvname)

            html_gen(grbdata, directory, savename, threshname, band, survey, ra, dec, thresh, regsurvname, regprimename)

            print(' Generated GRB data table & source DS9 regions!')
        else:
            print(
                ' GRB source at inputted coords %s and %s not found in PRIME catalog, perhaps increase photoDistThresh or '
                'alter sextractor params?' % (ra, dec))


# WLS fit calculations

def single_fit_calc(cleanPSFsources, band, good_cat_stars, idx_psfmass, idx_psfimage,
                         weights_noclip, clipped, sigma, with_plots=False, magtype=defaults['magtype']):
    """Calculates WLS fit line used for determining the goodness of the photometry, for 1 magtype"""

    # appropriate mag column
    colnames = good_cat_stars.colnames
    if len(colnames) > 6:
        magcol = f'{band}MAG_{magtype}'
        magerrcol = f'{band}MAG_{magtype}'
    else:
        magcol = colnames[2]
        magerrcol = colnames[3]

    # sigma residual fit
    x = good_cat_stars['%s' % magcol][idx_psfmass][~clipped.mask]
    y = cleanPSFsources[f'{band}MAG_{MAGTYPES[magtype]}'][idx_psfimage][~clipped.mask]

    x_const = sm.add_constant(x)
    model = sm.WLS(y, x_const, weights=weights_noclip[~clipped.mask]).fit()

    x_nc = good_cat_stars['%s' % magcol][idx_psfmass]
    y_nc = cleanPSFsources[f'{band}MAG_{MAGTYPES[magtype]}'][idx_psfimage]
    x_const_nc = sm.add_constant(x_nc)
    model2 = sm.WLS(y_nc, x_const_nc, weights=weights_noclip).fit()

    m = model.params[1]
    m_err = model.bse[1]
    b = model.params[0]
    b_err = model.bse[0]
    model_resid = model.resid

    print('Num of crossmatched sources used in %s %s sig fit: %i'
          % (MAGTYPES[magtype], sigma, len(cleanPSFsources[f'{band}MAG_{MAGTYPES[magtype]}'][idx_psfimage][~clipped.mask])))
    print(' %s sig fit: slope = %.4f +/- %.4f' % (sigma, m, m_err))
    print(' %s sig fit: y-int = %.4f +/- %.4f' % (sigma, b, b_err))

    if with_plots:
        return model, model2
    else:
        full_b_err = round(3 * b_err, 4)
        if full_b_err >= 0.25:
            print(f'3 sig int. error larger than expected.. {full_b_err}')
            print('Reducing error (3 sig -> 0.5 sig) threshold significantly, should investigate image / results')
            return m, b, round(0.5 * b_err, 4)
        return m, b, full_b_err


def photometric_fit_calc(cleanPSFsources, band, good_cat_stars, idx_psfmass, idx_psfimage,
                         psfweights_noclip, psf_clipped, sigma, aperweights_noclip, aper_clipped_all,
                         autoweights_noclip, auto_clipped, with_plots=False, magtype=defaults['magtype']):
    """Calculates WLS fit lines used for determining the goodness of the photometry"""

    # appropriate mag column
    colnames = good_cat_stars.colnames
    if len(colnames) > 6:
        magcol = f'{band}MAG_{magtype}'
        magerrcol = f'{band}MAG_{magtype}'
    else:
        magcol = colnames[2]
        magerrcol = colnames[3]

    # sigma residual fit for auto aperture photometry
    x_auto = good_cat_stars['%s' % magcol][idx_psfmass][~auto_clipped.mask]
    y_auto = cleanPSFsources['%sMAG_AUTO' % band][idx_psfimage][~auto_clipped.mask]
    x_const_auto = sm.add_constant(x_auto)
    model_auto = sm.WLS(y_auto, x_const_auto, weights=autoweights_noclip[~auto_clipped.mask]).fit()

    x_auto_nc = good_cat_stars['%s' % magcol][idx_psfmass]
    y_auto_nc = cleanPSFsources['%sMAG_AUTO' % band][idx_psfimage]
    x_const_auto_nc = sm.add_constant(x_auto_nc)
    model2 = sm.WLS(y_auto_nc, x_const_auto_nc, weights=autoweights_noclip).fit()

    m_auto = model_auto.params[1]
    m_autoerr = model_auto.bse[1]
    b_auto = model_auto.params[0]
    b_autoerr = model_auto.bse[0]
    model_auto_resid = model_auto.resid

    print('Num of crossmatched sources used in auto aperture %s sig fit: %i'
          % (sigma, len(cleanPSFsources['%sMAG_AUTO' % band][idx_psfimage][~auto_clipped.mask])))
    print(' %s sig fit: slope = %.4f +/- %.4f' % (sigma, m_auto, m_autoerr))
    print(' %s sig fit: y-int = %.4f +/- %.4f' % (sigma, b_auto, b_autoerr))

    aper_model_sigs = model_sig = model_sig_resid = 0
    aperweights_noclip = []

    if len(psfweights_noclip) > 0:
        model2 = 0
        x = cleanPSFsources['%sMAG_PSF' % band][idx_psfimage]
        y = good_cat_stars['%s' % magcol][idx_psfmass]
        x_const = sm.add_constant(x)
        # model = sm.WLS(y, x, weights=psfweights).fit()
        model2 = sm.WLS(y, x_const, weights=psfweights_noclip).fit()

        avg2 = np.average(model2.resid, weights=psfweights_noclip)
        var2 = np.average((model2.resid - avg2) ** 2, weights=psfweights_noclip)

        # residual fit - 3 sigma clip
        x_sig = good_cat_stars['%s' % magcol][idx_psfmass][~psf_clipped.mask]
        y_sig = cleanPSFsources['%sMAG_PSF' % band][idx_psfimage][~psf_clipped.mask]
        x_const_sig = sm.add_constant(x_sig)
        model_sig = sm.WLS(y_sig, x_const_sig, weights=psfweights_noclip[~psf_clipped.mask]).fit()
        m_sig = model_sig.params[1]
        m_sigerr = model_sig.bse[1]
        b_sig = model_sig.params[0]
        b_sigerr = model_sig.bse[0]
        model_sig_resid = model_sig.resid

        print('Num of crossmatched sources used in PSF %s sig fit: %i'
              % (sigma, len(cleanPSFsources['%sMAG_PSF' % band][idx_psfimage][~psf_clipped.mask])))
        print(' %s sig fit: slope = %.4f +/- %.4f' % (sigma, m_sig, m_sigerr))
        print(' %s sig fit: y-int = %.4f +/- %.4f' % (sigma, b_sig, b_sigerr))

    # sigma clipped resids for aper mags, disabled right now
    if len(aperweights_noclip) > 0:
        aperweights_noclip = aperweights_noclip.tolist()
        aper_model_sigs= []
        for idx, (aperweight_nc, aperclipped) in enumerate(zip(aperweights_noclip, aper_clipped_all)):
            x_sig_ap = good_cat_stars['%s' % magcol][idx_psfmass][~aperclipped.mask]
            y_sig_ap = cleanPSFsources['%sMAG_APER' % band][:,idx][idx_psfimage][~aperclipped.mask]
            x_const_sig_ap = sm.add_constant(x_sig_ap)
            model_sig_ap = sm.WLS(y_sig_ap, x_const_sig_ap, weights=aperweight_nc[~aperclipped.mask]).fit()
            aper_model_sigs.append(model_sig_ap)

    if with_plots:
        return model2, model_sig, model_sig_resid, model_auto, model_auto_resid, aper_model_sigs
    else:
        if magtype == 'AUTO':
            return m_auto, b_auto, round(3 * b_autoerr, 4)
        elif magtype == 'PSF':
            b_sig_err = round(3 * b_sigerr, 4)
            if b_sig_err > 0.1:
                print(' PSF intercept err > 0.1, reducing acc. err value from 3 sigma to 1.5')
                return m_sig, b_sig, round(1.5 * b_sigerr, 4)
            else:
                return m_sig, b_sig, round(3 * b_sigerr, 4)


# %% optional plots


# single photometry plots for int calibration

def single_plots(cleanPSFsources, PSFsources, data, imageName, survey, band, good_cat_stars, idx_mass,
                     idx_image, sigma, weights_noclip, clipped, crop):
    """
    For use in multiprocessing calibration, generates plots given any single magtype.
    """

    # appropriate mag column
    colnames = good_cat_stars.colnames
    if len(colnames) > 6:
        magcol = f'{band}MAG_{magtype}'
        magerrcol = f'{band}MAG_{magtype}'
    else:
        magcol = colnames[2]
        magerrcol = colnames[3]

    chip = imageName[-6]
    if len(imageName) <= 16:
        num = 'img'
    else:
        num = imageName[-16:-8]

    # PSF-specific naming for parallel running
    end_name = end_name_gen('png')

    def predict_y_for(x, m, b):
        return m * x + b

    # PHOTOMETRIC FIT LINE MODELS

    model, model2 = single_fit_calc(cleanPSFsources, band, good_cat_stars, idx_mass, idx_image,
                         weights_noclip, clipped, sigma, with_plots=True, magtype=magtype)

    m = model.params[1]
    m_err = model.bse[1]
    b = model.params[0]
    b_err = model.bse[0]
    rsquare = model.rsquared
    rss = model.ssr
    model_resid = model.resid

    m2 = model2.params[1]
    m2err = model2.bse[1]
    b2 = model2.params[0]
    b2err = model2.bse[0]
    avg2 = np.average(model2.resid, weights=weights_noclip)

    # BINNED RESIDUAL STATISTICS

    # residual fit 3 sig clip - bin errors and stats
    mags = range(10, 22)
    x_arr = np.arange(10 + 0.5, 22 + 0.5, 1)

    res_errs = []
    res_means = []
    res_ranges = []
    res_skews = []
    res_nums = []
    for i in mags:
        mask = ((cleanPSFsources[f'{band}MAG_{MAGTYPES[magtype]}'][idx_image][~clipped.mask] > i) &
                     (cleanPSFsources[f'{band}MAG_{MAGTYPES[magtype]}'][idx_image][~clipped.mask] < i + 1))
        mask = model_resid[mask]
        res_err = np.std(mask)  # spread
        res_range = '%s - %s' % (i, i + 1)
        res_mean = np.nanmean(mask)  # mean
        res_skew = skew(mask, bias=False)  # skew
        res_num = len(mask)
        if ma.is_masked(res_err):
            res_err = 0
        res_errs.append(res_err)
        res_means.append(res_mean)
        res_ranges.append(res_range)
        res_skews.append(res_skew)
        res_nums.append(res_num)
    res_errs = np.nan_to_num(np.array(res_errs))
    res_means = np.nan_to_num(np.array(res_means))
    res_skews = np.nan_to_num(np.array(res_skews))
    res_errs_min = np.min((res_errs[res_errs != 0]))
    res_errs_max = np.max(res_errs)

    bintable = Table()
    bintable['Bin Range (mag)'] = res_ranges
    bintable['Mean (mag)'] = res_means
    bintable['Spread (mag)'] = res_errs
    bintable['Skew'] = res_skews
    bintable['Source Number'] = res_nums
    bintable.write('Resid_%s-sig_Data_%s_C%s_%s%s' % (sigma, band, chip, survey, end_name_gen()), overwrite=True)

    # ALL PLOTS

    plt.close('all')
    # mag comparison plot
    plt.figure(1, figsize=(8, 8))
    plt.plot(cleanPSFsources[f'{band}MAG_{MAGTYPES[magtype]}'][idx_image], good_cat_stars['%s' % magcol][idx_mass],
             'r.', markersize=14, markeredgecolor='black')
    plt.xlim(10, 22)
    plt.ylim(10, 22)
    plt.title('PRIME Mags vs %s Mags' % survey)
    plt.xlabel('PRIME %s Mags' % band, fontsize=15)
    plt.ylabel('%s %s Mags' % (survey, band), fontsize=15)
    plt.grid()
    # plt.savefig('%s_C%s_mag_comp_plot_%s.png' % (survey, chip, num))
    plt.clf()
    # print('Saved mag comparison plot to dir!')

    # PRIME flux vs catalog AB mag for crossmatches
    x_flx_lin = cleanPSFsources[f'{MAGTYPES[magtype]}_FLUX_DENSITY'][idx_image][~clipped.mask]
    x_flx = np.log(x_flx_lin)
    y_flx = good_cat_stars['%s' % magcol][idx_mass][~clipped.mask]
    x_const_flx = sm.add_constant(x_flx)
    model_flx = sm.WLS(y_flx, x_const_flx, weights=weights_noclip[~clipped.mask]).fit()
    # print(model_flx.params)
    m_flx = model_flx.params[1]
    m_flxerr = model_flx.bse[1]
    b_flx = model_flx.params[0]
    b_flxerr = model_flx.bse[0]

    x_grid = np.linspace(x_flx_lin.min(), x_flx_lin.max(), 100)
    log_x_grid = np.log(x_grid)
    log_x_grid_const = sm.add_constant(log_x_grid)
    y_pred = model_flx.predict(log_x_grid_const)

    plt.figure(9, figsize=(8, 8))
    plt.plot(cleanPSFsources[f'{MAGTYPES[magtype]}_FLUX_DENSITY'][idx_image][~clipped.mask],
             good_cat_stars['%s' % magcol][idx_mass][~clipped.mask],
             'r.', markersize=14, markeredgecolor='black')
    plt.plot(x_grid, y_pred, c='b')

    flx_txt = ('slope = %.4f' % m_flx + '\nslope err = %.4f' % m_flxerr +
               '\nint = %.4f' % b_flx + '\nint err = %.4f' % b_flxerr)

    plt.xlim(10, 50000)
    plt.ylim(12, 20.5)
    plt.title('PRIME Flux Density vs %s AB mag - %s Sigma Clip' % (survey, sigma))
    plt.xlabel(r'PRIME Flux Density ($\mu$Jy)', fontsize=15)
    plt.ylabel('%s %s AB Mags' % (survey, band), fontsize=15)
    plt.grid()
    plt.xscale('log')
    flx_box = dict(facecolor='white')
    plt.text(1000, 19, flx_txt, fontsize=12, bbox=flx_box)
    plt.savefig('%s_C%s_flux_mag_plot_sig_%s%s' % (survey, chip, num, end_name))
    plt.clf()
    print('Saved flux v. mag plot to dir!')

    # res plot y int, histogram
    if len(idx_image) >= 75000:
        bin_num_int = round(len(idx_image) / 500)
    elif 50000 <= len(idx_image) <= 75000:
        bin_num_int = round(len(idx_image) / 400)
    elif 25000 <= len(idx_image) <= 50000:
        bin_num_int = round(len(idx_image) / 300)
    elif 5000 <= len(idx_image) <= 25000:
        bin_num_int = round(len(idx_image) / 75)
    elif 1000 <= len(idx_image) <= 5000:
        bin_num_int = round(len(idx_image) / 50)
    elif len(idx_image) <= 1000:
        bin_num_int = 50

    # magtype specific coloring, etc.
    alpha = 0.7 if model else None
    cmap = 'gist_heat_r'
    if magtype == MAGTYPES['PSF']:
        cmap = 'gist_heat_r'
        lm_colors = ['red', 'blue']
    elif magtype == MAGTYPES['AUTO']:
        cmap = 'gist_earth_r'
        lm_colors = ['blue', 'black']
    elif magtype == MAGTYPES['APER']:
        cmap = 'pink_r'
        lm_colors = ['green', 'black']


    fig, ax2 = plt.subplots(1, 1, figsize=(9, 8))
    fig.suptitle(f'{survey} Residuals - {sigma} Sigma Clip - Density Histogram')

    hist = ax2.hist2d(
        x=good_cat_stars['%s' % magcol][idx_mass][~clipped.mask], y=model_resid,
        bins=[bin_num_int, bin_num_int], range=[[10, 21], [-1, 1]],
        cmap=cmap
    )

    cbar2 = fig.colorbar(hist[3], ax=ax2, pad=0.03)
    cbar2.set_label(f'{MAGTYPES[magtype]} Density')

    sig = ax2.scatter(x_arr, res_errs, marker='_', s=1625, c='black',
                          label=fr'{MAGTYPES[magtype]} ap. photom 1 $\sigma$ range = [%.3f - %.3f]' % (
                              res_errs_min, res_errs_max))
    ax2.scatter(x_arr, -res_errs, marker='_', s=1625, c='black')

    ax2.axhline(y=0, color='black', linestyle='--', linewidth=1)
    ax2.set_xlim(10, 21)
    ax2.set_ylim(-1, 1)
    ax2.yaxis.set_tick_params(labelleft=True)
    ax2.set_title(f"{MAGTYPES[magtype]} Fit Residuals")
    ax2.set_xlabel(f"{survey} {band} Mags")
    ax2.set_ylabel("Residuals")

    info = (
            f'{MAGTYPES[magtype]} fit'
            + '\nslope = %.4f +/- %.4f' % (m, m_err)
            + '\nintercept = %.3f +/- %.3f' % (b, b_err)
            + '\nR$^{2}$ = %.3f' % rsquare
            + '\nRSS = %d' % rss
            + '\nn_sources = %i' % len(cleanPSFsources[f'{MAGTYPES[magtype]}_FLUX_DENSITY'][idx_image][~clipped.mask])
    )

    ax2.text(10.5, 0.5, info, fontsize=9,
             bbox=dict(facecolor='white', edgecolor='black', pad=5.0))

    ax2.legend([sig], [sig.get_label()], loc='lower left', markerscale=0.5)

    plt.savefig('%s_C%s_residual_plot_int_hist_%s%s' % (survey, chip, num, end_name), dpi=300)
    plt.close(fig)

    print('Saved y-int residual plots to dir!')

    # WLS fit line over data plot

    txt = ('slope = %.4f' % m2 + '\nslope err = %.4f' % m2err + '\nint = %.4f' % b2 + '\nint err = %.4f' % b2err +
           '\nn_sources = %i' % len(cleanPSFsources[idx_image]))

    # WLS hist density plot
    if len(idx_image) >= 75000:
        bin_num = round(len(idx_image) / 500)
    elif 50000 <= len(idx_image) <= 75000:
        bin_num = round(len(idx_image) / 350)
    elif 25000 <= len(idx_image) <= 50000:
        bin_num = round(len(idx_image) / 150)
    elif 5000 <= len(idx_image) <= 25000:
        bin_num = round(len(idx_image) / 50)
    elif 1000 <= len(idx_image) <= 5000:
        bin_num = round(len(idx_image) / 20)
    elif len(idx_image) <= 1000:
        bin_num = 100

    plt.figure(5, figsize=(10, 8))
    plt.clf()
    plt.xlim(10, 22)
    plt.ylim(10, 22)
    plt.title(f'{survey} vs PRIME {MAGTYPES[magtype]} Mag w/ Weighted Fit - Density Histogram')
    plt.grid()
    plt.ylabel('PRIME %s Mags' % band, fontsize=15)
    plt.xlabel('%s %s Mags' % (survey, band), fontsize=15)
    plt.hist2d(x=good_cat_stars['%s' % magcol][idx_mass], y=cleanPSFsources[f'{band}MAG_{MAGTYPES[magtype]}'][idx_image],
               bins=[bin_num, bin_num], range=[[10, 22], [10, 22]], cmap=cmap)
    plt.plot(good_cat_stars['%s' % magcol][idx_mass],
             predict_y_for(good_cat_stars['%s' % magcol][idx_mass], m2, b2), c='b')
    plt.colorbar(label='Density')
    box = dict(facecolor='white')
    plt.text(11, 18, txt, fontsize=12, bbox=box)
    plt.savefig('%s_C%s_WLS_fit_hist_plot_%s%s' % (survey, chip, num, end_name), dpi=300)
    plt.clf()

    print('Saved WLS fit plots to dir!')

    # WLS 3 sig hist density plot

    fig, ax2 = plt.subplots(1, 1, figsize=(9, 8))
    fig.suptitle(f'{survey} vs PRIME w/ Weighted Fit - {sigma} Sigma Clip - Density Histogram')

    hist = ax2.hist2d(
        x=good_cat_stars['%s' % magcol][idx_mass][~clipped.mask],
        y=cleanPSFsources[f'{band}MAG_{MAGTYPES[magtype]}'][idx_image][~clipped.mask],
        bins=[bin_num, bin_num], range=[[10, 22], [10, 22]],
        cmap=cmap
    )

    line = ax2.plot(good_cat_stars['%s' % magcol][idx_mass][~clipped.mask],
                        predict_y_for(good_cat_stars['%s' % magcol][idx_mass][~clipped.mask], m, b),
                        c='r')

    cbar2 = fig.colorbar(hist[3], ax=ax2, pad=0.03)
    cbar2.set_label(f'{MAGTYPES[magtype]} Density')

    ax2.grid()
    ax2.set_xlim(10, 22)
    ax2.set_ylim(10, 22)
    ax2.yaxis.set_tick_params(labelleft=True)
    ax2.set_title(f"{survey} vs PRIME {MAGTYPES[magtype]} Plot w/ WLS fit line")
    ax2.set_xlabel(f"{survey} {band} Mags")
    ax2.set_ylabel(f"PRIME {band} Mags")

    ax2.text(11, 18, info, fontsize=12,
             bbox=dict(facecolor='white', edgecolor='black'))

    plt.savefig('%s_C%s_WLS_fit_3sig_hist_plot_%s%s' % (survey, chip, num, end_name), dpi=300)
    plt.close(fig)

    print('Saved WLS 3 sig fit plots to dir!')

    # flux vs mag plot - histogram vers.
    plt.figure(10, figsize=(8, 8))
    plt.hist2d(x=cleanPSFsources[f'{MAGTYPES[magtype]}_FLUX_DENSITY'][idx_image], y=good_cat_stars['%s' % magcol][idx_mass],
               bins=[bin_num, bin_num], range=[[100, 50000], [12, 20.5]], cmap=cmap)
    plt.colorbar(label='Density')
    plt.xlim(10, 50000)
    plt.ylim(12, 20.5)
    plt.title('PRIME Flux Density vs %s AB mag - Histogram' % survey)
    plt.xlabel(r'PRIME Flux Density ($\mu$Jy)', fontsize=15)
    plt.ylabel('%s %s AB Mags' % (survey, band), fontsize=15)
    plt.grid()
    plt.xscale('log')
    # plt.savefig('%s_C%s_flux_mag_hist_plot_%s.png' % (survey, chip, num))
    plt.clf()

    # print('Saved flux v. mag hist plot to dir!')

    # Limiting Mag Plot

    def lim_mag_calc(bin_vals, all_mags):
        idxs = np.digitize(all_mags, bins=bin_vals)

        indices = {i: [] for i in range(len(bin_vals))}
        for idx, value in enumerate(idxs):
            indices[value].append(idx)
        sorted_indices_lists = list(indices.values())

        all_sources = []
        for i in sorted_indices_lists:
            number = len(i)
            all_sources.append(number)

        # plotting
        idx_arr = np.where(np.isclose(bin_vals, 12.5))
        min_x_idx = idx_arr[0][0]  # avoid saturated <12.5 mag sources from affecting maximum

        filtered_sources = all_sources[min_x_idx:]
        idxmax = filtered_sources.index(max(filtered_sources)) + min_x_idx
        split_sources = all_sources[idxmax:]

        halfmax = max(filtered_sources) / 2
        halfmaxpt = list(max(enumerate(split_sources), key=lambda x: -abs(halfmax - x[1])))
        halfmaxpt = [halfmaxpt[0] + idxmax, halfmaxpt[1]]

        limmag = round(bin_vals[halfmaxpt[0]], 1)

        return all_sources, halfmax, limmag

    max_x = data.shape[0]
    max_y = data.shape[1]
    all_mags_all = PSFsources[(PSFsources[f'{band}MAG_{MAGTYPES[magtype]}'] < 25) &
                              (PSFsources['XWIN_IMAGE'] < (max_x - crop)) & (PSFsources['XWIN_IMAGE'] > crop) &
                              (PSFsources['YWIN_IMAGE'] < (max_y) - crop) & (PSFsources['YWIN_IMAGE'] > crop)
                              ]
    all_mags = all_mags_all[f'{band}MAG_{MAGTYPES[magtype]}']

    bin_vals = np.array(np.arange(12, 25.5, 0.1))

    all_sources, halfmax, limmag = lim_mag_calc(bin_vals=bin_vals, all_mags=all_mags)

    print(f'{MAGTYPES[magtype]} Lim Mag = ', limmag)

    plt.figure(8, figsize=(24, 8))

    lmdata = plt.bar(bin_vals, height=all_sources, width=0.1, align='edge', color=lm_colors[0], edgecolor='black',
                       alpha=1, label=f'{MAGTYPES[magtype]} Binned Sources')
    lmhalf = plt.axhline(halfmax, linestyle='--', color=lm_colors[1],
                           label=f'{MAGTYPES[magtype]} Half Max = %s' % round(halfmax, 1))
    limmag_line = plt.axvline(limmag, color=lm_colors[1], linewidth=2,
                             label=f'{MAGTYPES[magtype]} Limiting Mag = %s' % round(limmag, 1))

    xticks = np.arange(12, 25.5, 0.5)
    plt.xticks(xticks, fontsize=10)
    plt.grid()
    plt.yscale('log')
    plt.title('PRIME Limiting Mag Plot')
    plt.ylabel('Number of Sources')
    plt.xlabel('%s Magnitude' % band)

    handles = [lmdata, lmhalf, limmag_line]
    labels = [h.get_label() for h in handles]
    plt.legend(handles, labels, fontsize=15, loc='upper right')
    plt.savefig('%s_C%s_lim_mag_plot_%s%s' % (survey, chip, num, end_name), dpi=300)
    print('Saved lim mag plot to dir!')
    plt.clf()

    # Crossmatch location check plot
    # if not os.path.isfile('%s_C%s_source_check_plot_%s.png' % (survey, chip, num)):
    # mean, median, sigma_plot = sigma_clipped_stats(data)
    #
    # fig = plt.figure(figsize=(10, 10))
    # ax = fig.gca()
    #
    # im = ax.imshow(
    #     data,
    #     vmin=median - 1.5 * sigma_plot,
    #     vmax=median + 1.5 * sigma_plot,
    #     origin='lower'
    # )
    #
    # # Draw circles
    # circles = [
    #     plt.Circle(
    #         (cleanPSFsources['X_IMAGE'][idx_image][i],
    #          cleanPSFsources['Y_IMAGE'][idx_image][i]),
    #         radius=5,
    #         edgecolor='r',
    #         facecolor='None'
    #     ) for i in range(len(cleanPSFsources['X_IMAGE'][idx_image]))
    # ]
    # for c in circles:
    #     ax.add_artist(c)
    #
    # cbar = fig.colorbar(im, ax=ax)
    # cbar.set_label("Pixel Value")
    #
    # plt.savefig('%s_C%s_source_check_plot_%s%s' % (survey, chip, num, end_name), dpi=150)
    # print('Saved source location check plot to dir!')
    # plt.clf()

    plt.close('all')

    # source location regions
    newtext = open('PRIME_%s_C%s_crsmtched_srcs_%s.reg' % (survey, chip, num), 'w+')
    newtext.write('fk5')
    for a, d, rad in zip(cleanPSFsources['ALPHA_J2000'][idx_image], cleanPSFsources['DELTA_J2000'][idx_image],
                         cleanPSFsources['FLUX_RADIUS'][idx_image]):
        newtext.write(f'\ncircle({a}, {d}, {rad}") # color=red')

    newtext_all = open('PRIME_%s_C%s_all_srcs_%s.reg' % (survey, chip, num), 'w+')
    newtext_all.write('fk5')
    for a, d, rad in zip(PSFsources['ALPHA_J2000'], PSFsources['DELTA_J2000'],
                         PSFsources['FLUX_RADIUS']):
        newtext_all.write(f'\ncircle({a}, {d}, {rad}") # color=green')

    print('Source location reg files saved!')

    print('Writing relevant plot info to image header...')
    with open_fits_robust(imageName) as hdul:
        hdr = hdul[0].header
        hdr.set(f'{MAGTYPES[magtype]}_M', m, 'WLS %s sig fit slope' % sigma, after='SURVEY')
        hdr.set(f'E_{MAGTYPES[magtype]}_M', m_err, 'Error in WLS %s sig fit slope' % sigma, after=f'{MAGTYPES[magtype]}_M')
        hdr.set(f'{MAGTYPES[magtype]}_B', b, 'WLS %s sig fit intercept' % sigma, after=f'E_{MAGTYPES[magtype]}_M')
        hdr.set(f'E_{MAGTYPES[magtype]}_B', b_err, 'Error in WLS %s sig fit intercept' % sigma, after=f'{MAGTYPES[magtype]}_B')
        hdr.set(f'LM_{MAGTYPES[magtype]}', limmag, 'Source Histogram FWHM Limiting Mag', after=f'E_{MAGTYPES[magtype]}_B')

    return m, b, round(3 * b_err, 4)


def photometry_plots(cleanPSFsources, PSFsources, data, imageName, survey, band, good_cat_stars, idx_psfmass, idx_psfimage,
                     psfweights_noclip, psf_clipped, sigma, aperweights_noclip, aper_clipped_all, autoweights_noclip, auto_clipped):

    # TODO temporary disabling of in progress aper mag plotting
    aperweights_noclip = []
    aper_clipped_all = []

    # appropriate mag column
    colnames = good_cat_stars.colnames
    if len(colnames) > 6:
        magcol = f'{band}MAG_{magtype}'
        magerrcol = f'{band}MAG_{magtype}'
    else:
        magcol = colnames[2]
        magerrcol = colnames[3]

    # aperture sizes
    if len(aperweights_noclip) > 0:
        apers_str = fits.getheader(imageName)['APERS']
        aper_arr = apers_str.split(',')

    chip = imageName[-6]
    if len(imageName) <= 16:
        num = 'img'
    else:
        num = imageName[-16:-8]

    # PSF-specific naming for parallel running
    end_name = end_name_gen('png')

    def predict_y_for(x, m, b):
        return m * x + b

    # PHOTOMETRIC FIT LINE MODELS

    model2, model_sig, model_sig_resid, model_auto, model_auto_resid, aper_model_sigs = (
        photometric_fit_calc(cleanPSFsources, band, good_cat_stars, idx_psfmass, idx_psfimage,
                         psfweights_noclip, psf_clipped, sigma, aperweights_noclip, aper_clipped_all,
                         autoweights_noclip, auto_clipped, with_plots=True))

    if model_auto:
        m_auto = model_auto.params[1]
        m_autoerr = model_auto.bse[1]
        b_auto = model_auto.params[0]
        b_autoerr = model_auto.bse[0]
        rsquare_auto = model_auto.rsquared
        rss_auto = model_auto.ssr
        model_auto_resid = model_auto.resid

        m2 = model2.params[1]
        m2err = model2.bse[1]
        b2 = model2.params[0]
        b2err = model2.bse[0]
        rsquare2 = model2.rsquared
        rss2 = model2.ssr
        avg2 = np.average(model2.resid, weights=autoweights_noclip)
        var2 = np.average((model2.resid - avg2) ** 2, weights=autoweights_noclip)

    if model_sig:
        m_sig = model_sig.params[1]
        m_sigerr = model_sig.bse[1]
        b_sig = model_sig.params[0]
        b_sigerr = model_sig.bse[0]
        rsquare_sig = model_sig.rsquared
        rss_sig = model_sig.ssr
        model_sig_resid = model_sig.resid

        m2 = model2.params[1]
        m2err = model2.bse[1]
        b2 = model2.params[0]
        b2err = model2.bse[0]
        rsquare2 = model2.rsquared
        rss2 = model2.ssr
        avg2 = np.average(model2.resid, weights=psfweights_noclip)
        var2 = np.average((model2.resid - avg2) ** 2, weights=psfweights_noclip)

    # BINNED RESIDUAL STATISTICS

    # auto aperture residual fit 3 sig clip - bin errors and stats
    mags = range(10, 22)
    x_arr = np.arange(10 + 0.5, 22 + 0.5, 1)

    auto_res_errs = []
    auto_res_means = []
    auto_res_ranges = []
    auto_res_skews = []
    auto_res_nums = []
    for i in mags:
        auto_mask = ((cleanPSFsources['%sMAG_AUTO' % band][idx_psfimage][~auto_clipped.mask] > i) &
                (cleanPSFsources['%sMAG_AUTO' % band][idx_psfimage][~auto_clipped.mask] < i+1))
        auto_mask = model_auto_resid[auto_mask]
        resauto_err = np.std(auto_mask)   # spread
        resauto_range = '%s - %s' % (i, i + 1)
        resauto_mean = np.nanmean(auto_mask)  # mean
        resauto_skew = skew(auto_mask, bias=False)    # skew
        resauto_num = len(auto_mask)
        if ma.is_masked(resauto_err):
            resauto_err = 0
        auto_res_errs.append(resauto_err)
        auto_res_means.append(resauto_mean)
        auto_res_ranges.append(resauto_range)
        auto_res_skews.append(resauto_skew)
        auto_res_nums.append(resauto_num)
    auto_res_errs = np.nan_to_num(np.array(auto_res_errs))
    auto_res_means = np.nan_to_num(np.array(auto_res_means))
    auto_res_skews = np.nan_to_num(np.array(auto_res_skews))
    auto_res_errs_min = np.min((auto_res_errs[auto_res_errs != 0]))
    auto_res_errs_max = np.max(auto_res_errs)

    # residual fit 3 sig clip - bin errors and stats
    if model_sig:
        res_errs = []
        res_means = []
        res_ranges = []
        res_skews = []
        res_nums = []
        for i in mags:
            mask = ((cleanPSFsources['%sMAG_PSF' % band][idx_psfimage][~psf_clipped.mask] > i) &
                    (cleanPSFsources['%sMAG_PSF' % band][idx_psfimage][~psf_clipped.mask] < i+1))
            res_mask = model_sig_resid[mask]
            # avgsig = np.average(res_mask, weights=psfweights_noclip[~psf_clipped.mask][mask])
            # varsig = np.average((res_mask - avgsig)**2, weights=psfweights_noclip[~psf_clipped.mask][mask])
            # ressig_err = np.sqrt(varsig)    # spread
            ressig_err = np.std(res_mask)   # spread
            ressig_range = '%s - %s' % (i, i + 1)
            ressig_mean = np.nanmean(res_mask)  # mean
            ressig_skew = skew(res_mask, bias=False)    # skew
            ressig_num = len(res_mask)
            if ma.is_masked(ressig_err):
                ressig_err = 0
            res_errs.append(ressig_err)
            res_means.append(ressig_mean)
            res_ranges.append(ressig_range)
            res_skews.append(ressig_skew)
            res_nums.append(ressig_num)
        res_errs = np.nan_to_num(np.array(res_errs))
        res_means = np.nan_to_num(np.array(res_means))
        res_skews = np.nan_to_num(np.array(res_skews))
        res_errs_min = np.min((res_errs[res_errs != 0]))
        res_errs_max = np.max(res_errs)

        # bin errors / stats table
        bintable = Table()
        bintable['Bin Range (mag)'] = [res_ranges, auto_res_ranges]
        bintable['Mean (mag)'] = [res_means, auto_res_means]
        bintable['Spread (mag)'] = [res_errs, auto_res_errs]
        bintable['Skew'] = [res_skews, auto_res_skews]
        bintable['Source Number'] = [res_nums, auto_res_nums]

        bintable.meta['Description'] = ('Table of data on residuals, binned by mag. Each column is a vector w/ the 1st idx '
                                        'being the PSF fit data and the 2nd idx being auto aperture photom fit data')

        bintable.write('Resid_%s-sig_Data_%s_C%s_%s%s' % (sigma, band, chip, survey, end_name_gen()), overwrite=True)
        print('Residual data table for PSF and auto photometry written!')

    else:
        bintable = Table()
        bintable['Bin Range (mag)'] = auto_res_ranges
        bintable['Mean (mag)'] = auto_res_means
        bintable['Spread (mag)'] = auto_res_errs
        bintable['Skew'] = auto_res_skews
        bintable['Source Number'] = auto_res_nums

        bintable.meta['Description'] = 'Table of data on residuals for PSF fit data, binned by mag'
        bintable.write('Resid_%s-sig_Data_%s_C%s_%s.ecsv' % (sigma, band, chip, survey), overwrite=True)
        print('Residual data table for %s sigma clip written!' % sigma)

    # aper photom residual fit sig clip - bin errors and stats
    if len(aperweights_noclip) > 0:
        aper_res_errs_all = []
        aper_res_means_all = []
        aper_res_ranges_all = []
        aper_res_skews_all = []
        aper_res_nums_all = []
        aper_res_max_min_all = []

        for idx, (aperclipped, apermodel) in enumerate(zip(aper_clipped_all, aper_model_sigs)):
            ap_res_errs = []
            ap_res_means = []
            ap_res_ranges = []
            ap_res_skews = []
            ap_res_nums = []
            for mag in mags:
                ap_mask = ((cleanPSFsources['%sMAG_APER' % band][idx_psfimage][~aperclipped.mask] > mag) &
                        (cleanPSFsources['%sMAG_PSF' % band][idx_psfimage][~aperclipped.mask] < mag+1))
                ap_res_mask = apermodel.resid[ap_mask]
                ap_ressig_err = np.std(ap_res_mask)   # spread
                ap_ressig_range = '%s - %s' % (mag, mag + 1)
                ap_ressig_mean = np.nanmean(ap_res_mask)  # mean
                ap_ressig_skew = skew(ap_res_mask, bias=False)    # skew
                ap_ressig_num = len(ap_res_mask)
                if ma.is_masked(ap_ressig_err):
                    ap_ressig_err = 0
                ap_res_errs.append(ap_ressig_err)
                ap_res_means.append(ap_ressig_mean)
                ap_res_ranges.append(ap_ressig_range)
                ap_res_skews.append(ap_ressig_skew)
                ap_res_nums.append(ap_ressig_num)
            ap_res_errs = np.nan_to_num(np.array(ap_res_errs))
            ap_res_means = np.nan_to_num(np.array(ap_res_means))
            ap_res_skews = np.nan_to_num(np.array(ap_res_skews))
            ap_res_errs_min = np.min((ap_res_errs[ap_res_errs != 0]))
            ap_res_errs_max = np.max(ap_res_errs)
            ap_res_max_min = (ap_res_errs_max, ap_res_errs_min)

            aper_res_errs_all.append(ap_res_errs)
            aper_res_means_all.append(ap_res_means)
            aper_res_ranges_all.append(ap_res_ranges)
            aper_res_skews_all.append(ap_res_skews)
            aper_res_nums_all.append(ap_res_nums)
            aper_res_max_min_all.append(ap_res_max_min)

    # ALL PLOTS

    plt.close('all')
    # mag comparison plot
    plt.figure(1, figsize=(8, 8))
    plt.plot(cleanPSFsources[f'{band}MAG_AUTO'][idx_psfimage], good_cat_stars['%s' % magcol][idx_psfmass],
             'r.', markersize=14, markeredgecolor='black')
    plt.xlim(10, 22)
    plt.ylim(10, 22)
    plt.title('PRIME Mags vs %s Mags' % survey)
    plt.xlabel('PRIME %s Mags' % band, fontsize=15)
    plt.ylabel('%s %s Mags' % (survey, band), fontsize=15)
    plt.grid()
    # plt.savefig('%s_C%s_mag_comp_plot_%s.png' % (survey, chip, num))
    plt.clf()
    # print('Saved mag comparison plot to dir!')

    # PRIME flux vs catalog AB mag for crossmatches
    x_flx_lin = cleanPSFsources[f'AUTO_FLUX_DENSITY'][idx_psfimage][~auto_clipped.mask]
    x_flx = np.log(x_flx_lin)
    y_flx = good_cat_stars['%s' % magcol][idx_psfmass][~auto_clipped.mask]
    x_const_flx = sm.add_constant(x_flx)
    model_flx = sm.WLS(y_flx, x_const_flx, weights=autoweights_noclip[~auto_clipped.mask]).fit()
    # print(model_flx.params)
    m_flx = model_flx.params[1]
    m_flxerr = model_flx.bse[1]
    b_flx = model_flx.params[0]
    b_flxerr = model_flx.bse[0]

    x_grid = np.linspace(x_flx_lin.min(), x_flx_lin.max(), 100)
    log_x_grid = np.log(x_grid)
    log_x_grid_const = sm.add_constant(log_x_grid)
    y_pred = model_flx.predict(log_x_grid_const)

    plt.figure(9, figsize=(8, 8))
    plt.plot(cleanPSFsources['AUTO_FLUX_DENSITY'][idx_psfimage][~auto_clipped.mask],
             good_cat_stars['%s' % magcol][idx_psfmass][~auto_clipped.mask],
             'r.', markersize=14, markeredgecolor='black')
    plt.plot(x_grid, y_pred, c='b')

    flx_txt = ('slope = %.4f' % m_flx + '\nslope err = %.4f' % m_flxerr +
           '\nint = %.4f' % b_flx + '\nint err = %.4f' % b_flxerr)

    plt.xlim(10, 50000)
    plt.ylim(12, 20.5)
    plt.title('PRIME Flux Density vs %s AB mag - %s Sigma Clip' % (survey, sigma))
    plt.xlabel(r'PRIME Flux Density ($\mu$Jy)', fontsize=15)
    plt.ylabel('%s %s AB Mags' % (survey, band), fontsize=15)
    plt.grid()
    plt.xscale('log')
    flx_box = dict(facecolor='white')
    plt.text(1000, 19, flx_txt, fontsize=12, bbox=flx_box)
    plt.savefig('%s_C%s_flux_mag_plot_sig_%s%s' % (survey, chip, num, end_name))
    plt.clf()
    print('Saved flux v. mag plot to dir!')


    # residual plot - y int forced to zero
    """
    plt.figure(2, figsize=(8, 6))
    plt.scatter(cleanPSFsources['%sMAG_PSF' % band][idx_psfimage], model.resid, color='red')
    plt.ylim(-1.5, 1.5)
    plt.xlim(10,21)
    plt.title('PRIME vs %s Residuals' % survey)
    plt.ylabel('Residuals')
    plt.xlabel('Mags')
    plt.axhline(y=res_err, color='tab:orange', linestyle='--', linewidth=1)
    plt.axhline(y=res_err * 2, color='green', linestyle='--', linewidth=1)
    plt.axhline(y=-res_err, color='tab:orange', linestyle='--', linewidth=1)
    plt.axhline(y=-res_err * 2, color='green', linestyle='--', linewidth=1)
    plt.axhline(y=0, color='black', linestyle='--', linewidth=1)
    plt.legend(['Residuals', r'1 $\sigma$ = %.3f' % res_err, r'2 $\sigma$ = %.3f' % (2*res_err)], loc='lower left')
    info = ('eqn: y = mx'+'\nslope = %.5f +/- %.5f' % (m,merr))+('\nR$^{2}$ = %.3f' % rsquare)+('\nRSS = %d' % rss)
    #+('\nintercept = %.3f +/- %.3f' % (b,berr))
    plt.text(15,-1.25,info,bbox=dict(facecolor='white',edgecolor='black',alpha=1,pad=5.0))
    plt.savefig('%s_C%s_residual_plot_%s.png' % (survey,chip,num),dpi=300)
    print('Saved residual plot to dir!')
    """

    # res plot, y int include, PSF and aperture fits

    if len(aperweights_noclip) > 0:
        plt.figure(2, figsize=(8, 6))
        psf_sc = plt.scatter(good_cat_stars['%s' % magcol][idx_psfmass][~psf_clipped.mask], model_sig_resid, color='red', alpha=0.5,
                    label='PSF Residuals')
        plt.ylim(-1, 1)
        plt.xlim(10, 21)
        plt.title('%s Residuals - %s Sigma Clip' % (survey, sigma))
        plt.ylabel('Residuals')
        plt.xlabel('%s %s Mags' % (survey, band))
        # plt.axhline(y=ressig_err, color='blue', linestyle='--', linewidth=1)
        # plt.axhline(y=ressig_err * 2, color='green', linestyle='--', linewidth=1)
        # plt.axhline(y=-ressig_err, color='blue', linestyle='--', linewidth=1)
        # plt.axhline(y=-ressig_err * 2, color='green', linestyle='--', linewidth=1)
        plt.scatter(x_arr, res_errs, marker='_', s=1625, c='black')
        plt.scatter(x_arr, -res_errs, marker='_', s=1625, c='black')
        plt.axhline(y=0, color='black', linestyle='--', linewidth=1)
        info2 = ('PSF eqn: y = mx+b' + '\nslope = %.4f +/- %.4f' % (m_sig, m_sigerr)) + (
                    '\nintercept = %.3f +/- %.3f' % (b_sig, b_sigerr)) + ('\nR$^{2}$ = %.3f' % rsquare_sig) + ('\nRSS = %d' % rss_sig)
        plt.text(15, -0.9, info2, fontsize=9, bbox=dict(facecolor='white', edgecolor='black', pad=5.0))

        scatters = []
        color_arr = ['blue', 'green', 'magenta', 'yellow', 'tab:orange']
        for aperclipped, apermodel, aperminmaxes, apererrs, apersize, color in (
                zip(aper_clipped_all, aper_model_sigs, aper_res_max_min_all, aper_res_errs_all, aper_arr, color_arr)):
            sc = plt.scatter(good_cat_stars['%s' % magcol][idx_psfmass][~aperclipped.mask], apermodel.resid, color=color,
                        alpha=0.4, label=f'{apersize * 0.498}" Aperture Residuals')
            scatters.append(sc)
            plt.scatter(x_arr, apererrs, marker='_', s=1625, c='black', alpha=0.4)
            plt.scatter(x_arr, -apererrs, marker='_', s=1625, c='black', alpha=0.4)

        handles = scatters + [psf_sc]
        labels = [h.get_label() for h in handles]
        plt.legend(handles, labels, loc='lower left',
                   markerscale=0.5)

        plt.savefig('%s_C%s_residual_plot_all_%s%s' % (survey, chip, num, end_name), dpi=300)
        plt.clf()

    # res plot y int, histogram
    if len(idx_psfimage) >= 75000:
        bin_num_int = round(len(idx_psfimage) / 500)
    elif 50000 <= len(idx_psfimage) <= 75000:
        bin_num_int = round(len(idx_psfimage) / 400)
    elif 25000 <= len(idx_psfimage) <= 50000:
        bin_num_int = round(len(idx_psfimage) / 300)
    elif 5000 <= len(idx_psfimage) <= 25000:
        bin_num_int = round(len(idx_psfimage) / 75)
    elif 1000 <= len(idx_psfimage) <= 5000:
        bin_num_int = round(len(idx_psfimage) / 50)
    elif len(idx_psfimage) <= 1000:
        bin_num_int = 50

    alpha = 0.7 if model_auto else None

    # PSF PANEL
    if model_sig:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 8), sharey=True)
        fig.subplots_adjust(wspace=0.12)
        fig.suptitle(f'{survey} Residuals - {sigma} Sigma Clip - Density Histogram')

        psfhist = ax1.hist2d(x=good_cat_stars['%s' % magcol][idx_psfmass][~psf_clipped.mask],y=model_sig_resid,
            bins=[bin_num_int, bin_num_int], range=[[10, 21], [-1, 1]],
            cmap='gist_heat_r'
        )
        cbar1 = fig.colorbar(psfhist[3], ax=ax1, pad=0.03)
        cbar1.set_label('PSF Density')

        psfsig = ax1.scatter(x_arr, res_errs, marker='_', s=1625, c='blue',
                             label=r'PSF photom 1 $\sigma$ range = [%.3f - %.3f]' % (res_errs_min, res_errs_max))
        ax1.scatter(x_arr, -res_errs, marker='_', s=1625, c='blue')

        ax1.axhline(y=0, color='black', linestyle='--', linewidth=1)
        ax1.set_xlim(10, 21)
        ax1.set_ylim(-1, 1)
        ax1.set_title(f"PSF Fit Residuals")
        ax1.set_xlabel(f"{survey} {band} Mags")
        ax1.set_ylabel("Residuals")

        infohist = (
                'PSF fit'
                + '\nslope = %.4f +/- %.4f' % (m_sig, m_sigerr)
                + '\nintercept = %.3f +/- %.3f' % (b_sig, b_sigerr)
                + '\nR$^{2}$ = %.3f' % rsquare_sig
                + '\nRSS = %d' % rss_sig
                + '\nn_sources = %i' % len(cleanPSFsources['PSF_FLUX_DENSITY'][idx_psfimage][~psf_clipped.mask])
        )

        ax1.text(10.5, 0.5, infohist, fontsize=9,
                 bbox=dict(facecolor='white', edgecolor='black', pad=5.0))

        ax1.legend([psfsig], [psfsig.get_label()], loc='lower left', markerscale=0.5)

    else:
        fig, ax2 = plt.subplots(1, 1, figsize=(9, 8))
        fig.suptitle(f'{survey} Residuals - {sigma} Sigma Clip - Density Histogram')

    # AUTO PANEL
    autohist = ax2.hist2d(
        x=good_cat_stars['%s' % magcol][idx_psfmass][~auto_clipped.mask],y=model_auto_resid,
        bins=[bin_num_int, bin_num_int], range=[[10, 21], [-1, 1]],
        cmap='gist_earth_r'
    )

    cbar2 = fig.colorbar(autohist[3], ax=ax2, pad=0.03)
    cbar2.set_label('Auto Ap. Density')

    autosig = ax2.scatter(x_arr, auto_res_errs, marker='_', s=1625, c='black',
                          label=r'Auto ap. photom 1 $\sigma$ range = [%.3f - %.3f]' % (
                          auto_res_errs_min, auto_res_errs_max))
    ax2.scatter(x_arr, -auto_res_errs, marker='_', s=1625, c='black')

    ax2.axhline(y=0, color='black', linestyle='--', linewidth=1)
    ax2.set_xlim(10, 21)
    ax2.set_ylim(-1, 1)
    ax2.yaxis.set_tick_params(labelleft=True)
    ax2.set_title(f"Auto Aperture Fit Residuals")
    ax2.set_xlabel(f"{survey} {band} Mags")
    ax2.set_ylabel("Residuals")

    infoauto = (
            'Auto Ap. fit'
            + '\nslope = %.4f +/- %.4f' % (m_auto, m_autoerr)
            + '\nintercept = %.3f +/- %.3f' % (b_auto, b_autoerr)
            + '\nR$^{2}$ = %.3f' % rsquare_auto
            + '\nRSS = %d' % rss_auto
            + '\nn_sources = %i' % len(cleanPSFsources['AUTO_FLUX_DENSITY'][idx_psfimage][~auto_clipped.mask])
    )

    ax2.text(10.5, 0.5, infoauto, fontsize=9,
             bbox=dict(facecolor='white', edgecolor='black', pad=5.0))

    ax2.legend([autosig], [autosig.get_label()], loc='lower left', markerscale=0.5)

    plt.savefig('%s_C%s_residual_plot_int_hist_%s%s' % (survey, chip, num, end_name), dpi=300)
    plt.close(fig)

    print('Saved y-int residual plots to dir!')

    # WLS fit line over data plot

    txt = ('slope = %.4f' % m2 + '\nslope err = %.4f' % m2err + '\nint = %.4f' % b2 + '\nint err = %.4f' % b2err +
           '\nn_sources = %i' % len(cleanPSFsources[idx_psfimage]))
    #
    # plt.figure(4, figsize=(8, 8))
    # plt.xlim(10, 22)
    # plt.ylim(10, 22)
    # plt.title('%s vs PRIME w/ Weighted Fit' % survey)
    # plt.grid()
    # plt.ylabel('PRIME %s Mags' % band, fontsize=15)
    # plt.xlabel('%s %s Mags' % (survey, band), fontsize=15)
    # plt.scatter(good_cat_stars['%s' % magcol][idx_psfmass], cleanPSFsources['%sMAG_AUTO' % band][idx_psfimage])
    # plt.plot(good_cat_stars['%s' % magcol][idx_psfmass],
    #          predict_y_for(good_cat_stars['%s' % magcol][idx_psfmass], m2, b2), c='r')
    # box = dict(facecolor='white')
    # plt.text(11, 18, txt, fontsize=12, bbox=box)
    # plt.savefig('%s_C%s_WLS_fit_plot_%s.png' % (survey, chip, num), dpi=300)

    # WLS hist density plot
    if len(idx_psfimage) >= 75000:
        bin_num = round(len(idx_psfimage) / 500)
    elif 50000 <= len(idx_psfimage) <= 75000:
        bin_num = round(len(idx_psfimage) / 350)
    elif 25000 <= len(idx_psfimage) <= 50000:
        bin_num = round(len(idx_psfimage) / 150)
    elif 5000 <= len(idx_psfimage) <= 25000:
        bin_num = round(len(idx_psfimage) / 50)
    elif 1000 <= len(idx_psfimage) <= 5000:
        bin_num = round(len(idx_psfimage) / 20)
    elif len(idx_psfimage) <= 1000:
        bin_num = 100

    plt.figure(5, figsize=(10, 8))
    plt.clf()
    plt.xlim(10, 22)
    plt.ylim(10, 22)
    plt.title('%s vs PRIME w/ Weighted Fit - Density Histogram' % survey)
    plt.grid()
    plt.ylabel('PRIME %s Mags' % band, fontsize=15)
    plt.xlabel('%s %s Mags' % (survey, band), fontsize=15)
    plt.hist2d(x=good_cat_stars['%s' % magcol][idx_psfmass], y=cleanPSFsources['%sMAG_AUTO' % band][idx_psfimage],
               bins=[bin_num, bin_num], range=[[10, 22],[10, 22]], cmap='gist_heat_r')
    plt.plot(good_cat_stars['%s' % magcol][idx_psfmass],
             predict_y_for(good_cat_stars['%s' % magcol][idx_psfmass], m2, b2), c='b')
    plt.colorbar(label='Density')
    box = dict(facecolor='white')
    plt.text(11, 18, txt, fontsize=12, bbox=box)
    plt.savefig('%s_C%s_WLS_fit_hist_plot_%s%s' % (survey, chip, num, end_name), dpi=300)
    plt.clf()

    print('Saved WLS fit plots to dir!')

    # WLS fit for 3 sig clip of data
    # sigtxt = ('slope = %.4f' % m_sig + '\nslope err = %.4f' % m_sigerr + '\nint = %.4f' % b_sig +
    #           '\nint err = %.4f' % b_sigerr + '\nn_sources = %i' % len(cleanPSFsources[idx_psfimage][~psf_clipped.mask]))
    #
    # plt.figure(6, figsize=(8, 8))
    # plt.clf()
    # plt.xlim(10, 22)
    # plt.ylim(10, 22)
    # plt.title('%s vs PRIME w/ Weighted Fit - %s Sigma Clip' % (survey, sigma))
    # plt.grid()
    # plt.ylabel('PRIME %s Mags' % band, fontsize=15)
    # plt.xlabel('%s %s Mags' % (survey, band), fontsize=15)
    # plt.scatter(good_cat_stars['%s' % magcol][idx_psfmass][~psf_clipped.mask],
    #             cleanPSFsources['%sMAG_AUTO' % band][idx_psfimage][~psf_clipped.mask])
    # plt.plot(good_cat_stars['%s' % magcol][idx_psfmass][~psf_clipped.mask],
    #          predict_y_for(good_cat_stars['%s' % magcol][idx_psfmass][~psf_clipped.mask], m_sig, b_sig), c='r')
    # box = dict(facecolor='white')
    # plt.text(11, 18, sigtxt, fontsize=12, bbox=box)
    # plt.savefig('%s_C%s_WLS_fit_3sig_plot_%s.png' % (survey, chip, num), dpi=300)

    # WLS 3 sig hist density plot

    # PSF PANEL
    if model_sig:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 8), sharey=True)
        fig.subplots_adjust(wspace=0.12)
        fig.suptitle(f'{survey} vs PRIME w/ Weighted Fit - {sigma} Sigma Clip - Density Histogram')

        psfhist = ax1.hist2d(x=good_cat_stars['%s' % magcol][idx_psfmass][~psf_clipped.mask],
                   y=cleanPSFsources['%sMAG_PSF' % band][idx_psfimage][~psf_clipped.mask],
                   bins=[bin_num, bin_num], range=[[10, 22],[10, 22]], cmap='gist_heat_r')

        psfline = ax1.plot(good_cat_stars['%s' % magcol][idx_psfmass][~psf_clipped.mask],
                 predict_y_for(good_cat_stars['%s' % magcol][idx_psfmass][~psf_clipped.mask], m_sig, b_sig), c='b')

        cbar1 = fig.colorbar(psfhist[3], ax=ax1, pad=0.03)
        cbar1.set_label('PSF Density')

        ax1.grid()
        ax1.set_xlim(10, 22)
        ax1.set_ylim(10, 22)
        ax1.set_title(f"{survey} vs PRIME PSF Mag Plot w/ WLS fit line")
        ax1.set_xlabel(f"{survey} {band} Mags")
        ax1.set_ylabel(f"PRIME {band} Mags")

        ax1.text(11, 18, infohist, fontsize=12,
                 bbox=dict(facecolor='white', edgecolor='black'))

    else:
        fig, ax2 = plt.subplots(1, 1, figsize=(9, 8))
        fig.suptitle(f'{survey} vs PRIME w/ Weighted Fit - {sigma} Sigma Clip - Density Histogram')

    # AUTO PANEL
    autohist = ax2.hist2d(
        x=good_cat_stars['%s' % magcol][idx_psfmass][~auto_clipped.mask],
        y=cleanPSFsources['%sMAG_AUTO' % band][idx_psfimage][~auto_clipped.mask],
        bins=[bin_num, bin_num], range=[[10, 22],[10, 22]],
        cmap='gist_earth_r'
    )

    autoline = ax2.plot(good_cat_stars['%s' % magcol][idx_psfmass][~auto_clipped.mask],
                       predict_y_for(good_cat_stars['%s' % magcol][idx_psfmass][~auto_clipped.mask], m_auto, b_auto),
                       c='r')

    cbar2 = fig.colorbar(autohist[3], ax=ax2, pad=0.03)
    cbar2.set_label('Auto Ap. Density')

    ax2.grid()
    ax2.set_xlim(10, 22)
    ax2.set_ylim(10, 22)
    ax2.yaxis.set_tick_params(labelleft=True)
    ax2.set_title(f"{survey} vs PRIME Auto Ap. Plot w/ WLS fit line")
    ax2.set_xlabel(f"{survey} {band} Mags")
    ax2.set_ylabel(f"PRIME {band} Mags")

    ax2.text(11, 18, infoauto, fontsize=12,
             bbox=dict(facecolor='white', edgecolor='black'))

    plt.savefig('%s_C%s_WLS_fit_3sig_hist_plot_%s%s' % (survey, chip, num, end_name), dpi=300)
    plt.close(fig)

    print('Saved WLS 3 sig fit plots to dir!')

    # flux vs mag plot - histogram vers.
    plt.figure(10, figsize=(8, 8))
    plt.hist2d(x=cleanPSFsources['AUTO_FLUX_DENSITY'][idx_psfimage], y=good_cat_stars['%s' % magcol][idx_psfmass],
              bins=[bin_num, bin_num], range=[[100, 50000],[12, 20.5]], cmap='gist_heat_r')
    plt.colorbar(label='Density')
    plt.xlim(10, 50000)
    plt.ylim(12, 20.5)
    plt.title('PRIME Flux Density vs %s AB mag - Histogram' % survey)
    plt.xlabel(r'PRIME Flux Density ($\mu$Jy)', fontsize=15)
    plt.ylabel('%s %s AB Mags' % (survey, band), fontsize=15)
    plt.grid()
    plt.xscale('log')
    # plt.savefig('%s_C%s_flux_mag_hist_plot_%s.png' % (survey, chip, num))
    plt.clf()
    # print('Saved flux v. mag hist plot to dir!')

    # Limiting Mag Plot

    def lim_mag_calc(bin_vals, all_mags):
        idxs = np.digitize(all_mags, bins=bin_vals)

        indices = {i: [] for i in range(len(bin_vals))}
        for idx, value in enumerate(idxs):
            indices[value].append(idx)
        sorted_indices_lists = list(indices.values())

        all_sources = []
        for i in sorted_indices_lists:
            number = len(i)
            all_sources.append(number)

        # plotting
        idx_arr = np.where(np.isclose(bin_vals, 12.5))
        min_x_idx = idx_arr[0][0]   # avoid saturated <12.5 mag sources from affecting maximum

        filtered_sources = all_sources[min_x_idx:]
        idxmax = filtered_sources.index(max(filtered_sources)) + min_x_idx
        split_sources = all_sources[idxmax:]

        halfmax = max(filtered_sources) / 2
        halfmaxpt = list(max(enumerate(split_sources), key=lambda x: -abs(halfmax - x[1])))
        halfmaxpt = [halfmaxpt[0] + idxmax, halfmaxpt[1]]

        limmag = round(bin_vals[halfmaxpt[0]], 1)

        return all_sources, halfmax, limmag

    all_mags_all = PSFsources[PSFsources['%sMAG_AUTO' % band] < 25]
    all_mags = all_mags_all['%sMAG_AUTO' % band]

    bin_vals = np.array(np.arange(12, 25.5, 0.1))

    all_sources, halfmax, limmag = lim_mag_calc(bin_vals=bin_vals, all_mags=all_mags)

    print('Auto Ap. Lim Mag = ', limmag)

    plt.figure(8, figsize=(24, 8))

    autodata = plt.bar(bin_vals, height=all_sources, width=0.1, align='edge', color='blue', edgecolor='black',
                       alpha=1, label='Auto Ap. Binned Sources')
    autohalf = plt.axhline(halfmax, linestyle='--', color='black',
                           label='Auto ap. Half Max = %s' % round(halfmax, 1))
    autolimmag = plt.axvline(limmag, color='black', linewidth=2,
                             label='Auto ap. Limiting Mag = %s' % round(limmag, 1))

    psfhalf = psf_limmag = []
    if model_sig:
        psf_mags = PSFsources[PSFsources['%sMAG_PSF' % band] < 25]
        all_psf_mags = psf_mags['%sMAG_PSF' % band]

        all_psf_sources, psf_halfmax, psf_limmag = lim_mag_calc(bin_vals=bin_vals, all_mags=all_psf_mags)

        psfdata = plt.bar(bin_vals, height=all_psf_sources, width=0.1, align='edge', color='red', edgecolor='black',
                          alpha=alpha,
                          label='PSF Binned Sources')
        psfhalf = plt.axhline(psf_halfmax, linestyle='--', label='PSF Half Max = %s' % round(psf_halfmax, 1))
        psflimmag = plt.axvline(psf_limmag, color='b', linewidth=2, label='PSF Limiting Mag = %s' % round(psf_limmag, 1))
        print('PSF Lim Mag = ', psf_limmag)

    xticks = np.arange(12, 25.5, 0.5)
    plt.xticks(xticks, fontsize=10)
    plt.grid()
    plt.yscale('log')
    plt.title('PRIME Limiting Mag Plot')
    plt.ylabel('Number of Sources')
    plt.xlabel('%s Magnitude' % band)

    if psfhalf:
        handles = [psfdata, psfhalf, psflimmag, autodata, autohalf, autolimmag]
    else:
        handles = [autohalf, autolimmag]
    labels = [h.get_label() for h in handles]
    plt.legend(handles, labels, fontsize=15, loc='upper right')
    plt.savefig('%s_C%s_lim_mag_plot_%s%s' % (survey, chip, num, end_name), dpi=300)
    print('Saved lim mag plot to dir!')
    plt.clf()

    # Crossmatch location check plot
    # if not os.path.isfile('%s_C%s_source_check_plot_%s.png' % (survey, chip, num)):
    mean, median, sigma_plot = sigma_clipped_stats(data)

    fig = plt.figure(figsize=(10, 10))
    ax = fig.gca()

    im = ax.imshow(
        data,
        vmin=median - 1.5 * sigma_plot,
        vmax=median + 1.5 * sigma_plot,
        origin='lower'
    )

    # Draw circles
    circles = [
        plt.Circle(
            (cleanPSFsources['X_IMAGE'][idx_psfimage][i],
             cleanPSFsources['Y_IMAGE'][idx_psfimage][i]),
            radius=5,
            edgecolor='r',
            facecolor='None'
        ) for i in range(len(cleanPSFsources['X_IMAGE'][idx_psfimage]))
    ]
    for c in circles:
        ax.add_artist(c)

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Pixel Value")

    plt.savefig('%s_C%s_source_check_plot_%s%s' % (survey, chip, num, end_name), dpi=150)
    print('Saved source location check plot to dir!')
    plt.clf()

    plt.close('all')

    # source location regions
    newtext = open('%s_C%s_crsmtched_srcs_%s%s' % (survey, chip, num, end_name_gen('reg')), 'w+')
    newtext.write('fk5')
    for a, d, rad in zip(cleanPSFsources['ALPHA_J2000'][idx_psfimage], cleanPSFsources['DELTA_J2000'][idx_psfimage],
                    cleanPSFsources['FLUX_RADIUS'][idx_psfimage]):
        newtext.write(f'\ncircle({a}, {d}, {rad}") # color=red')

    newtext = open('%s_C%s_all_srcs_%s.reg' % (survey, chip, num), 'w+')
    newtext.write('fk5')
    for a, d, rad in zip(PSFsources['ALPHA_J2000'], PSFsources['DELTA_J2000'],
                    PSFsources['FLUX_RADIUS']):
        newtext.write(f'\ncircle({a}, {d}, {rad}") # color=green')

    print('Source location reg files saved!')

    print('Writing relevant plot info to image header...')
    with open_fits_robust(imageName) as hdul:
        hdr = new_hdu_gen_or_set(hdul)
        hdr.set('AUTO_M', m_auto, 'WLS %s sig fit slope' % sigma, after='SURVEY')
        hdr.set('E_AUTO_M', m_autoerr, 'Error in WLS %s sig fit slope' % sigma, after='AUTO_M')
        hdr.set('AUTO_B', b_auto, 'WLS %s sig fit intercept' % sigma, after='E_AUTO_M')
        hdr.set('E_AUTO_B', b_autoerr, 'Error in WLS %s sig fit intercept' % sigma, after='AUTO_B')
        hdr.set('LM_AUTO', limmag, 'Source Histogram FWHM Limiting Mag', after='E_AUTO_B')
        if model_sig:
            try:
                hdr.set('PSF_M', m_sig, 'WLS %s sig fit slope' % sigma, after='auto_fit_m')
            except KeyError:
                hdr.set('PSF_M', m_sig, 'WLS %s sig fit slope' % sigma, after='SURVEY')
            hdr.set('E_PSF_M', m_sigerr, 'Error in WLS %s sig fit slope' % sigma, after='PSF_M')
            hdr.set('PSF_B', b_sig, 'WLS %s sig fit intercept' % sigma, after='E_PSF_M')
            hdr.set('E_PSF_B', b_sigerr, 'Error in WLS %s sig fit intercept' % sigma, after='PSF_B')
            hdr.set('LM_PSF', psf_limmag, 'Source Histogram FWHM Limiting Mag', after='E_PSF_B')

    # return m_sig, b_sig, round(3 * m_sigerr, 4)
    # TODO to run calibration based on auto aperture photom, uncomment line below, comment above
    return m_auto, b_auto, round(3 * b_autoerr, 4)

#%% automatic grb threshold calculation


def grb_rad_convert(rad):
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


#%%
def comb_mag_catalogs(directory, name, chosen_survey):
    """
    Combines separate ecsv cols (w/ diff magtypes) into 1 for ease of data access / grb functions
    """
    all_magtypes = set(MAGTYPES.keys())

    main_ecsvs = []
    ecsvs = [file for file in os.listdir(directory) if name in file and chosen_survey in file and file.endswith('.ecsv')]
    for file in ecsvs:
        if any(f in file for f in all_magtypes):
            main_ecsvs.append(file)
    main_ecsvs = sorted(main_ecsvs)
    print(f'Combining catalogs: {main_ecsvs}')

    if len(main_ecsvs) == 1:
        prime_ecsv_name = ''.join(main_ecsvs)
        new_prime_ecsv_name = '%s.%s.ecsv' % (name, chosen_survey)
        os.rename(prime_ecsv_name, new_prime_ecsv_name)

    elif len(main_ecsvs) > 1:
        prime_ecsv = Table.read(sorted(main_ecsvs)[0])
        comb_check = [col for col in prime_ecsv.colnames if any(magtype in col for magtype in MAGTYPES)]

        psf_check = any('PSF' in ecsv for ecsv in comb_check)
        if psf_check:
            acc_col_num = len(main_ecsvs) * 4
        else:
            acc_col_num = len(main_ecsvs) * 6

        if len(comb_check) <= acc_col_num:
            for cat in main_ecsvs[1:]:
                try:
                    cat_magtype = cat.split('.')[-2]
                    tbl = Table.read(cat)
                    cat_mag_cols = [col for col in tbl.colnames[-4:] if cat_magtype in col]
                    for col_name in cat_mag_cols:
                        prime_ecsv.add_column(tbl[col_name], name=col_name)
                except ValueError:
                    pass
            prime_ecsv.write('%s.%s.ecsv' % (name, chosen_survey), overwrite=True)
            for old_ecsv in main_ecsvs:
                os.remove(old_ecsv)
        else:
            raise Exception('Could not combine .ecsvs!')


#%% automated y int fit calibration


def int_calibration(
        name, directory, band, chip, crop, sigma, Q, chosen_survey, given_catalog, survey,
        mag_low_lim, mag_high_lim, grb_ra, grb_dec,
        grb_coordlist, grb_radius, grb_name, max_int, comp_lvl, sync_signal, sync_queue, continue_queue, make_plots=False, magtype=defaults['magtype'],
        parallel=defaults['parallel']
):
    print('\n# Starting new intercept calibration interation!')
    print('3 sigma fit y-intercept > %s! Redoing photometry w/ sigma = %s, mag low cutoff = %s\n' % (max_int, sigma, mag_low_lim))
    data, header, w, raImage, decImage, bulge, det_thresh, chip = img(directory, name, crop)
    # Q, chosen_survey, mag_low_cutoff = query(raImage, decImage, band, w, data, crop, comp_lvl, given_catalog_path=given_catalog, mag_lower_lim=mag_low_lim,
    #                                          mag_upper_lim=mag_high_lim, bulge=bulge, no_check=True)
    magcol = Q[0].colnames[2]
    key0 = list(Q.keys())[0]
    Q = type(Q)([(key0, Q[0][Q[0][magcol] > mag_low_lim])])

    head_name = os.path.split(name)[0]
    if grb_ra:
        psfcatalogName = [f for f in os.listdir(directory) if f.endswith(f'.photom.cat') and f'C{chip}' in f
                          and '.fits.photom' not in f and head_name in f]
    else:
        psfcatalogName = [f for f in os.listdir(directory) if f.endswith(f'.fits.photom.cat') and f'C{chip}' in f and
                          head_name in f]

    psfcatalogName = ''.join(psfcatalogName)
    good_cat_stars, cleanPSFSources, PSFSources, idx_psfmass, idx_psfimage, massCatCoords, crop_cat_stars = (
        tables(Q, data, w, psfcatalogName, crop, given_catalog))

    cleanPSFSources, PSFSources, weights_noclip, clipped, ab_cat_stars = single_zeropt(
        good_cat_stars, crop_cat_stars, cleanPSFSources, PSFSources, idx_psfmass, idx_psfimage, name, band, survey, sigma,
        data, crop)

    if make_plots:
        slope, intercept, int_err = single_plots(cleanPSFSources, PSFSources, data, name, survey, band, good_cat_stars,
                                                 idx_psfmass, idx_psfimage, sigma, weights_noclip, clipped, crop)
        if sync_queue:
            sync_queue.put("Calib. done, combining .ecsv files")
            continue_queue.get()

        if grb_ra:
            if grb_radius > 60:
                newsourcesearch(grb_ra, grb_dec, name, chosen_survey, band, grb_radius, massCatCoords, ab_cat_stars, directory,
                    chip, grbname=grb_name, mag_low_lim=mag_low_lim)
                # GRB(grb_ra, grb_dec, name, chosen_survey, band, grb_radius, massCatCoords, ab_cat_stars, directory)
            else:
                GRB(grb_ra, grb_dec, name, chosen_survey, band, grb_radius, massCatCoords, ab_cat_stars, directory,
                    chip, grbname=grb_name)
        elif grb_coordlist:
            GRB(grb_ra, grb_dec, name, chosen_survey, band, grb_radius, massCatCoords, ab_cat_stars, directory, chip,
                grb_coordlist, grbname=grb_name)

        return intercept, int_err, True

    else:
        slope, intercept, int_err = single_fit_calc(cleanPSFSources, band, good_cat_stars, idx_psfmass, idx_psfimage,
                                                    weights_noclip, clipped, sigma, with_plots=False, magtype=MAGTYPES[magtype])

        plots_args_dict = dict(
            cleanPSFsources=cleanPSFSources,
            PSFsources=PSFSources,
            data=data,
            imageName=name,
            survey=chosen_survey,
            band=band,
            good_cat_stars=ab_cat_stars,
            idx_mass=idx_psfmass,
            idx_image=idx_psfimage,
            sigma=sigma,
            weights_noclip=weights_noclip,
            clipped=clipped,
            crop=crop
        )

        if grb_ra:
            grb_args_dict = dict(
                ra=grb_ra,
                dec=grb_dec,
                imageName=name,
                survey=chosen_survey,
                band=band,
                thresh=grb_radius,
                massCatCoords=massCatCoords,
                good_cat_stars=ab_cat_stars,
                directory=directory,
                chip=chip,
                coordlist=grb_coordlist,
                grbname=grb_name,
                mag_low_lim=PHOTOMETRY_MAG_LOWER_LIMIT
            )
        else:
            grb_args_dict = {}

        return intercept, int_err, False, plots_args_dict, grb_args_dict


# full intercept calibration loop
def full_int_calibration(
        name, directory, band, chip, crop, data, sigma, Q, chosen_survey, given_catalog, good_cat_stars, crop_cat_stars, cleanPSFSources, PSFSources, massCatCoords,
        idx_psfmass, idx_psfimage, mag_high_lim, mag_low_cutoff, grb_ra, grb_dec, grb_coordlist, grb_name, grb_thresh, comp_lvl, sync_queue=None, continue_queue=None, no_plots=False,
        input_magtype=defaults['magtype'], input_parallel=False
):

    global magtype
    global parallel
    magtype = input_magtype
    parallel = input_parallel

    with (log_output(parallel, f"{name}_int_cal_{MAGTYPES[magtype]}.log")):

        sync_signal = False
        synced = False

        if MAGTYPES[magtype] == 'PSF':
            img_hdr = fits.getheader(name)
            tot_exptime = float(img_hdr['EXPTIME'])
            if tot_exptime < 130:
                multiplier = 1
                sigma = sigma * multiplier
                print(f'Sigma value increased for PSF in lower exp. time field: sigma = {sigma}')

        cleanPSFSources, PSFSources, weights_noclip, clipped, ab_cat_stars = single_zeropt(
                                                    good_cat_stars, crop_cat_stars, cleanPSFSources, PSFSources, idx_psfmass,
                                                    idx_psfimage, name, band, chosen_survey, sigma, data, crop
                                                )

        slope, intercept, int_err = single_fit_calc(cleanPSFSources, band, good_cat_stars, idx_psfmass, idx_psfimage,
                                                    weights_noclip, clipped, sigma, with_plots=False, magtype=MAGTYPES[magtype])
        grb_arg_dict = dict(
            ra=grb_ra,
            dec=grb_dec,
            imageName=name,
            survey=chosen_survey,
            band=band,
            thresh=grb_thresh,
            massCatCoords=massCatCoords,
            good_cat_stars=ab_cat_stars,
            directory=directory,
            chip=chip,
            coordlist=grb_coordlist,
            grbname=grb_name,
            mag_low_lim=PHOTOMETRY_MAG_LOWER_LIMIT
        )

        plots_arg_dict = dict(
            cleanPSFsources=cleanPSFSources,
            PSFsources=PSFSources,
            data=data,
            imageName=name,
            survey=chosen_survey,
            band=band,
            good_cat_stars=ab_cat_stars,
            idx_mass=idx_psfmass,
            idx_image=idx_psfimage,
            sigma=sigma,
            weights_noclip=weights_noclip,
            clipped=clipped,
            crop=crop
        )

        if abs(intercept) >= 6:
            print('Significant photometric intercept value!: %s' % intercept)
            print('Photometric calibration likely unreliable! Is there an issue with the image, catalog, '
                  'or band? Moving on...')
            print('*RECOMMEND DOUBLE-CHECKING THIS FIELD*')
        else:
            prev_intercept = intercept
            revert_flag = False

            while abs(intercept) > int_err:
                print('\nAdjusted minimum acc. intercept: ', int_err)
                print('\nIntercept = %.4f\n' % intercept)
                # sigma -= 0.5
                mag_low_cutoff += 0.5
                new_intercept, new_int_err, synced, plots_arg_dict, grb_arg_dict = int_calibration(name, directory, band, chip,
                                                                                           crop, sigma, Q,
                                                                                           chosen_survey, given_catalog,
                                                                                           chosen_survey,
                                                                                           mag_low_cutoff, mag_high_lim,
                                                                                           grb_ra, grb_dec, grb_coordlist,
                                                                                           grb_thresh, grb_name,
                                                                                           max_int=int_err,
                                                                                           comp_lvl=comp_lvl,
                                                                                           magtype=magtype,
                                                                                           parallel=parallel,
                                                                                           sync_signal=sync_signal,
                                                                                           sync_queue=sync_queue,
                                                                                           continue_queue=continue_queue
                                                                                           )
                if abs(new_intercept) > abs(prev_intercept):
                    print("\nNew intercept: %.4f is higher than previous: %.4f! Reverting and "
                          "redoing...\n" % (new_intercept, prev_intercept))
                    intercept = prev_intercept
                    mag_low_cutoff -= 0.5
                    new_intercept, new_int_err, synced, _, _ = int_calibration(name, directory, band, chip, crop, sigma, Q,
                                                                       chosen_survey, given_catalog,
                                                                       chosen_survey,
                                                                       mag_low_cutoff, mag_high_lim, grb_ra, grb_dec,
                                                                       grb_coordlist,
                                                                       grb_thresh, grb_name, max_int=int_err,
                                                                       comp_lvl=comp_lvl, magtype=magtype, parallel=parallel,
                                                                       sync_signal=sync_signal, sync_queue=sync_queue,
                                                                       continue_queue=continue_queue)
                    revert_flag = True
                    break
                else:
                    intercept = new_intercept
                    prev_intercept = intercept
                    int_err = new_int_err
                    revert_flag = False
            if revert_flag:
                print("Loop stopped due to intercept reverting to the previous value: %.4f" % intercept)
            else:
                print(f"Final intercept below {int_err}: %.4f" % intercept)

            prev_intercept = intercept
            revert_flag = False
            success_flag = False
            leniency_val = 1.15  # 1.1    # val above prev_intercept the new int. can be and still be accepted

            while abs(intercept) > int_err and sigma > 1:  # ensure sigma doesn't go negative
                print('\nAdjusted minimum acc. intercept: ', int_err)
                print('\nIntercept = %.4f\n' % intercept)
                step = 0.125 if sigma <= 1.5 else 0.5
                sigma -= step
                new_intercept, new_int_err, synced, plots_arg_dict, grb_arg_dict = int_calibration(name, directory, band, chip,
                                                                                           crop, sigma, Q,
                                                                                           chosen_survey, given_catalog,
                                                                                           chosen_survey,
                                                                                           mag_low_cutoff, mag_high_lim,
                                                                                           grb_ra, grb_dec, grb_coordlist,
                                                                                           grb_thresh, grb_name,
                                                                                           max_int=int_err,
                                                                                           comp_lvl=comp_lvl,
                                                                                           magtype=magtype,
                                                                                           parallel=parallel,
                                                                                           sync_signal=sync_signal,
                                                                                           sync_queue=sync_queue,
                                                                                           continue_queue=continue_queue
                                                                                           )
                if abs(new_intercept) > abs(leniency_val * prev_intercept):
                    print("\nNew intercept: %.4f is higher than previous: %.4f! Reverting and "
                          "redoing...\n" % (new_intercept, prev_intercept))
                    # revert
                    intercept = prev_intercept
                    sigma += step
                    new_intercept, new_int_err, synced = int_calibration(name, directory, band, chip, crop, sigma, Q,
                                                                 chosen_survey, given_catalog,
                                                                 chosen_survey,
                                                                 mag_low_cutoff, mag_high_lim, grb_ra, grb_dec,
                                                                 grb_coordlist,
                                                                 grb_thresh, grb_name, max_int=int_err,
                                                                 comp_lvl=comp_lvl, make_plots=True, magtype=magtype,
                                                                 parallel=parallel, sync_signal=sync_signal, sync_queue=sync_queue,
                                                                continue_queue=continue_queue)
                    revert_flag = True
                    break
                else:
                    if abs(leniency_val * prev_intercept) > abs(new_intercept) > abs(prev_intercept):
                        print(f'New int. slightly higher than prev.: {round(new_intercept, 3)} > '
                              f'{round(prev_intercept, 3)}, but '
                              f'w/in {leniency_val} leniency val, so accepted')
                    intercept = new_intercept
                    prev_intercept = intercept
                    int_err = new_int_err
                    revert_flag = False

            if revert_flag:
                print("Sigma loop stopped due to intercept reverting to the previous value: %.4f" % intercept)
            else:
                print(f"Final intercept after sigma tuning: %.4f" % intercept)
                success_flag = True
            if success_flag and not no_plots:
                print('\nSuccess! Generating plots and GRB data (if applicable)...\n')

                if not no_plots:
                    single_plots(**plots_arg_dict)

                if not synced:
                    if input_parallel:
                        sync_queue.put("Calib. done, combining .ecsv files")
                        continue_queue.get()
                    else:
                        comb_mag_catalogs(directory, name, chosen_survey)

                if not parallel:
                    if grb_ra:
                        keep = True
                        if grb_dec is None:
                            print('Only GRB RA is found, GRB Dec is None!  Cant conduct grb analysis, '
                                  'make sure the -grb_dec flag is correctly formatted!')
                        if grb_thresh > 60:
                            newsourcesearch(**grb_arg_dict)
                        else:
                            GRB(**grb_arg_dict)
                    elif grb_coordlist:
                        GRB(**grb_arg_dict)


# %% optional removal of intermediate files

def removal(directory):
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

#%%


def photometry(
        full_filename=defaults['filepath'], band=defaults['band'], crop=defaults['crop'], sigma=defaults['sigma_photom'], given_catalog=defaults['catalog'], survey=defaults['survey'],
        mag_low_lim=defaults['mag_low'], mag_high_lim=defaults['mag_high'], no_plots=defaults['no_plots'],
        keep=defaults['keep'], grb_only=defaults['grb_only'], grb_ra=defaults['grb_ra'], grb_dec=defaults['grb_dec'], grb_coordlist=defaults['grb_coordlist'],
        grb_radius=defaults['grb_radius'], grb_name=defaults['grb_name'], no_int_cal=defaults['no_int_cal'], det_cut=defaults['det_cut']
):
    global magtype
    global parallel

    start_time = dt.now()

    set_vizier_mirror()

    magtype = defaults['magtype']
    parallel = defaults['parallel']
    grb_flag = False

    comp_lvl = 0.3

    try:
        directory = os.path.dirname(full_filename)
    except TypeError:
        raise Exception('-filepath not specified!')
    if directory == '':
        directory = '.'
    directory = directory + os.path.sep
    name = os.path.basename(full_filename)

    if grb_ra or grb_dec or grb_radius != defaults['grb_radius']:
        if not grb_dec or not grb_ra:
            print('\n**WARNING: either -grb_ra or -grb_dec parameters are NOT recognized!**  '
                      'Photometry will continue... but you must specify both to use w/ GRB functionality!\n')
            if grb_only:
                raise ValueError('\n**-grb_only flag used w/o correct parameters! Need -grb_ra and -grb_dec!**')

        grb_flag = True

    grb_thresh = grb_rad_convert(grb_radius)

    if grb_only:
        os.chdir(directory)
        data, header, w, raImage, decImage, bulge, det_thresh, chip = img(directory, name, crop)
        Q, chosen_survey, mag_low_cutoff = query(raImage, decImage, band, w, data, crop, comp_lvl, survey, given_catalog, mag_low_lim,
                                                 mag_high_lim, bulge)

        psf_check = name.replace('.fits', '.fits.psf')
        if not os.path.isfile(psf_check):
            raise FileNotFoundError(f'\n{psf_check} file not found! This indicates PSF photom data products '
                                    f'have not been kept or PSF photom has not been run! Run GRB photometry again *W/O* '
                                    f'the -grb_only flag!')

        psfcatalogName = name.replace('.fits', '.photom.cat')
        good_cat_stars, cleanPSFSources, PSFSources, idx_psfmass, idx_psfimage, massCatCoords, crop_cat_stars = (
            tables(Q, data, w, psfcatalogName, crop, given_catalog))
        colnames = good_cat_stars.colnames
        magcolname = colnames[2]
        good_cat_stars[magcolname] = ab_convert(good_cat_stars[magcolname], band=band, survey=chosen_survey)
        ab_cat_stars = good_cat_stars
        if grb_coordlist:
            GRB(grb_ra, grb_dec, name, chosen_survey, band, grb_thresh, massCatCoords, ab_cat_stars, directory, chip, grb_coordlist, grbname=grb_name)
        else:
            if grb_thresh > 60:
                # newsourcesearch(grb_ra, grb_dec, w, name, chosen_survey, band, massCatCoords, grb_thresh)
                newsourcesearch(grb_ra, grb_dec, name, chosen_survey, band, grb_thresh, massCatCoords, ab_cat_stars, directory, chip,
                                grbname=grb_name, mag_low_lim=mag_low_cutoff)
            else:
                GRB(grb_ra, grb_dec, name, chosen_survey, band, grb_thresh, massCatCoords, ab_cat_stars, directory, chip, grbname=grb_name)
    else:
        data, header, w, raImage, decImage, bulge, det_thresh, chip = img(directory, name, crop)
        Q, chosen_survey, mag_low_cutoff = query(raImage, decImage, band, w, data, crop, comp_lvl, survey, given_catalog, mag_low_lim,
                                                 mag_high_lim, bulge)
        if not grb_flag:
            psfcatalogName = sex1(name, det_cut=det_thresh)
        else:
            catalogName = sex1(name, det_cut=det_thresh, grb_flag=grb_flag)
            psfex(catalogName, band, data, crop)
            psfcatalogName = sex2(name, det_cut=det_thresh, catalogName=catalogName)
        good_cat_stars, cleanPSFSources, PSFSources, idx_psfmass, idx_psfimage, massCatCoords, crop_cat_stars = (
            tables(Q, data, w, psfcatalogName, crop, given_catalog))
        if len(idx_psfimage) == 0:
            raise ValueError('No crossmatches found!  Cannot continue with photometry!  Is there something wrong with the image, '
                  'source catalogs, or psf model?  If those all seem normal, perhaps the image has had pixel values scaled'
                  'to uJy.  The current setup only applies the uJy/ADU conv factor as a header card, rerun the image '
                  'stacking and try again!')
        else:
            if no_int_cal:
                (cleanPSFSources, PSFsources, psfweights_noclip, psf_clipped, ab_cat_stars,
                 aperweights_noclip, aper_clipped_all, autoweights_noclip, auto_clipped) = (
                    zeropt(good_cat_stars, crop_cat_stars, cleanPSFSources, PSFSources, idx_psfmass, idx_psfimage, name, band,
                           chosen_survey,
                           sigma, data, crop))

                if no_plots:
                    photometric_fit_calc(cleanPSFSources, band, ab_cat_stars, idx_psfmass, idx_psfimage,
                                             psfweights_noclip, psf_clipped, sigma, aperweights_noclip,
                                             aper_clipped_all,
                                             autoweights_noclip, auto_clipped)
                else:
                    photometry_plots(cleanPSFSources, PSFsources, data, name, chosen_survey,
                                                                 band, ab_cat_stars, idx_psfmass,
                                                                 idx_psfimage, psfweights_noclip, psf_clipped, sigma,
                                                                 aperweights_noclip, aper_clipped_all
                                                                 , autoweights_noclip, auto_clipped)

                if grb_ra:
                    keep = True
                    if grb_dec is None:
                        print('Only GRB RA is found, GRB Dec is None!  Cant conduct grb analysis, '
                              'make sure the -grb_dec flag is correctly formatted!')
                    if grb_thresh > 60:
                        newsourcesearch(grb_ra, grb_dec, name, chosen_survey, band, grb_thresh, massCatCoords, ab_cat_stars,
                            directory, chip, grbname=grb_name, mag_low_lim=mag_low_cutoff)
                    else:
                        GRB(grb_ra, grb_dec, name, chosen_survey, band, grb_thresh, massCatCoords, ab_cat_stars,
                            directory, chip, grbname=grb_name)
                elif grb_coordlist:
                    GRB(grb_ra, grb_dec, name, chosen_survey, band, grb_thresh, massCatCoords, ab_cat_stars, directory,
                        chip, grb_coordlist, grbname=grb_name)

            else:
                if grb_flag:
                    parallel = True

                    print('\nMultiprocessing mag intercept calibration!'
                          '\nResults will be written to log files!')

                    sync_queue = Queue()
                    continue_queue = Queue()

                    auto_calib = Process(target=full_int_calibration,
                                         args=(
                                             name, directory, band, chip, crop, data, sigma, Q, chosen_survey,
                                             given_catalog, good_cat_stars, crop_cat_stars, cleanPSFSources, PSFSources, massCatCoords,
                                             idx_psfmass, idx_psfimage, mag_high_lim, mag_low_cutoff, grb_ra, grb_dec,
                                             grb_coordlist, grb_name, grb_thresh, comp_lvl, sync_queue, continue_queue
                                         ),
                                         kwargs=dict(
                                             no_plots=no_plots, input_magtype="AUTO", input_parallel=parallel
                                         )
                    )

                    psf_calib = Process(target=full_int_calibration,
                                        args=(
                                            name, directory, band, chip, crop, data, sigma, Q, chosen_survey,
                                            given_catalog, good_cat_stars, crop_cat_stars, cleanPSFSources, PSFSources, massCatCoords,
                                            idx_psfmass, idx_psfimage, mag_high_lim, mag_low_cutoff, grb_ra, grb_dec,
                                            grb_coordlist, grb_name, grb_thresh, comp_lvl, sync_queue, continue_queue
                                        ),
                                        kwargs=dict(
                                            no_plots=no_plots, input_magtype="PSF", input_parallel=parallel
                                        )
                    )

                    aper_calib = Process(target=full_int_calibration,
                                        args=(
                                            name, directory, band, chip, crop, data, sigma, Q, chosen_survey,
                                            given_catalog, good_cat_stars, crop_cat_stars, cleanPSFSources, PSFSources, massCatCoords,
                                            idx_psfmass, idx_psfimage, mag_high_lim, mag_low_cutoff, grb_ra, grb_dec,
                                            grb_coordlist, grb_name, grb_thresh, comp_lvl, sync_queue, continue_queue
                                        ),
                                        kwargs=dict(
                                            no_plots=no_plots, input_magtype="APER", input_parallel=parallel
                                        )
                    )

                    processes = {
                        'auto': auto_calib,
                        'psf': psf_calib,
                        'aper': aper_calib
                    }

                    auto_calib.start()
                    time.sleep(1)
                    psf_calib.start()
                    time.sleep(1)
                    aper_calib.start()

                    try:
                        # collect sync signals, only count signals from surviving processes
                        received = 0
                        deadline = time.time() + 120
                        while received < len(processes):
                            remaining = deadline - time.time()
                            if remaining <= 0:
                                break
                            try:
                                sync_queue.get(timeout=remaining)
                                received += 1
                            except Empty:
                                break

                        # dead & alive process filtering
                        alive_processes = {name: proc for name, proc in processes.items() if proc.is_alive()}
                        dead_processes = {name: proc for name, proc in processes.items() if not proc.is_alive()}

                        if dead_processes:
                            print(
                                f"Warning: the following calibration(s) did not reach sync point: {list(dead_processes.keys())}")

                        if not alive_processes:
                            raise RuntimeError(
                                "All calibration processes died before reaching sync point! Check log files!")

                        # combine catalogs into single ecsv w/ surviving mag types
                        comb_mag_catalogs(directory, name, chosen_survey)

                        # only signal and join surviving processes
                        for proc in alive_processes.values():
                            continue_queue.put("continue")
                        for proc in alive_processes.values():
                            proc.join()

                        # cull the dead
                        for proc in dead_processes.values():
                            proc.terminate()
                            proc.join()

                        # check exit codes of processes that ran to completion
                        failed = {name: proc.exitcode for name, proc in alive_processes.items() if proc.exitcode != 0}
                        if failed:
                            raise RuntimeError(
                                f"The following calibration(s) failed with non-zero exit codes: {failed}")

                    except (TimeoutError, RuntimeError) as e:
                        print(f"Error in calibration processes: {e}")
                        for proc in processes.values():
                            proc.terminate()
                            proc.join()
                        raise
                    except KeyboardInterrupt:
                        print(f"Calibration processes keyboard interrupted!")
                        for proc in processes.values():
                            proc.terminate()
                            proc.join()
                        raise

                    colnames = good_cat_stars.colnames
                    if len(colnames) > 6:
                        magcolname = f'{band}MAG_{MAGTYPES[magtype]}'
                    else:
                        magcolname = colnames[2]

                    if grb_flag:
                        print(' Regenerating survey star catalog for GRB functions.')
                        good_cat_stars[magcolname] = ab_convert(good_cat_stars[magcolname], band=band,
                                                                survey=chosen_survey)
                        ab_cat_stars = good_cat_stars

                        grb_arg_dict = dict(
                            ra=grb_ra,
                            dec=grb_dec,
                            imageName=name,
                            survey=chosen_survey,
                            band=band,
                            thresh=grb_thresh,
                            massCatCoords=massCatCoords,
                            good_cat_stars=ab_cat_stars,
                            directory=directory,
                            chip=chip,
                            coordlist=grb_coordlist,
                            grbname=grb_name,
                            mag_low_lim=mag_low_cutoff
                        )

                        keep = True
                        if grb_dec is None:
                            print('Only GRB RA is found, GRB Dec is None!  Cant conduct grb analysis, '
                                  'make sure the -grb_dec flag is correctly formatted!')
                        if grb_thresh > 60:
                            newsourcesearch(**grb_arg_dict)
                        else:
                            GRB(**grb_arg_dict)

                else:
                    parallel = True

                    print('\nMultiprocessing mag intercept calibration (only AUTO & APER)!'
                          '\nResults will be written to log files!')

                    sync_queue = Queue()
                    continue_queue = Queue()

                    auto_calib = Process(target=full_int_calibration,
                                         args=(
                                             name, directory, band, chip, crop, data, sigma, Q, chosen_survey,
                                             given_catalog, good_cat_stars, crop_cat_stars, cleanPSFSources, PSFSources, massCatCoords,
                                             idx_psfmass, idx_psfimage, mag_high_lim, mag_low_cutoff, grb_ra, grb_dec,
                                             grb_coordlist, grb_name, grb_thresh, comp_lvl, sync_queue, continue_queue
                                         ),
                                         kwargs=dict(
                                             no_plots=no_plots, input_magtype="AUTO", input_parallel=parallel
                                         )
                    )

                    aper_calib = Process(target=full_int_calibration,
                                        args=(
                                            name, directory, band, chip, crop, data, sigma, Q, chosen_survey,
                                            given_catalog, good_cat_stars, crop_cat_stars, cleanPSFSources, PSFSources, massCatCoords,
                                            idx_psfmass, idx_psfimage, mag_high_lim, mag_low_cutoff, grb_ra, grb_dec,
                                            grb_coordlist, grb_name, grb_thresh, comp_lvl, sync_queue, continue_queue
                                        ),
                                        kwargs=dict(
                                            no_plots=no_plots, input_magtype="APER", input_parallel=parallel
                                        )
                    )

                    processes = {
                        'auto': auto_calib,
                        'aper': aper_calib
                    }

                    auto_calib.start()
                    time.sleep(1)
                    aper_calib.start()

                    try:
                        # collect sync signals, only count signals from surviving processes
                        received = 0
                        deadline = time.time() + 120
                        while received < len(processes):
                            remaining = deadline - time.time()
                            if remaining <= 0:
                                break
                            try:
                                sync_queue.get(timeout=remaining)
                                received += 1
                            except Empty:
                                break

                        # dead & alive process filtering
                        alive_processes = {name: proc for name, proc in processes.items() if proc.is_alive()}
                        dead_processes = {name: proc for name, proc in processes.items() if not proc.is_alive()}

                        if dead_processes:
                            print(
                                f"Warning: the following calibration(s) did not reach sync point: {list(dead_processes.keys())}")

                        if not alive_processes:
                            raise RuntimeError(
                                "All calibration processes died before reaching sync point! Check log files!")

                        # combine catalogs into single ecsv w/ surviving mag types
                        comb_mag_catalogs(directory, name, chosen_survey)

                        # only signal and join surviving processes
                        for proc in alive_processes.values():
                            continue_queue.put("continue")
                        for proc in alive_processes.values():
                            proc.join()

                        # cull the dead
                        for proc in dead_processes.values():
                            proc.terminate()
                            proc.join()

                        # check exit codes of processes that ran to completion
                        failed = {name: proc.exitcode for name, proc in alive_processes.items() if proc.exitcode != 0}
                        if failed:
                            raise RuntimeError(
                                f"The following calibration(s) failed with non-zero exit codes: {failed}")

                    except (TimeoutError, RuntimeError) as e:
                        print(f"Error in calibration processes: {e}")
                        for proc in processes.values():
                            proc.terminate()
                            proc.join()
                        raise
                    except KeyboardInterrupt:
                        print(f"Calibration processes keyboard interrupted!")
                        for proc in processes.values():
                            proc.terminate()
                            proc.join()
                        raise

                    colnames = good_cat_stars.colnames
                    if len(colnames) > 6:
                        magcolname = f'{band}MAG_{MAGTYPES[magtype]}'
                    else:
                        magcolname = colnames[2]

                    if grb_flag:
                        print(' Regenerating survey star catalog for GRB functions.')
                        good_cat_stars[magcolname] = ab_convert(good_cat_stars[magcolname], band=band,
                                                                survey=chosen_survey)
                        ab_cat_stars = good_cat_stars

                        grb_arg_dict = dict(
                            ra=grb_ra,
                            dec=grb_dec,
                            imageName=name,
                            survey=chosen_survey,
                            band=band,
                            thresh=grb_thresh,
                            massCatCoords=massCatCoords,
                            good_cat_stars=ab_cat_stars,
                            directory=directory,
                            chip=chip,
                            coordlist=grb_coordlist,
                            grbname=grb_name,
                            mag_low_lim=mag_low_cutoff
                        )

                        keep = True
                        if grb_dec is None:
                            print('Only GRB RA is found, GRB Dec is None!  Cant conduct grb analysis, '
                                  'make sure the -grb_dec flag is correctly formatted!')
                        if grb_thresh > 60:
                            newsourcesearch(**grb_arg_dict)
                        else:
                            GRB(**grb_arg_dict)

                    # full_int_calibration(
                    #     name, directory, band, chip, crop, data, sigma, Q, chosen_survey,
                    #     given_catalog, good_cat_stars, crop_cat_stars, cleanPSFSources, PSFSources, massCatCoords,
                    #     idx_psfmass, idx_psfimage, mag_high_lim, mag_low_cutoff, grb_ra, grb_dec,
                    #     grb_coordlist, grb_name, grb_thresh, comp_lvl,
                    #     no_plots=no_plots, input_magtype="AUTO", input_parallel=False
                    # )

            if not keep:
                removal(directory)

    end_time = dt.now()
    print('\nFull photometric processing time:', (end_time - start_time).total_seconds())


def main():
    parser = argparse.ArgumentParser(
        description='runs sextractor and psfex on swarped img to get psf fit photometry, then '
                    'outputs ecsv w/ corrected mags')
    parser.add_argument('-exp_query', action='store_true', help='optional flag, exports ecsv of astroquery results '
                                                                'along with photometry', default=defaults['exp_query'])
    # parser.add_argument('-exp_query_only', action='store_true',
    #                     help='optional flag, use if photom is run already to only generate query results')
    parser.add_argument('-no_plots', action='store_true',
                        help='optional flag, stops creation of mag comparison plot betw. PRIME and survey, '
                             'along with residual plot w/ statistics, lim mag plot', default=defaults['no_plots'])
    parser.add_argument('-keep', action='store_true',
                        help='optional flag, use if you DONT want to remove intermediate products after getting photom,'
                             ' i.e. the ".cat" and ".psf" files', default=defaults['keep'])
    parser.add_argument('-grb_only', action='store_true',
                        help='optional flag, use if running -grb again on already created catalog',
                        default=defaults['grb_only'])
    parser.add_argument('-filepath', type=str, help='[str], full file path of stacked image, can also place just'
                                                    'filename and it will default to current directory',
                        default=defaults['filepath'])
    parser.add_argument('-band', type=str, help='[str], band, ex. "J"', default=defaults['band'])
    parser.add_argument('-survey', type=str,
                        help='[str], *NOW OPTIONAL* manually specify which survey to query, choose from VHS, 2MASS'
                             ', VIKING, Skymapper, SDSS, UKIDSS, & DES.  If you leave out this arg, it will automatically'
                             'pick a survey from the above list depending on the area and coverage.',
                        default=defaults['survey'])
    parser.add_argument('-crop', type=int, help='[int], # of pixels from edge of image to crop, default = 300',
                        default=defaults["crop"])
    parser.add_argument('-sigma', type=float, help='[float], # of sigma w/ which to sigma clip for '
                                                  'zero point calculation, default = 3',
                        default=defaults["sigma_photom"])
    parser.add_argument('-catalog', type=str, help='[str], optional field to supply an already generated'
                                                   'catalog for photometry INSTEAD of querying, put in full file path.',
                        default=defaults["catalog"])
    parser.add_argument('-mag_low', type=float, help='[float], Lower mag cutoff for survey query & crossmatch'
                                                        ' settings default = 12.5',
                        default=defaults["mag_low"])
    parser.add_argument('-mag_high', type=float, help='[float], Higher mag cutoff for survey query & crossmatch'
                                                        ' settings default = 21, currently only applies to DES & Skymapper',
                        default=defaults["mag_high"])
    parser.add_argument('-grb_ra', type=str, help='[str], RA for GRB source, either in hh:mm:ss or decimal'
                                                  '*NOTE* When using sexagesimal, use "-grb_ra=value_here" NOT "-grb_ra '
                                                  'value_here", as argparse doesnt like negative sexagesimals',
                        default=defaults["grb_ra"])
    parser.add_argument('-grb_dec', type=str, help='[str], DEC for GRB source, either in dd:mm:ss or decimal'
                                                   '*NOTE* When using sexagesimal, use "-grb_dec=value_here" NOT "-grb_dec '
                                                  'value_here", as argparse doesnt like negative sexagesimals',
                        default=defaults["grb_dec"])
    parser.add_argument('-grb_coordlist', type=str, nargs='+',
                        help='[float] Used to check multiple GRB locations.  Input RA and DECs of locations '
                             'with the format: -coordlist 123,45 -123,-45 etc..  *DONT USE -RA '
                             '& -DEC BUT INCLUDE -grb_radius*', default=defaults["grb_coordlist"])
    parser.add_argument('-grb_radius', type=str,
                        help='[str], # of arcsec diameter to search for GRB, default = 4.0".  You can specify arcsec,'
                             ' arcmin, or deg w/ an underscore.  Ex. "-grb_radius 3_arcmin" will specify an area of 3 '
                             'arcminutes.  If just a number is applied, it defaults to arcsec.',
                        default=defaults["grb_radius"])
    parser.add_argument('-grb_name', type=str,
                        help='[str] optional name for grb-related data products, default = "GRB"',
                        default=defaults["grb_name"])
    parser.add_argument('-no_int_cal', action='store_true',
                        help='optional flag, use to STOP photometric fit intercept calibration, taking only the initial calc.',
                        default=defaults["no_int_cal"])
    parser.add_argument('-det_cut', type=float, help='[float], num of median image sigma to cut off sources'
                                                     ' (ex. det_thresh of 2 => cutoff = med - 2*sigma',
                        default=defaults["det_cut"])

    args, unknown = parser.parse_known_args()
    # print(args)
    # print(unknown)

    photometry(args.filepath, args.band, args.crop, args.sigma, args.catalog, args.survey, args.mag_low,
               args.mag_high, args.no_plots, args.keep,
               args.grb_only, args.grb_ra, args.grb_dec, args.grb_coordlist, args.grb_radius, args.grb_name,
               args.no_int_cal, args.det_cut)


if __name__ == "__main__":
    main()
