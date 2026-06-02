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
import psycopg2
from scipy.stats import skew
from scipy import odr
import warnings
from datetime import datetime as dt

from photometrus.settings import (gen_config_file_name, bulge_checker, PHOTOMETRY_MAG_LOWER_LIMIT, PHOTOMETRY_MAG_UPPER_LIMIT,
                                  PHOTOMETRY_QUERY_WIDTH, PHOTOMETRY_QUERY_CATALOGS, PHOTOMETRY_LIM_MAGS,
                                  AB_OFFSET_DICT, PRIME_FILTERS_DICT, get_weight_thresh, set_vizier_mirror, CHIP_ZPS,
                                  MAGTYPES, local_query_box)

from photometrus.utils.defaults import PROCESSING_DEFAULTS as defaults
from photometrus.query import query
from photometrus.photometry.photom_utils import (open_fits_robust, log_output, get_table_from_ldac, get_band,
                                                 grb_rad_convert, ab_convert, removal)
from photometrus.photometry.GRB import GRB, newsourcesearch
from photometrus.photometry.plots import single_plots, photometry_plots

# %%
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings(action="ignore", module="scipy", message="^One or more")
warnings.filterwarnings(action="ignore", module="numpy", message="Warning: 'partition' will ignore the 'mask' of the MaskedColumn.")

magtype = defaults['magtype']
parallel = defaults['parallel']


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
# run sextractor on swarped img to find sources


def sex1(imageName, det_cut, grb_flag=False, sx_cfg=defaults['sx_cfg']):
    print('Running sextractor on img to initially find sources...')
    aper_str = ''

    configFile = gen_config_file_name(sx_cfg)
    if grb_flag:
        paramName = gen_config_file_name('tempsource.param')
        catalogName = imageName + '.cat'
    else:
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

    if sx_cfg != defaults['sx_cfg']:
        print(f' Custom sxtrctr cfg specified: {sx_cfg}')

    weightName = 'weight'+imageName[5:]
    if os.path.isfile(weightName):
        # imghdr = fits.getheader(imageName)
        # if 'BUNIT' in imghdr:
        #     scale_fac = imghdr['CONV_FAC']
        # else:
        scale_fac = 1
        weightdata = fits.getdata(weightName)
        # weight map quality check (catch for all 0 or nan image)
        if not np.any(np.nan_to_num(weightdata)):
            print(' *WARNING* ISSUE W/ WEIGHT MAP!  Input weight map is all zero or nan, cannot use!\n'
                  ' Continuing w/o weight map...')
            weight_med = -100
            detect_cutoff = -100
            weight_map_str = ''
        else:
            weightdata = weightdata / scale_fac**2
            weight_med = np.nanmedian(weightdata)
            weight_std = np.nanstd(weightdata)
            detect_cutoff = weight_med - (weight_std * det_cut)
            weight_map_str = f'-WEIGHT_TYPE MAP_WEIGHT -WEIGHT_THRESH {detect_cutoff} -WEIGHT_IMAGE {weightName}'
            print(' Including weight map!')

        with fits.open(weightName, mode='update') as hdu:
            whdr = hdu[0].header
            whdr.set('MEDIAN', weight_med, 'Median of weight image', after='EQUINOX')
            whdr.set('DET_CUT', detect_cutoff, 'Pix value cutoff for source detection', after='MEDIAN')
            hdu.close()
        try:
            command = (f'sex %s -c %s -CATALOG_NAME %s {weight_map_str} -PARAMETERS_NAME %s {aper_str}'
                       %
                       (f'{imageName}', configFile, catalogName, paramName))
            # print('Executing command: %s' % command)
            subprocess.run(command.split(), check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as err:
            print(f"SExtractor failed with exit code {err.returncode}")
            print(f"STDERR:\n{err.stderr}")
            print(f"STDOUT:\n{err.stdout}")
            raise Exception('Sextractor failed to run, is the stacked image quality adequate?')
    else:
        try:
            command = (f'sex %s -c %s -CATALOG_NAME %s -PARAMETERS_NAME %s {aper_str}' %
                       (f'{imageName}', configFile, catalogName, paramName))
            # print('Executing command: %s' % command)
            subprocess.run(command.split(), check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as err:
            print(f"SExtractor failed with exit code {err.returncode}")
            print(f"STDERR:\n{err.stderr}")
            print(f"STDOUT:\n{err.stdout}")
            raise Exception('Sextractor failed to run, is the stacked image quality adequate?')
    return catalogName


# %%
# run psfex on sextractor LDAC from previous step


def psfex(catalogName, band, data, crop):
    print('Running PSFex on sextrctr catalogue to generate psf for stars in the img...')
    psfConfigFile = gen_config_file_name('default.psfex')

    if catalogName.startswith(f'coadd.Open-{band}'):
        psfImageName = 'PSF' + catalogName[5:-4]
    else:
        psfImageName = 'PSF' + catalogName

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
        subprocess.run(command.split(), check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as err:
        print(f"PSFeX failed with exit code {err.returncode}")
        print(f"STDERR:\n{err.stderr}")
        print(f"STDOUT:\n{err.stdout}")
        raise Exception('Is there a problem with the PSF model or input sextractor catalog?')
    os.rename('PSF_' + catalogName[:-4] + '.fits', psfImageName)


# %%
# feed generated psf model back into sextractor w/ diff param (or could use that param from the start but its slower)


def sex2(imageName, det_cut, catalogName, sx_cfg=defaults['sx_cfg']):
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
    psfcatalogName = imageName.replace(os.path.splitext(imageName)[1], '.photom.cat')
    configFile = gen_config_file_name(sx_cfg)
    psfparamName = gen_config_file_name('photomPSF.param')
    if sx_cfg != defaults['sx_cfg']:
        print(f' Custom sxtrctr cfg specified: {sx_cfg}')

    weightName = 'weight' + imageName[5:]
    if os.path.isfile(weightName):
        # imghdr = fits.getheader(imageName)
        # if 'BUNIT' in imghdr:
        #     scale_fac = imghdr['CONV_FAC']
        # else:
        scale_fac = 1
        weightdata = fits.getdata(weightName)

        if not np.any(np.nan_to_num(weightdata)):
            print(' *WARNING* ISSUE W/ WEIGHT MAP!  Input weight map is all zero or nan, cannot use!\n'
                  ' Continuing w/o weight map...')
            weight_map_str = ''
        else:
            weightdata = weightdata / scale_fac**2
            weight_med = np.nanmedian(weightdata)
            weight_std = np.nanstd(weightdata)
            detect_cutoff = weight_med - (weight_std * det_cut)
            weight_map_str = f'-WEIGHT_TYPE MAP_WEIGHT -WEIGHT_THRESH {detect_cutoff} -WEIGHT_IMAGE {weightName}'
            print(' Including weight map!')

        try:
            # We are supplying SExtactor with the PSF model with the PSF_NAME option
            command = (f'sex {imageName} -c {configFile} -CATALOG_NAME {psfcatalogName} {weight_map_str} -PSF_NAME {psfName} '
                       f'-PARAMETERS_NAME {psfparamName} {aper_str}')
            # print("Executing command: %s" % command)
            subprocess.run(command.split(), check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as err:
            print(f"SExtractor failed with exit code {err.returncode}")
            print(f"STDERR:\n{err.stderr}")
            print(f"STDOUT:\n{err.stdout}")
            raise Exception('Is there a problem with the sextractor configs or PSF model?')
    else:
        try:
            # We are supplying SExtactor with the PSF model with the PSF_NAME option
            command = (f'sex {imageName} -c {configFile} -CATALOG_NAME {psfcatalogName} -PSF_NAME {psfName} '
                       f'-PARAMETERS_NAME {psfparamName} {aper_str}')
            # print("Executing command: %s" % command)
            subprocess.run(command.split(), check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as err:
            print(f"SExtractor failed with exit code {err.returncode}")
            print(f"STDERR:\n{err.stderr}")
            print(f"STDOUT:\n{err.stdout}")
            raise Exception('Is there a problem with the sextractor configs or PSF model?')
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

        if 'Mclass' in colnames or 'mergedClass' in colnames:
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


#%% old fit calculation function, MOSTLY DEPRECIATED (only used w/ -no_int_cal flag)

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

    if grb_ra:
        psfcatalogName = name.replace(os.path.splitext(name)[1], '.photom.cat')
    else:
        psfcatalogName = f"{name}.photom.cat"

    good_cat_stars, cleanPSFSources, PSFSources, idx_psfmass, idx_psfimage, massCatCoords, crop_cat_stars = (
        tables(Q, data, w, psfcatalogName, crop, given_catalog))

    cleanPSFSources, PSFSources, weights_noclip, clipped, ab_cat_stars = single_zeropt(
        good_cat_stars, crop_cat_stars, cleanPSFSources, PSFSources, idx_psfmass, idx_psfimage, name, band, survey, sigma,
        data, crop)

    if make_plots:
        slope, intercept, int_err = single_plots(cleanPSFSources, PSFSources, data, name, survey, band, good_cat_stars,
                                                 idx_psfmass, idx_psfimage, sigma, weights_noclip, clipped, crop,
                                                 magtype=MAGTYPES[magtype])
        if sync_queue:
            sync_queue.put("Calib. done, combining .ecsv files")
            continue_queue.get()

        if grb_ra:
            if grb_radius > 60:
                newsourcesearch(grb_ra, grb_dec, name, chosen_survey, band, grb_radius, massCatCoords, ab_cat_stars, directory,
                    chip, grbname=grb_name, mag_low_lim=mag_low_lim, magtype=MAGTYPES[magtype])
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
            crop=crop,
            magtype=MAGTYPES[magtype]
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
                mag_low_lim=PHOTOMETRY_MAG_LOWER_LIMIT,
                magtype=MAGTYPES[magtype]
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
            mag_low_lim=PHOTOMETRY_MAG_LOWER_LIMIT,
            magtype=MAGTYPES[magtype]
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
            crop=crop,
            magtype=MAGTYPES[magtype]
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
                                                                                           magtype=MAGTYPES[magtype],
                                                                                           parallel=parallel,
                                                                                           sync_signal=sync_signal,
                                                                                           sync_queue=sync_queue,
                                                                                           continue_queue=continue_queue,
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
                                                                       comp_lvl=comp_lvl, magtype=MAGTYPES[magtype], parallel=parallel,
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
                                                                                           magtype=MAGTYPES[magtype],
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
                                                                 comp_lvl=comp_lvl, make_plots=True, magtype=MAGTYPES[magtype],
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

#%%


def photometry(
        full_filename=defaults['filepath'], band=defaults['band'], crop=defaults['crop'], sigma=defaults['sigma_photom'], given_catalog=defaults['catalog'], survey=defaults['survey'],
        mag_low_lim=defaults['mag_low'], mag_high_lim=defaults['mag_high'], no_plots=defaults['no_plots'],
        keep=defaults['keep'], grb_only=defaults['grb_only'], grb_ra=defaults['grb_ra'], grb_dec=defaults['grb_dec'], grb_coordlist=defaults['grb_coordlist'],
        grb_radius=defaults['grb_radius'], grb_name=defaults['grb_name'], no_int_cal=defaults['no_int_cal'], det_cut=defaults['det_cut'],
        sx_cfg=defaults['sx_cfg']
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

    # try:
    #     band = get_band(full_filename)
    # except KeyError as e:
    #     print(f'Band could not be found in image header, going with default = {band}'
    #           f'\n{e}')

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

        psf_check = name.replace(os.path.splitext(name)[1], f'{os.path.splitext(name)[1]}.psf')
        if not os.path.isfile(psf_check):
            raise FileNotFoundError(f'\n{psf_check} file not found! This indicates PSF photom data products '
                                    f'have not been kept or PSF photom has not been run! Run GRB photometry again *W/O* '
                                    f'the -grb_only flag!')

        psfcatalogName = name.replace(os.path.splitext(name)[1], '.photom.cat')
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
            psfcatalogName = sex1(name, det_cut=det_thresh, sx_cfg=sx_cfg)
        else:
            catalogName = sex1(name, det_cut=det_thresh, grb_flag=grb_flag, sx_cfg=sx_cfg)
            psfex(catalogName, band, data, crop)
            psfcatalogName = sex2(name, det_cut=det_thresh, catalogName=catalogName, sx_cfg=sx_cfg)
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

                    removal(directory=directory, end_names='.lock')

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
                            mag_low_lim=mag_low_cutoff,
                            magtype=MAGTYPES[magtype]
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

                    removal(directory=directory, end_names='.lock')

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
                            mag_low_lim=mag_low_cutoff,
                            magtype=MAGTYPES[magtype]
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
    parser.add_argument('-filepath', type=str, help='[str], *REQUIRED* full file path of stacked image, can also place just'
                                                    'filename and it will default to current directory',
                        default=defaults['filepath'])
    parser.add_argument('-band', type=str, help='[str], *REQUIRED* band of observation, default = J',
                        default=defaults['band'])
    parser.add_argument('-survey', type=str,
                        help='[str], *NOW OPTIONAL* manually specify which survey to query, choose from VHS, 2MASS'
                             ', VIKING, Skymapper, SDSS, UKIDSS, & DES.  If you leave out this arg, it will automatically'
                             ' pick a survey from the above list depending on the area and coverage.',
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
    parser.add_argument('-sx_cfg', type=str,
                        help='[str] optionally specify different sxtrctr config file to use for main source extraction,'
                             'must be in .prime/config/ directory, '
                             'default = sex2.config',
                        default=defaults["sx_cfg"])

    args, unknown = parser.parse_known_args()
    # print(args)
    # print(unknown)

    photometry(args.filepath, args.band, args.crop, args.sigma, args.catalog, args.survey, args.mag_low,
               args.mag_high, args.no_plots, args.keep,
               args.grb_only, args.grb_ra, args.grb_dec, args.grb_coordlist, args.grb_radius, args.grb_name,
               args.no_int_cal, args.det_cut, args.sx_cfg)


if __name__ == "__main__":
    main()
