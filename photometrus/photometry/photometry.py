"""
Calibrates photometry for stacked image
"""

import os
import sys
import re
import numpy as np
import numpy.ma as ma
import argparse
import astropy.units as u
from astroquery.vizier import Vizier
from astroquery.ipac.ned import Ned
from astropy.coordinates import Angle, SkyCoord
from astropy.wcs import WCS
from astropy.wcs import utils
from astropy.stats import sigma_clip, sigma_clipped_stats
from astropy.io import fits
from astropy.io import ascii
from astropy.nddata import Cutout2D
from regions import CircleSkyRegion
import base64
from astropy.table import Column
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
                                  AB_OFFSET_DICT)

from photometrus.utils.defaults import PROCESSING_DEFAULTS as defaults

# %%
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings(action="ignore", module="scipy", message="^One or more")


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
    # jy zero points
    # from 2MASS
    offset_dict = AB_OFFSET_DICT

    if survey == '2MASS':
        zp_dict = {'zp_J': 1594, 'zp_H': 1024}
        pick_zp = 'zp_' + band

        for k, v in zp_dict.items():
            if k == pick_zp:
                zp = v

        flx = zp * 10 ** (-mag / 2.5)
        ab_mag = -2.5 * np.log10(flx / 3631)
        if revert:
            # print('Temporarily reverting PRIME mags to Vega to compare to survey (for GRB or plotting)')
            v_flx = 3631 * 10 ** (-mag / 2.5)
            vega_mag = -2.5 * np.log10(v_flx / zp)
            return vega_mag

    elif 'UKIDSS' in survey:
        # for ukirt conversion: https://adsabs.harvard.edu/full/2006MNRAS.367..454H
        for k, v in offset_dict['UKIDSS'].items():
            if k == band:
                offset = v
        ab_mag = mag+offset
        if revert:
            # print('Temporarily reverting PRIME mags to Vega to compare to survey (for GRB or plotting)')
            vega_mag = mag - offset
            return vega_mag

    elif (survey == 'DES_Z' or survey == 'DES_Y' or survey == 'Skymapper' or survey == 'SDSS' or survey == 'PanSTARRS'
          or survey == 'PanSTARRS_Z'):
        print('Survey %s is already reported in AB mag, no offset required.' % survey)
        ab_mag = mag

    else:
        # for vista (AB-Vega offsets): https://www.aanda.org/articles/aa/full_html/2015/03/aa24973-14/T3.html

        for k, v in offset_dict['VISTA'].items():
            if k == band:
                offset = v
        ab_mag = mag+offset
        if revert:
            # print('Temporarily reverting PRIME mags to Vega to compare to survey (for GRB or plotting)')
            vega_mag = mag - offset
            return vega_mag

    return ab_mag

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

    return data, header, w, raImage, decImage, bulge


# %%
# Use astroquery to get catalog search


def query(raImage, decImage, band, survey=None, given_catalog_path=None, mag_lower_lim=None, mag_upper_lim=None, bulge=False):
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
        Q = Q[Q['%sMAG_PSF' % band] > mag_low_cutoff]
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
                                '%spmag' % band.lower(), 'e_%spmag' % band.lower()])
            try:
                result = v.query_region(coords, width=str(checkwidth) + 'm', catalog=[f[1] for f in catalogs])
                test = result[0]
            except IndexError:
                print('Sadly, no current surveys available in current area in %s band' % band)
                sys.exit('No surveys available.')

            keys = result.format_table_list()

            print('%s, %s, Box Width: %s arcmin... '
                  '\n%s band initial query resulting in: \n%s' % (
                  frame_long_str, frame_lat_str, checkwidth, band, keys))

            keycheck = result.keys()

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

                            print('\nQuerying Vizier %s around %s, %s, boxwidth %.2f arcmin, mag lim of %s'
                                  % (catNum, frame_long_str, frame_lat_str, width, mag_lims))
                            try:
                                v = Vizier(columns=['%s' % cols[0], '%s' % cols[1], '%s' % cols[2], '%s' % cols[3]],
                                           column_filters={"%s" % cols[2]: mag_lims,
                                                           "%sFlag" % band.lower(): "<4",
                                                           "%sflags1" % band.lower(): "<16"
                                                           # "%sperrbits" % band: '<=16',
                                                           }, row_limit=-1)
                                Q = v.query_region(SkyCoord(coords, unit=(u.deg, u.deg)),
                                                   width=str(width) + 'm'
                                                   , catalog=k, cache=False, frame=chosen_frame)
                                if Q and len(Q[0]) > 0:
                                    print('Queried source total = ', len(Q[0]))
                                    break  # success
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
                    sys.exit('DES only supports Z and Y band!')
            else:
                print('No supported survey found, currently use either 2MASS, VHS, VIKING, Skymapper, SDSS, UKIDSS, or DES')
                sys.exit('No surveys found.')

    return Q, chosen_survey, mag_low_cutoff


# %%
# run sextractor on swarped img to find sources


def sex1(imageName, det_cut):
    print('Running sextractor on img to initially find sources...')
    configFile = gen_config_file_name('sex2.config')
    paramName = gen_config_file_name('tempsource.param')
    catalogName = imageName + '.cat'
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
        detect_cutoff = weight_med * det_cut
        with fits.open(weightName, mode='update') as hdu:
            whdr = hdu[0].header
            whdr.set('MEDIAN', weight_med, 'Median of weight image', after='EQUINOX')
            whdr.set('DET_CUT', detect_cutoff, 'Pix value cutoff for source detection', after='MEDIAN')
            hdu.close()
        try:
            print('Including weight map!')
            command = ('sex %s -c %s -CATALOG_NAME %s -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_THRESH %s -WEIGHT_IMAGE %s -PARAMETERS_NAME %s' %
                       (imageName, configFile, catalogName, detect_cutoff, weightName, paramName))
            # print('Executing command: %s' % command)
            subprocess.run(command.split(), check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except subprocess.CalledProcessError as err:
            print('Could not run sextractor with exit error %s' % err)
    else:
        try:
            command = ('sex %s -c %s -CATALOG_NAME %s -PARAMETERS_NAME %s' %
                       (imageName, configFile, catalogName, paramName))
            # print('Executing command: %s' % command)
            subprocess.run(command.split(), check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except subprocess.CalledProcessError as err:
            print('Could not run sextractor with exit error %s' % err)
    return catalogName


# %%
# run psfex on sextractor LDAC from previous step


def psfex(catalogName):
    print('Running PSFex on sextrctr catalogue to generate psf for stars in the img...')
    psfConfigFile = gen_config_file_name('default.psfex')
    psfImageName = 'PSF' + catalogName[5:-4]
    try:
        command = 'psfex %s -c %s -CHECKIMAGE_TYPE SNAPSHOTS -CHECKIMAGE_NAME %s' % (catalogName, psfConfigFile,
                                                                                     'PSF.fits')
        # print('Executing command: %s' % command)
        subprocess.run(command.split(), check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError as err:
        print('Could not run psfex with exit error %s' % err)
    os.rename('PSF_' + catalogName[:-4] + '.fits', psfImageName)


# %%
# feed generated psf model back into sextractor w/ diff param (or could use that param from the start but its slower)


def sex2(imageName, det_cut):
    print('Feeding psf model back into sextractor for fitting and flux calculation...')
    psfName = imageName + '.psf'
    psfcatalogName = imageName.replace('.fits', '.psf.cat')
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
        detect_cutoff = weight_med * det_cut
        try:
            # We are supplying SExtactor with the PSF model with the PSF_NAME option
            command = 'sex %s -c %s -CATALOG_NAME %s -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_THRESH %s -WEIGHT_IMAGE %s -PSF_NAME %s -PARAMETERS_NAME %s' % (
            imageName, configFile, psfcatalogName, detect_cutoff, weightName, psfName, psfparamName)
            # print("Executing command: %s" % command)
            subprocess.run(command.split(), check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except subprocess.CalledProcessError as err:
            print('Could not run sextractor with exit error %s' % err)
    else:
        try:
            # We are supplying SExtactor with the PSF model with the PSF_NAME option
            command = 'sex %s -c %s -CATALOG_NAME %s -PSF_NAME %s -PARAMETERS_NAME %s' % (
            imageName, configFile, psfcatalogName, psfName, psfparamName)
            # print("Executing command: %s" % command)
            subprocess.run(command.split(), check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except subprocess.CalledProcessError as err:
            print('Could not run sextractor with exit error %s' % err)
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

        good_cat_stars = given_cat[(given_cat['FLAGS'] == 0) & (given_cat['FLAGS_MODEL'] == 0) &
                               (given_cat['XMODEL_IMAGE'] < (max_x - crop)) & (given_cat['XMODEL_IMAGE'] > crop) &
                               (given_cat['YMODEL_IMAGE'] < (max_y) - crop) & (given_cat['YMODEL_IMAGE'] > crop)]

        psfsourceTable = get_table_from_ldac(psfcatalogName)
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
            (psfsourceTable['XMODEL_IMAGE'] < (max_x - crop)) & (psfsourceTable['XMODEL_IMAGE'] > crop)
            & (psfsourceTable['YMODEL_IMAGE'] < (max_y) - crop) & (psfsourceTable['YMODEL_IMAGE'] > crop) &
            (flux_radius >= 1 / 0.498)]

        cleanPSFSources = psfsourceTable[
            (psfsourceTable['FLAGS'] == 0) & (psfsourceTable['FLAGS_MODEL'] == 0) & (psfsourceTable['XMODEL_IMAGE']
            < (max_x - crop)) & (psfsourceTable['XMODEL_IMAGE'] > crop) & (psfsourceTable['YMODEL_IMAGE']
            < (max_y) - crop) & (psfsourceTable['YMODEL_IMAGE'] > crop) & (flux_radius >= 1 / 0.498)]

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
        good_cat_stars = Q[0][np.where(
            (mass_imCoords[0] > crop) & (mass_imCoords[0] < (max_x - crop)) & (mass_imCoords[1] > crop) & (
                        mass_imCoords[1] < (max_y - crop)))]
        print('Catalogue cropped, source total = ', len(good_cat_stars))

        psfsourceTable = get_table_from_ldac(psfcatalogName)
        if isinstance(psfsourceTable['FLUX_RADIUS'][0], np.ndarray):
            r50 = psfsourceTable['FLUX_RADIUS'][:, 0]
            r90 = psfsourceTable['FLUX_RADIUS'][:, 1]

            idx_r = psfsourceTable.colnames.index('FLUX_RADIUS')
            psfsourceTable['FLUX_RADIUS'] = r50
            psfsourceTable.add_column(r90, name='FLUX_RADIUS_90', index=idx_r + 1)
            flux_radius = psfsourceTable['FLUX_RADIUS_90']
        else:
            flux_radius = psfsourceTable['FLUX_RADIUS']

        PSFSources = psfsourceTable[
            (psfsourceTable['XMODEL_IMAGE'] < (max_x - crop)) & (psfsourceTable['XMODEL_IMAGE'] > crop)
            & (psfsourceTable['YMODEL_IMAGE'] < (max_y) - crop) & (psfsourceTable['YMODEL_IMAGE'] > crop)
            & (flux_radius >= 1 / 0.498)]

        cleanPSFSources = psfsourceTable[
            (psfsourceTable['FLAGS'] == 0) & (psfsourceTable['FLAGS_MODEL'] == 0) & (psfsourceTable['XMODEL_IMAGE']
                                                                                     < (max_x - crop)) & (
                        psfsourceTable['XMODEL_IMAGE'] > crop) & (psfsourceTable['YMODEL_IMAGE'] < (max_y) - crop)
            & (psfsourceTable['YMODEL_IMAGE'] > crop) & (flux_radius >= 1 / 0.498)]

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

        print('Found %d good cross-matches' % len(idx_psfmass))
    return good_cat_stars, cleanPSFSources, PSFSources, idx_psfmass, idx_psfimage, massCatCoords


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


def zeropt(good_cat_stars, cleanPSFSources, PSFSources, idx_psfmass, idx_psfimage, imageName, band, survey, sigma):
    colnames = good_cat_stars.colnames
    magcolname = colnames[2]
    magerrcolname = colnames[3]

    caterr = good_cat_stars[magerrcolname][idx_psfmass]
    primeerr = cleanPSFSources['MAGERR_POINTSOURCE'][idx_psfimage]

    comberr = np.sqrt(caterr ** 2 + primeerr ** 2)

    psfweights_noclip = 1 / (comberr ** 2)

    psfoffsets = ma.array(good_cat_stars[magcolname][idx_psfmass] - cleanPSFSources['MAG_POINTSOURCE'][idx_psfimage])
    psfoffsets = psfoffsets.data

    # 3 sigma clip
    psf_clipped = sigma_clip(psfoffsets, sigma=sigma)
    psfoffsets = psfoffsets[~psf_clipped.mask]
    psfweights = np.array(psfweights_noclip[~psf_clipped.mask])
    print('Zero point source offsets clipped by %s sigma, total clipped offset # = %s' % (sigma, len(psfoffsets)))

    # Compute statistics
    # zero_psfmean = np.average(psfoffsets, weights=psfweights)
    zero_psfmean = sum(psfoffsets * psfweights) / sum(psfweights)
    # zero_psfvar = np.average((psfoffsets - zero_psfmean) ** 2, weights=psfweights)
    # zero_psfstd = np.sqrt(zero_psfvar)
    zero_psfstd = np.sqrt(1 / sum(psfweights))

    # zero_psfmean, zero_psfmed, zero_psfstd = sigma_clipped_stats(psfoffsets)
    # print('PSF Mean ZP: %.2f\nPSF Median ZP: %.2f\nPSF STD ZP: %.2f'%(zero_psfmean, zero_psfmed, zero_psfstd))

    print('zp = %.4f, zp err = %.6f' % (zero_psfmean, zero_psfstd))
    # writing zp to header
    print('Writing ZP info to image header...')
    with fits.open(imageName, mode='update') as hdul:
        hdr = hdul[0].header
        hdr.set('ZP', zero_psfmean, 'Zero Point Offset', after='NINT')
        hdr.set('e_ZP', zero_psfstd, 'Zero Point Offset Error', after='ZP')
        hdr.set('N_CRSMCH', len(idx_psfimage), 'Number of Crossmatches', after='e_ZP')
        hdr.set('N_SRCS', len(PSFSources), 'Total PRIME Sources Number', after='N_CRSMCH')
        hdr.set('Survey', survey, 'Chosen Survey for Crossmatch', after='N_SRCS')
        hdul.close()

    # catalog for all detected sources
    psfmag = zero_psfmean + PSFSources['MAG_POINTSOURCE']
    psfmagerr = np.sqrt(PSFSources['MAGERR_POINTSOURCE'] ** 2 + zero_psfstd ** 2)
    # ab mag conversion
    print('Converting mags from Vega to AB for all sources! (if not already in AB)')
    psfmag = ab_convert(psfmag, band=band, survey=survey)

    psfmagcol = Column(psfmag, name='%sMAG_PSF' % band, unit=u.ABmag)
    psfmagerrcol = Column(psfmagerr, name='e_%sMAG_PSF' % band, unit=u.ABmag)
    PSFSources.add_column(psfmagcol)
    PSFSources.add_column(psfmagerrcol)
    PSFSources.remove_column('VIGNET')
    PSFSources['FLUX_RADIUS'] = PSFSources['FLUX_RADIUS'] * 0.498
    PSFSources['FLUX_RADIUS'].unit = u.arcsec
    if 'FLUX_RADIUS_90' in PSFSources.colnames:
        PSFSources['FLUX_RADIUS_90'] = PSFSources['FLUX_RADIUS_90'] * 0.498
        PSFSources['FLUX_RADIUS_90'].unit = u.arcsec
    print('Total PRIME source # = ', len(PSFSources))

    # real unit flux conversion
    psfflux = psfmagcol.to(u.microjansky)
    psffluxcol = Column(psfflux, name='FLUX_DENSITY', unit=u.microjansky)
    PSFSources.add_column(psffluxcol)

    # image / col data conversion to uJy

    conv_factor_all = psffluxcol / PSFSources['FLUX_POINTSOURCE']   # u = uJy / adu
    conv_factor = np.nanmedian(conv_factor_all)

    fluxerr_ujy = PSFSources['FLUXERR_POINTSOURCE'] * conv_factor
    fluxerrujycol = Column(fluxerr_ujy, name='E_FLUX_DENSITY', unit=u.microjansky)
    PSFSources.add_column(fluxerrujycol)

    # PSFSources['%sMAG_PSF' % band] = -2.5 * np.log10(PSFSources['FLUX_DENSITY']) + 23.9
    # PSFSources['e_%sMAG_PSF' % band] = 1.086 * (PSFSources['E_FLUX_DENSITY'] / PSFSources['FLUX_DENSITY'])

    # catalog for just clean sources (no flags)
    psfmag_clean = zero_psfmean + cleanPSFSources['MAG_POINTSOURCE']
    psfmagerr_clean = np.sqrt(cleanPSFSources['MAGERR_POINTSOURCE'] ** 2 + zero_psfstd ** 2)
    # ab mag conversion
    psfmag_clean = ab_convert(psfmag_clean, band=band, survey=survey)

    psfmagcol_clean = Column(psfmag_clean, name='%sMAG_PSF' % band, unit=u.ABmag)
    psfmagerrcol_clean = Column(psfmagerr_clean, name='e_%sMAG_PSF' % band, unit=u.ABmag)
    cleanPSFSources.add_column(psfmagcol_clean)
    cleanPSFSources.add_column(psfmagerrcol_clean)
    cleanPSFSources.remove_column('VIGNET')
    cleanPSFSources['FLUX_RADIUS'] = cleanPSFSources['FLUX_RADIUS'] * 0.498
    cleanPSFSources['FLUX_RADIUS'].unit = u.arcsec
    if 'FLUX_RADIUS_90' in cleanPSFSources.colnames:
        cleanPSFSources['FLUX_RADIUS_90'] = cleanPSFSources['FLUX_RADIUS_90'] * 0.498
        cleanPSFSources['FLUX_RADIUS_90'].unit = u.arcsec

    cleanpsfflux = psfmagcol_clean.to(u.microjansky)
    cleanpsffluxcol = Column(cleanpsfflux, name='FLUX_DENSITY', unit=u.microjansky)
    cleanPSFSources.add_column(cleanpsffluxcol)

    cleanfluxerr_ujy = cleanPSFSources['FLUXERR_POINTSOURCE'] * conv_factor
    cleanfluxerrujycol = Column(cleanfluxerr_ujy, name='E_FLUX_DENSITY', unit=u.microjansky)
    cleanPSFSources.add_column(cleanfluxerrujycol)

    cleanPSFSources['%sMAG_PSF' % band] = -2.5 * np.log10(cleanPSFSources['FLUX_DENSITY']) + 23.9
    cleanPSFSources['e_%sMAG_PSF' % band] = 1.086 * (cleanPSFSources['E_FLUX_DENSITY'] / cleanPSFSources['FLUX_DENSITY'])

    with fits.open(imageName, mode='update') as imagehdu:
        imagehdr = imagehdu[0].header
        # TODO if you replace the pix values, put the if statement back in
        # imagehdu[0].data = imagehdu[0].data * conv_factor  # adu * (uJy / adu) = uJy
        imagehdr.set('BUNIT', 'uJy', 'Physical units of the array values', after='EXTEND')
        imagehdr.set('CONV_FAC', conv_factor, 'uJy / ADU Conversion Factor, multiply img by this to get in uJy', after='BUNIT')
        print('Conversion of ADU to uJy calculated, med conversion factor: %.4f' % conv_factor)

        # print('BUNIT found already in header, skipping conversion.')
        imagehdu.close()

    PSFSources.write('%s.%s.ecsv' % (imageName, survey), overwrite=True)
    print('%s.%s.ecsv written, CSV w/ corrected mags' % (imageName, survey))

    # catalog conversion to AB

    print('Converting Vega surveys to AB to ensure correctness!')
    good_cat_stars[magcolname] = ab_convert(good_cat_stars[magcolname], band=band, survey=survey)
    ab_cat_stars = good_cat_stars

    return cleanPSFSources, PSFSources, psfweights_noclip, psf_clipped, ab_cat_stars


#%%


# New Source Search
def newsourcesearch(source_ra, source_dec, thresh, w, imageName, survey, band, Q, data, crop):
    # crossmatch for all detected sources
    print('Large grb radius inputted, using new source search to find'
          ' detected sources brighter than survey lim mag w/ no crossmatch')

    # sexigesimal conversion
    try:
        float(source_ra)
        source_ra = source_ra
        source_dec = source_dec
    except ValueError:
        coords = source_ra + ' ' + source_dec
        print('Sexagesimal RA = %s & Dec = %s' % (source_ra, source_dec))

        deci_coords = SkyCoord(coords, frame='icrs', unit=(u.hourangle, u.deg)).to_string()
        deci_coords = deci_coords.split(' ')

        source_ra = deci_coords[0]
        source_dec = deci_coords[1]

    # regenerate massCatCoords
    crop = int(crop)
    max_x = data.shape[0]
    max_y = data.shape[1]

    colnames = Q[0].colnames
    RA = colnames[0]
    DEC = colnames[1]

    mass_imCoords = w.all_world2pix(Q[0][RA], Q[0][DEC], 1)
    good_cat_stars = Q[0][np.where(
        (mass_imCoords[0] > crop) & (mass_imCoords[0] < (max_x - crop)) & (mass_imCoords[1] > crop) & (
                mass_imCoords[1] < (max_y - crop)))]

    massCatCoords = SkyCoord(ra=good_cat_stars[RA], dec=good_cat_stars[DEC], frame='icrs', unit='degree')
    print('Catalog cropped #:', len(massCatCoords))

    # initial crossmatch
    mag_ecsvname = '%s.%s.ecsv' % (imageName, survey)
    mag_ecsvtable = ascii.read(mag_ecsvname)
    mag_ecsvSources = mag_ecsvtable[(mag_ecsvtable['FLAGS'] == 0) & (mag_ecsvtable['FLAGS_MODEL'] == 0)]
    # print(len(mag_ecsvSources))
    mag_ecsvsourceCatCoords = SkyCoord(ra=mag_ecsvSources['ALPHA_J2000'], dec=mag_ecsvSources['DELTA_J2000'],
                                       frame='icrs',
                                       unit='degree')

    mag_ecsvsourceCatCoords = utils.pixel_to_skycoord(mag_ecsvSources['X_IMAGE'], mag_ecsvSources['Y_IMAGE'], w, origin=1)

    photoDistThresh = 1.0   # set higher to combat offset wcs in certain sources, maybe change for denser fields?
    idx_psfimage_noclean, idx_psfmass_noclean, d2d, d3d = massCatCoords.search_around_sky(mag_ecsvsourceCatCoords,
                                                                          photoDistThresh * u.arcsec)

    idx_psfimage_noclean_set = set(idx_psfimage_noclean)

    # rad / cross-match pruning
    PSFsources_nomatch = mag_ecsvSources[[i for i in range(len(mag_ecsvSources)) if i not in idx_psfimage_noclean_set]]   # removing previous crossmatched sources
    print('# of sources found after removing crossmatches: %i' % len(PSFsources_nomatch))


    sourcecoords = SkyCoord(ra=[source_ra], dec=[source_dec], frame='icrs', unit='degree')
    PSFsources_nomatchCatCoords = SkyCoord(ra=PSFsources_nomatch['ALPHA_J2000'], dec=PSFsources_nomatch['DELTA_J2000'],
                                       frame='icrs',
                                       unit='degree')
    idx_inputcoords, idx_PSFsources_nomatch, d2dd, d3dd = PSFsources_nomatchCatCoords.search_around_sky(sourcecoords,
                                                                                                      thresh * u.arcsec)
    PSFsources_nomatch = PSFsources_nomatch[idx_PSFsources_nomatch]     # implementing error radius

    # lim mag pruning
    for f in PHOTOMETRY_LIM_MAGS.keys():
        if survey == f:
            lim_mag = PHOTOMETRY_LIM_MAGS[f]
            break
    else:
        print('Error in finding lim mag for catalog, no matching catalog found?')
        sys.exit('No catalog match for lim mag.')

    PSFsources_new = PSFsources_nomatch[(PSFsources_nomatch['%sMAG_PSF' % band] < lim_mag) &
                                        (PSFsources_nomatch['%sMAG_PSF' % band] > 12.5)]
    print('# of sources found after all pruning: %i' % len(PSFsources_new))

    PSFsources_new.write('New_Sources.%s.%s.ecsv' % (imageName, survey), overwrite=True)
    print('New source catalog written!')


# %% optional GRB-specific photom


def GRB(ra, dec, imageName, survey, band, thresh, massCatCoords, good_cat_stars, directory, coordlist=None):
    mag_ecsvname = '%s.%s.ecsv' % (imageName, survey)
    mag_ecsvtable = ascii.read(mag_ecsvname)
    mag_ecsvcleanSources = mag_ecsvtable  # [(mag_ecsvtable['FLAGS'] == 0) & (mag_ecsvtable['FLAGS_MODEL'] == 0)]
    # mag_ecsvcleanSources['%sMAG_PSF' % band] = (
    #     ab_convert(mag_ecsvcleanSources['%sMAG_PSF' % band], band=band, survey=survey, revert=True))
    mag_ecsvsourceCatCoords = SkyCoord(ra=mag_ecsvcleanSources['ALPHA_J2000'], dec=mag_ecsvcleanSources['DELTA_J2000'],
                                       frame='icrs',
                                       unit='degree')

    # survey AB conversion for plot correctness
    # print('Converting Vega surveys to AB to ensure GRB reporting correctness!')
    # colnames = good_cat_stars.colnames
    # magcol = colnames[2]
    # good_cat_stars[magcol] = ab_convert(good_cat_stars[magcol], band=band, survey=survey)


    # postage stamp cutout fctn (png & fits)
    def grb_cutout(imageName, GRBcoords, photoDistThresh, loc=None, append=False, regprimename=None, regsurvname=None):
        imgdata = fits.getdata(imageName)
        img = fits.open(imageName)
        head = img[0].header
        w = WCS(head)

        size = 4 * photoDistThresh * u.arcsec
        cutout = Cutout2D(imgdata, GRBcoords, size, wcs=w, copy=True)
        region = CircleSkyRegion(center=GRBcoords[0], radius=Angle(thresh, unit='arcsec'))
        pix_region = region.to_pixel(cutout.wcs)

        mean, median, sigma_cut = sigma_clipped_stats(cutout.data)
        plt.figure(10, figsize=(8, 8))
        plt.imshow(cutout.data, vmin=median - 3 * sigma_cut, vmax=median + 3 * sigma_cut, origin='lower', cmap='viridis')
        pix_region.plot(color='cyan', ls='--', label='Input GRB threshold')

        if regprimename:
            primeregs = open('GRB_PRIME_srcs.reg', 'r')
            plt_primeregs = []
            primeallregs = [reg for reg in primeregs if reg != 'fk5\n']
            for reg in primeallregs:
                nums = re.findall(r'[-+]?\d*\.?\d+', reg)
                srcra = nums[0]
                srcdec = nums[1]
                srcrad = nums[2]
                srccoords = SkyCoord(ra=[srcra], dec=[srcdec], frame='icrs', unit='degree')
                srcreg = CircleSkyRegion(center=srccoords[0], radius=Angle(srcrad, unit='arcsec'))
                plt_primeregs.append(srcreg)

            for reg in plt_primeregs:
                pix_reg = reg.to_pixel(cutout.wcs)
                pix_reg.plot(color='red', ls='-', label='PRIME Source')

        if regsurvname:
            survregs = open('GRB_%s_srcs.reg' % survey, 'r')
            plt_survregs = []
            survallregs = [reg for reg in survregs if reg != 'fk5\n']
            for reg in survallregs:
                nums = re.findall(r'[-+]?\d*\.?\d+', reg)
                srcra = nums[0]
                srcdec = nums[1]
                srcrad = nums[2]
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

        if loc:
            savename = 'GRB_%s_Cutout_%s_loc_%d' % (band, survey, loc)
            threshname = 'GRB_queries_thresh.reg'
        else:
            savename = 'GRB_%s_Cutout_%s' % (band, survey)
            threshname = 'GRB_query_thresh.reg'
        plt.savefig(savename + '.png', dpi=300)
        plt.clf()

        fits.writeto(savename + '.fits', cutout.data, cutout.wcs.to_header(), overwrite=True)

        # ds9 regions
        if append:
            newtext = open(threshname, 'a')  # input threshold
            newtext.write(f'\ncircle({ra}, {dec}, {photoDistThresh}") # color=cyan width=2 text={{Query Thresh}}')
        else:
            newtext = open(threshname, 'w+')       # input threshold
            newtext.write('fk5')
            newtext.write(f'\ncircle({ra}, {dec}, {photoDistThresh}") # color=cyan width=2 text={{Query Thresh}}')

        print(' Exported cutout & ds9 region of GRB search area!')

        return savename, threshname

    # source ds9 region writing
    def source_reg_gen(src_ra=0, src_dec=0, rad=2, src_survey=None, append=False):
        if not src_survey:
            name = 'GRB_PRIME_srcs.reg'
            color = 'red'
        else:
            name = 'GRB_%s_srcs.reg' % survey
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
    def mag_diff_calc(survey_cat, prime_cat, survey_idx, prime_idx, band):
        colnames = survey_cat.colnames
        magcolname = colnames[2]
        magerrcolname = colnames[3]

        if len(prime_idx) > 1:
            mag_diff = float(survey_cat[magcolname][survey_idx]) - float(prime_cat['%sMAG_PSF' % band][prime_idx])

            # errors in quad
            comb_err = np.sqrt(survey_cat[magerrcolname][survey_idx] ** 2 +
                               prime_cat['e_%sMAG_PSF' % band][prime_idx] ** 2)

        else:
            # survey mag - prime mag
            mag_diff = float(survey_cat[magcolname][survey_idx]) - float(prime_cat['%sMAG_PSF' % band])

            # errors in quad
            comb_err = np.sqrt(survey_cat[magerrcolname][survey_idx] ** 2 +
                               prime_cat['e_%sMAG_PSF' % band] ** 2)

        err_3_sig = 3 * comb_err

        if abs(mag_diff) > err_3_sig:
            flg = 2
            # print(' At least 1 source in threshold has a significant mag difference to the catalog!')
        else:
            flg = 1
        return mag_diff, np.int16(flg)

    # generation of html file
    def html_gen(data, directory, savename, threshname, band, survey, ra, dec, thresh, survname=None, primename=None):
        df = data.to_pandas()
        tbl_html = df.to_html(index=False, classes="my-table")

        with open(savename+'.png', "rb") as img_file:
            encoded = base64.b64encode(img_file.read()).decode("utf-8")

        fits_items = [savename+'.fits', threshname, survname, primename]
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

        with open('GRB_Multisource_%s_Data_%s.html' % (band, survey), "w") as f:
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

                grb_mag = mag_ecsvcleanSources[idx_GRBcleanpsf]['%sMAG_PSF' % band][0]
                grb_magerr = mag_ecsvcleanSources[idx_GRBcleanpsf]['e_%sMAG_PSF' % band][0]

                grb_ra = mag_ecsvcleanSources[idx_GRBcleanpsf]['ALPHA_J2000'][0]
                grb_dec = mag_ecsvcleanSources[idx_GRBcleanpsf]['DELTA_J2000'][0]
                grb_rad = mag_ecsvcleanSources[idx_GRBcleanpsf]['FLUX_RADIUS'][0]
                grb_snr = mag_ecsvcleanSources[idx_GRBcleanpsf]['SNR_WIN'][0]
                grb_dist = d2d[0].to(u.arcsec)
                grb_dist = grb_dist / u.arcsec

                print(' Detected GRB ra = %.6f, dec = %.6f, with 50 percent flux radius (HWHM) = %.3f arcsec and SNR = %.3f' % (
                    grb_ra, grb_dec, grb_rad, grb_snr))
                print(' %s magnitude of GRB is %.2f +/- %.2f' % (band, grb_mag, grb_magerr))

                # survey crsmtch check
                idx_both, idx_bothcleanpsf, d2d_crs, d3d_crs = massCatCoords.search_around_sky(
                    mag_ecsvsourceCatCoords[idx_GRBcleanpsf],
                    grb_rad * u.arcsec)
                if len(idx_bothcleanpsf) > 0:
                    # print(' Detected source crossmatched to existing %s source!' % survey)
                    prime_crs_cat = mag_ecsvcleanSources[idx_GRBcleanpsf]
                    mag_diff_crs, survey_flg = mag_diff_calc(survey_cat=good_cat_stars, prime_cat=prime_crs_cat,
                                                             survey_idx=idx_bothcleanpsf, prime_idx=[idx_GRBcleanpsf][idx_both],
                                                             band=band)
                else:
                    mag_diff_crs = 99
                    survey_flg = np.int16(0)

                grbdata = Table()
                grbdata['RA'] = np.round(np.array([grb_ra]), decimals=5) * u.deg
                grbdata['DEC'] = np.round(np.array([grb_dec]), decimals=5) * u.deg
                grbdata['%sMag' % band] = np.round(np.array([grb_mag]), decimals=3) * u.ABmag
                grbdata['%sMag_Err' % band] = np.round(np.array([grb_magerr]), decimals=3) * u.ABmag
                grbdata['%sMag_Err_Crsmtch' % band] = np.round(np.array([mag_diff_crs]), decimals=3) * u.ABmag
                grbdata['Radius'] = np.round(np.array([grb_rad]), decimals=2) * u.arcsec
                grbdata['SNR'] = np.round(np.array([grb_snr]), decimals=2)
                grbdata['Distance'] = np.round(np.array([grb_dist]), decimals=5) * u.arcsec
                grbdata['Survey Crsmtch'] = survey_flg


                grbdata.write('GRB_%s_Data_%s_loc_%d.ecsv' % (band, survey, key), overwrite=True)

                source_reg_gen(grb_ra, grb_dec, rad=grb_rad, append=True)
                print(' Generated GRB data table & source DS9 regions!')
            elif len(idx_GRBcleanpsf) > 1:
                print(' Multiple sources detected in search radius (ra = %s, dec = %s, rad = %s arcsec)'
                      ', refer to .ecsv file for source info!' % (ra, dec, photoDistThresh))
                mag_ar = []
                mag_err_ar = []
                ra_ar = []
                dec_ar = []
                rad_ar = []
                snr_ar = []
                dist_ar = []
                diff_ar = []
                crsmtch_ar = []
                idx_GRBcleanpsflist = idx_GRBcleanpsf.tolist()
                for i in idx_GRBcleanpsflist:
                    grb_mag = mag_ecsvcleanSources[idx_GRBcleanpsf]['%sMAG_PSF' % band][idx_GRBcleanpsflist.index(i)]
                    mag_ar.append(grb_mag)
                    grb_magerr = mag_ecsvcleanSources[idx_GRBcleanpsf]['e_%sMAG_PSF' % band][
                        idx_GRBcleanpsflist.index(i)]
                    mag_err_ar.append(grb_magerr)
                    grb_ra = mag_ecsvcleanSources[idx_GRBcleanpsf]['ALPHA_J2000'][idx_GRBcleanpsflist.index(i)]
                    ra_ar.append(grb_ra)
                    grb_dec = mag_ecsvcleanSources[idx_GRBcleanpsf]['DELTA_J2000'][idx_GRBcleanpsflist.index(i)]
                    dec_ar.append(grb_dec)
                    grb_rad = mag_ecsvcleanSources[idx_GRBcleanpsf]['FLUX_RADIUS'][idx_GRBcleanpsflist.index(i)]
                    rad_ar.append(grb_rad)
                    grb_snr = mag_ecsvcleanSources[idx_GRBcleanpsf]['SNR_WIN'][idx_GRBcleanpsflist.index(i)]
                    snr_ar.append(grb_snr)
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
                        mag_diff_crs, survey_flg = mag_diff_calc(survey_cat=good_cat_stars,
                                                                 prime_cat=prime_crs_cat,
                                                                 survey_idx=idx_bothcleanpsf, prime_idx=[i][idx_both],
                                                                 band=band)
                    else:
                        mag_diff_crs = 99
                        survey_flg = np.int16(0)

                    diff_ar.append(mag_diff_crs)
                    crsmtch_ar.append(survey_flg)

                    source_reg_gen(grb_ra, grb_dec, rad=grb_rad, append=True)

                grbdata = Table()
                grbdata['RA'] = np.round(np.array(ra_ar), decimals=5) * u.deg
                grbdata['DEC'] = np.round(np.array(dec_ar), decimals=5) * u.deg
                grbdata['%sMag' % band] = np.round(np.array(mag_ar), decimals=3) * u.ABmag
                grbdata['%sMag_Err' % band] = np.round(np.array(mag_err_ar), decimals=3) * u.ABmag
                grbdata['%sMag_Err_Crsmtch' % band] = np.round(np.array(diff_ar), decimals=3) * u.ABmag
                grbdata['Radius'] = np.round(np.array(rad_ar), decimals=2) * u.arcsec
                grbdata['SNR'] = np.round(np.array(snr_ar), decimals=2)
                grbdata['Distance'] = np.round(np.array(dist_ar), decimals=5) * u.arcsec
                grbdata['Survey_Crsmtch'] = np.array(crsmtch_ar)
                grbdata.write('GRB_Multisource_%s_Data_%s_loc_%d.ecsv' % (band, survey, key), overwrite=True)
                print(' Generated GRB data table & source DS9 regions!')
            else:
                print(' GRB source at inputted coords %s and %s not found in PRIME catalog, perhaps increase photoDistThresh or '
                      'alter sextractor params?' % (ra, dec))
    else:
        print(' idx size = %d' % len(idx_GRBcleanpsf))
        if len(idx_GRBcleanpsf) == 1:
            print(' GRB source at inputted coords %s and %s, rad = %s arcsec found!' % (ra, dec, photoDistThresh))

            grb_mag = mag_ecsvcleanSources[idx_GRBcleanpsf]['%sMAG_PSF' % band][0]
            grb_magerr = mag_ecsvcleanSources[idx_GRBcleanpsf]['e_%sMAG_PSF' % band][0]

            grb_ra = mag_ecsvcleanSources[idx_GRBcleanpsf]['ALPHA_J2000'][0]
            grb_dec = mag_ecsvcleanSources[idx_GRBcleanpsf]['DELTA_J2000'][0]
            grb_rad = mag_ecsvcleanSources[idx_GRBcleanpsf]['FLUX_RADIUS'][0]
            grb_snr = mag_ecsvcleanSources[idx_GRBcleanpsf]['SNR_WIN'][0]
            grb_dist = d2d[0].to(u.arcsec)
            grb_dist = grb_dist / u.arcsec

            print(' Detected GRB ra = %.6f, dec = %.6f, with 50 percent flux radius (HWHM) = %.3f arcsec and SNR = %.3f' % (
                grb_ra, grb_dec, grb_rad, grb_snr))
            print(' %s magnitude of GRB is %.2f +/- %.2f' % (band, grb_mag, grb_magerr))

            # survey crsmtch check
            idx_both, idx_bothcleanpsf, d2d_crs, d3d_crs = (
                massCatCoords.search_around_sky(mag_ecsvsourceCatCoords[idx_GRBcleanpsf], grb_rad * u.arcsec))
            if len(idx_bothcleanpsf) > 0:
                # print(' Detected source crossmatched to existing %s source!' % survey)
                prime_crs_cat = mag_ecsvcleanSources[idx_GRBcleanpsf]
                mag_diff_crs, survey_flg = mag_diff_calc(survey_cat=good_cat_stars, prime_cat=prime_crs_cat,
                                                         survey_idx=idx_bothcleanpsf, prime_idx=[idx_GRBcleanpsf][idx_both],
                                                         band=band)
            else:
                mag_diff_crs = 99
                survey_flg = np.int16(0)

            grbdata = Table()
            grbdata['RA'] = np.round(np.array([grb_ra]), decimals=5) * u.deg
            grbdata['DEC'] = np.round(np.array([grb_dec]), decimals=5) * u.deg
            grbdata['%sMag' % band] = np.round(np.array([grb_mag]), decimals=3) * u.ABmag
            grbdata['%sMag_Err' % band] = np.round(np.array([grb_magerr]), decimals=3) * u.ABmag
            grbdata['%sMag_Err_Crsmtch' % band] = np.round(np.array([mag_diff_crs]), decimals=3) * u.ABmag
            grbdata['Radius'] = np.round(np.array([grb_rad]), decimals=2) * u.arcsec
            grbdata['SNR'] = np.round(np.array([grb_snr]), decimals=2)
            grbdata['Distance'] = np.round(np.array([grb_dist]), decimals=5) * u.arcsec
            grbdata['Survey Crsmtch'] = survey_flg

            grbdata.write('GRB_%s_Data_%s.ecsv' % (band, survey), overwrite=True)

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
            ra_ar = []
            dec_ar = []
            rad_ar = []
            snr_ar = []
            dist_ar = []
            diff_ar = []
            crsmtch_ar = []
            idx_GRBcleanpsflist = idx_GRBcleanpsf.tolist()
            source_reg_gen()
            for i in idx_GRBcleanpsflist:
                grb_mag = mag_ecsvcleanSources[idx_GRBcleanpsf]['%sMAG_PSF' % band][idx_GRBcleanpsflist.index(i)]
                mag_ar.append(grb_mag)
                grb_magerr = mag_ecsvcleanSources[idx_GRBcleanpsf]['e_%sMAG_PSF' % band][idx_GRBcleanpsflist.index(i)]
                mag_err_ar.append(grb_magerr)
                grb_ra = mag_ecsvcleanSources[idx_GRBcleanpsf]['ALPHA_J2000'][idx_GRBcleanpsflist.index(i)]
                ra_ar.append(grb_ra)
                grb_dec = mag_ecsvcleanSources[idx_GRBcleanpsf]['DELTA_J2000'][idx_GRBcleanpsflist.index(i)]
                dec_ar.append(grb_dec)
                grb_rad = mag_ecsvcleanSources[idx_GRBcleanpsf]['FLUX_RADIUS'][idx_GRBcleanpsflist.index(i)]
                rad_ar.append(grb_rad)
                grb_snr = mag_ecsvcleanSources[idx_GRBcleanpsf]['SNR_WIN'][idx_GRBcleanpsflist.index(i)]
                snr_ar.append(grb_snr)
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
                    mag_diff_crs, survey_flg = mag_diff_calc(survey_cat=good_cat_stars, prime_cat=prime_crs_cat,
                                                             survey_idx=idx_bothcleanpsf, prime_idx=idx_both, band=band)
                else:
                    mag_diff_crs = 99
                    survey_flg = np.int16(0)

                diff_ar.append(mag_diff_crs)
                crsmtch_ar.append(survey_flg)

                regprimename = source_reg_gen(grb_ra, grb_dec, rad=grb_rad, append=True)

            grbdata = Table()
            grbdata['RA'] = np.round(np.array(ra_ar), decimals=5) * u.deg
            grbdata['DEC'] = np.round(np.array(dec_ar), decimals=5) * u.deg
            grbdata['%sMag' % band] = np.round(np.array(mag_ar), decimals=3) * u.ABmag
            grbdata['%sMag_Err' % band] = np.round(np.array(mag_err_ar), decimals=3) * u.ABmag
            grbdata['%sMag_Err_Crsmtch' % band] = np.round(np.array(diff_ar), decimals=3) * u.ABmag
            grbdata['Radius'] = np.round(np.array(rad_ar), decimals=2) * u.arcsec
            grbdata['SNR'] = np.round(np.array(snr_ar), decimals=2)
            grbdata['Distance'] = np.round(np.array(dist_ar), decimals=5) * u.arcsec
            grbdata['Survey_Crsmtch'] = np.array(crsmtch_ar)
            # grbdata.write('GRB_Multisource_%s_Data_%s.ecsv' % (band, survey), overwrite=True)
            grbdata.write('GRB_Multisource_%s_Data_%s.ecsv' % (band, survey), overwrite=True)

            savename, threshname = grb_cutout(imageName, GRBcoords, photoDistThresh,
                                              regprimename=regprimename, regsurvname=regsurvname)

            html_gen(grbdata, directory, savename, threshname, band, survey, ra, dec, thresh, regsurvname, regprimename)

            print(' Generated GRB data table & source DS9 regions!')
        else:
            print(' GRB source at inputted coords %s and %s not found in PRIME catalog, perhaps increase photoDistThresh or '
                  'alter sextractor params?' % (ra, dec))


# %% optional plots


def photometry_plots(cleanPSFsources, PSFsources, data, imageName, survey, band, good_cat_stars, idx_psfmass, idx_psfimage,
                     psfweights_noclip, psf_clipped, sigma):

    # appropriate mag column
    colnames = good_cat_stars.colnames
    magcol = colnames[2]
    magerrcol = colnames[3]

    # survey AB conversion for plot correctness
    # print('Converting Vega surveys to AB to ensure plot correctness!')
    # good_cat_stars[magcol] = ab_convert(good_cat_stars[magcol], band=band, survey=survey)

    chip = imageName[-6]
    if len(imageName) <= 16:
        num = 'img'
    else:
        num = imageName[-16:-8]

    def predict_y_for(x, m, b):
        return m * x + b

    plt.close('all')
    # mag comparison plot
    plt.figure(1, figsize=(8, 8))
    plt.plot(cleanPSFsources['%sMAG_PSF' % band][idx_psfimage], good_cat_stars['%s' % magcol][idx_psfmass],
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
    x_flx_lin = cleanPSFsources['FLUX_DENSITY'][idx_psfimage][~psf_clipped.mask]
    x_flx = np.log(x_flx_lin)
    y_flx = good_cat_stars['%s' % magcol][idx_psfmass][~psf_clipped.mask]
    x_const_flx = sm.add_constant(x_flx)
    model_flx = sm.WLS(y_flx, x_const_flx, weights=psfweights_noclip[~psf_clipped.mask]).fit()
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
    plt.plot(cleanPSFsources['FLUX_DENSITY'][idx_psfimage][~psf_clipped.mask],
             good_cat_stars['%s' % magcol][idx_psfmass][~psf_clipped.mask],
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
    plt.savefig('%s_C%s_flux_mag_plot_sig_%s.png' % (survey, chip, num))
    plt.clf()
    print('Saved flux v. mag plot to dir!')

    # residual fits
    x = cleanPSFsources['%sMAG_PSF' % band][idx_psfimage]
    y = good_cat_stars['%s' % magcol][idx_psfmass]
    x_const = sm.add_constant(x)
    # model = sm.WLS(y, x, weights=psfweights).fit()
    model2 = sm.WLS(y, x_const, weights=psfweights_noclip).fit()
    # m = model.params[0]
    # merr = model.bse[0]
    m2 = model2.params[1]
    m2err = model2.bse[1]
    b2 = model2.params[0]
    b2err = model2.bse[0]
    # rsquare = model.rsquared
    # rss = model.ssr
    # res_err = np.std(model.resid)
    rsquare2 = model2.rsquared
    rss2 = model2.ssr

    avg2 = np.average(model2.resid, weights=psfweights_noclip)
    var2 = np.average((model2.resid - avg2)**2, weights=psfweights_noclip)
    res2_err = np.sqrt(var2)
    # res2_err = np.std(model2.resid)

    # residual fit - 3 sigma clip
    x_sig = good_cat_stars['%s' % magcol][idx_psfmass][~psf_clipped.mask]
    y_sig = cleanPSFsources['%sMAG_PSF' % band][idx_psfimage][~psf_clipped.mask]
    x_const_sig = sm.add_constant(x_sig)
    model_sig = sm.WLS(y_sig, x_const_sig, weights=psfweights_noclip[~psf_clipped.mask]).fit()
    m_sig = model_sig.params[1]
    m_sigerr = model_sig.bse[1]
    b_sig = model_sig.params[0]
    b_sigerr = model_sig.bse[0]
    rsquare_sig = model_sig.rsquared
    rss_sig = model_sig.ssr
    model_sig_resid = model_sig.resid

    print('# of crossmatched sources used in %s sig fit: %i'
          % (sigma, len(cleanPSFsources['%sMAG_PSF' % band][idx_psfimage][~psf_clipped.mask])))
    print('%s sig fit: slope = %.4f +/- %.4f' % (sigma, m_sig, m_sigerr))
    print('%s sig fit: y-int = %.4f +/- %.4f' % (sigma, b_sig, b_sigerr))

    # sigma-clipped data
    # x_sig = cleanPSFsources['%sMAG_PSF' % band][idx_psfimage][~psf_clipped.mask]
    # y_sig = good_cat_stars['%s' % magcol][idx_psfmass][~psf_clipped.mask]
    #
    # # errors
    # x_err = cleanPSFsources['e_%sMAG_PSF' % band][idx_psfimage][~psf_clipped.mask]
    # y_err = good_cat_stars['%s' % magerrcol][idx_psfmass][~psf_clipped.mask]
    #
    # # define linear function for ODR
    # def linear_func(B, x):
    #     return B[0] * x + B[1]  # B[0]=slope, B[1]=intercept
    #
    # linear_model = odr.Model(linear_func)
    #
    # # create RealData object with x, y and their errors
    # odr_data = odr.RealData(x_sig, y_sig, sx=x_err, sy=y_err)
    #
    # # set up ODR with initial guess [slope, intercept] = [1, 0]
    # odr_obj = odr.ODR(odr_data, linear_model, beta0=[1.0, 0.0])
    #
    # # run ODR fit
    # odr_out = odr_obj.run()
    #
    # m_sig, b_sig = odr_out.beta
    # m_sigerr, b_sigerr = odr_out.sd_beta
    #
    # # residuals
    # model_sig_resid = y_sig - (m_sig * x_sig + b_sig)
    # model_sig_resid = np.array(model_sig_resid, dtype=float)
    #
    # # RSS (residual sum of squares)
    # rss_sig = np.sum(model_sig_resid ** 2)
    #
    # # R squared
    # y_mean = np.mean(y_sig)
    # ss_tot = np.sum((y_sig - y_mean) ** 2)
    # rsquare_sig = 1 - rss_sig / ss_tot
    #
    # print('# of crossmatched sources used in ODR fit:', len(x_sig))
    # print('ODR fit: slope = %.4f +/- %.4f' % (m_sig, m_sigerr))
    # print('ODR fit: y-int = %.4f +/- %.4f' % (b_sig, b_sigerr))

    # residual fit 3 sig clip - bin errors and stats
    mags = range(10,22)
    res_errs = []
    res_means = []
    res_ranges = []
    res_skews = []
    res_nums = []
    x_arr = np.arange(10+0.5, 22+0.5, 1)
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
    bintable['Bin Range (mag)'] = res_ranges
    bintable['Mean (mag)'] = res_means
    bintable['Spread (mag)'] = res_errs
    bintable['Skew'] = res_skews
    bintable['Source Number'] = res_nums

    bintable.write('Resid_%s-sig_Data_%s_C%s_%s.ecsv' % (sigma, band, chip, survey), overwrite=True)
    print('Residual data table for %s sigma clip written!' % sigma)

    # t = Table()
    # t['Residuals'] = model_sig.resid
    # t['mags'] = x_sig
    # t['Weights'] = psfweights_noclip[~psf_clipped.mask]
    # t.write('residual_3sig.ecsv', overwrite=True)

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

    # res plot, y int include
    plt.figure(2, figsize=(8, 6))
    plt.scatter(good_cat_stars['%s' % magcol][idx_psfmass][~psf_clipped.mask], model_sig_resid, color='red')
    plt.ylim(-1, 1)
    plt.xlim(10, 21)
    plt.title('%s Residuals - %s Sigma Clip' % (survey, sigma))
    plt.ylabel('Residuals')
    plt.xlabel('%s %s Mags' % (survey, band))
    # plt.axhline(y=ressig_err, color='blue', linestyle='--', linewidth=1)
    # plt.axhline(y=ressig_err * 2, color='green', linestyle='--', linewidth=1)
    # plt.axhline(y=-ressig_err, color='blue', linestyle='--', linewidth=1)
    # plt.axhline(y=-ressig_err * 2, color='green', linestyle='--', linewidth=1)
    plt.scatter(x_arr, res_errs, marker='_', s=1625, c='blue')
    plt.scatter(x_arr, -res_errs, marker='_', s=1625, c='blue')
    plt.axhline(y=0, color='black', linestyle='--', linewidth=1)
    plt.legend(['Residuals',r'1 $\sigma$ range = [%.3f - %.3f]' % (res_errs_min, res_errs_max)], loc='lower left',
               markerscale=0.5)
    info2 = ('eqn: y = mx+b' + '\nslope = %.4f +/- %.4f' % (m_sig, m_sigerr)) + (
                '\nintercept = %.3f +/- %.3f' % (b_sig, b_sigerr)) + ('\nR$^{2}$ = %.3f' % rsquare_sig) + ('\nRSS = %d' % rss_sig)
    plt.text(15, -0.9, info2, fontsize=9, bbox=dict(facecolor='white', edgecolor='black', pad=5.0))
    # plt.savefig('%s_C%s_residual_plot_int_%s.png' % (survey, chip, num), dpi=300)
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

    plt.figure(3, figsize=(10, 6))
    plt.hist2d(x=good_cat_stars['%s' % magcol][idx_psfmass][~psf_clipped.mask], y=model_sig_resid,
               bins=[bin_num_int, bin_num_int], range=[[10, 21],[-1, 1]], cmap='gist_heat_r')
    plt.colorbar(label='Density')
    plt.ylim(-1, 1)
    plt.xlim(10, 21)
    plt.title('%s Residuals - %s Sigma Clip - Density Histogram' % (survey, sigma))
    plt.ylabel('Residuals')
    plt.xlabel('%s %s Mags' % (survey, band))
    # plt.axhline(y=ressig_err, color='blue', linestyle='--', linewidth=1)
    # plt.axhline(y=ressig_err * 2, color='green', linestyle='--', linewidth=1)
    # plt.axhline(y=-ressig_err, color='blue', linestyle='--', linewidth=1)
    # plt.axhline(y=-ressig_err * 2, color='green', linestyle='--', linewidth=1)
    plt.scatter(x_arr, res_errs, marker='_', s=1625, c='blue')
    plt.scatter(x_arr, -res_errs, marker='_', s=1625, c='blue')
    plt.axhline(y=0, color='black', linestyle='--', linewidth=1)
    plt.legend([r'1 $\sigma$ range = [%.3f - %.3f]' % (res_errs_min, res_errs_max)], loc='lower left',
               markerscale=0.5)
    infohist = ('eqn: y = mx+b' + '\nslope = %.4f +/- %.4f' % (m_sig, m_sigerr)) + (
                '\nintercept = %.3f +/- %.3f' % (b_sig, b_sigerr)) + ('\nR$^{2}$ = %.3f' % rsquare_sig) + ('\nRSS = %d' % rss_sig)
    plt.text(15, -0.9, infohist, fontsize=9, bbox=dict(facecolor='white', edgecolor='black', pad=5.0))
    plt.savefig('%s_C%s_residual_plot_int_hist_%s.png' % (survey, chip, num), dpi=300)
    plt.clf()

    print('Saved y-int residual plots to dir!')

    # WLS fit line over data plot

    txt = ('slope = %.4f' % m2 + '\nslope err = %.4f' % m2err + '\nint = %.4f' % b2 + '\nint err = %.4f' % b2err)

    plt.figure(4, figsize=(8, 8))
    plt.xlim(10, 22)
    plt.ylim(10, 22)
    plt.title('%s vs PRIME w/ Weighted Fit' % survey)
    plt.grid()
    plt.ylabel('PRIME %s Mags' % band, fontsize=15)
    plt.xlabel('%s %s Mags' % (survey, band), fontsize=15)
    plt.scatter(good_cat_stars['%s' % magcol][idx_psfmass], cleanPSFsources['%sMAG_PSF' % band][idx_psfimage])
    plt.plot(good_cat_stars['%s' % magcol][idx_psfmass],
             predict_y_for(good_cat_stars['%s' % magcol][idx_psfmass], m2, b2), c='r')
    box = dict(facecolor='white')
    plt.text(11, 18, txt, fontsize=12, bbox=box)
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
    plt.hist2d(x=good_cat_stars['%s' % magcol][idx_psfmass], y=cleanPSFsources['%sMAG_PSF' % band][idx_psfimage],
               bins=[bin_num, bin_num], range=[[10, 22],[10, 22]], cmap='gist_heat_r')
    plt.plot(good_cat_stars['%s' % magcol][idx_psfmass],
             predict_y_for(good_cat_stars['%s' % magcol][idx_psfmass], m2, b2), c='b')
    plt.colorbar(label='Density')
    box = dict(facecolor='white')
    plt.text(11, 18, txt, fontsize=12, bbox=box)
    plt.savefig('%s_C%s_WLS_fit_hist_plot_%s.png' % (survey, chip, num), dpi=300)
    plt.clf()

    print('Saved WLS fit plots to dir!')

    # WLS fit for 3 sig clip of data
    sigtxt = ('slope = %.4f' % m_sig + '\nslope err = %.4f' % m_sigerr + '\nint = %.4f' % b_sig +
              '\nint err = %.4f' % b_sigerr)

    plt.figure(6, figsize=(8, 8))
    plt.xlim(10, 22)
    plt.ylim(10, 22)
    plt.title('%s vs PRIME w/ Weighted Fit - %s Sigma Clip' % (survey, sigma))
    plt.grid()
    plt.ylabel('PRIME %s Mags' % band, fontsize=15)
    plt.xlabel('%s %s Mags' % (survey, band), fontsize=15)
    plt.scatter(good_cat_stars['%s' % magcol][idx_psfmass][~psf_clipped.mask],
                cleanPSFsources['%sMAG_PSF' % band][idx_psfimage][~psf_clipped.mask])
    plt.plot(good_cat_stars['%s' % magcol][idx_psfmass][~psf_clipped.mask],
             predict_y_for(good_cat_stars['%s' % magcol][idx_psfmass][~psf_clipped.mask], m_sig, b_sig), c='r')
    box = dict(facecolor='white')
    plt.text(11, 18, sigtxt, fontsize=12, bbox=box)
    # plt.savefig('%s_C%s_WLS_fit_3sig_plot_%s.png' % (survey, chip, num), dpi=300)

    # WLS 3 sig hist density plot
    plt.figure(7, figsize=(10, 8))
    plt.clf()
    plt.xlim(10, 22)
    plt.ylim(10, 22)
    plt.title('%s vs PRIME w/ Weighted Fit - %s Sigma Clip - Density Histogram' % (survey, sigma))
    plt.grid()
    plt.ylabel('PRIME %s Mags' % band, fontsize=15)
    plt.xlabel('%s %s Mags' % (survey, band), fontsize=15)
    plt.hist2d(x=good_cat_stars['%s' % magcol][idx_psfmass][~psf_clipped.mask],
               y= cleanPSFsources['%sMAG_PSF' % band][idx_psfimage][~psf_clipped.mask],
               bins=[bin_num, bin_num], range=[[10, 22],[10, 22]], cmap='gist_heat_r')
    plt.plot(good_cat_stars['%s' % magcol][idx_psfmass][~psf_clipped.mask],
             predict_y_for(good_cat_stars['%s' % magcol][idx_psfmass][~psf_clipped.mask], m_sig, b_sig), c='b')
    plt.colorbar(label='Density')
    box = dict(facecolor='white')
    plt.text(11, 18, sigtxt, fontsize=12, bbox=box)
    plt.savefig('%s_C%s_WLS_fit_3sig_hist_plot_%s.png' % (survey, chip, num), dpi=300)
    plt.clf()

    print('Saved WLS 3 sig fit plots to dir!')

    # flux vs mag plot - histogram vers.
    plt.figure(10, figsize=(8, 8))
    plt.hist2d(x=cleanPSFsources['FLUX_DENSITY'][idx_psfimage], y=good_cat_stars['%s' % magcol][idx_psfmass],
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
    # binning
    all_mags_all = PSFsources[PSFsources['%sMAG_PSF' % band] < 25]
    all_mags = all_mags_all['%sMAG_PSF' % band]

    bin_vals = np.array(np.arange(12, 25.5, 0.1))
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
    idxmax = all_sources.index(max(all_sources))

    split_sources = all_sources[idxmax:]

    halfmax = max(all_sources) / 2
    halfmaxpt = list(min(enumerate(split_sources), key=lambda x: abs(halfmax - x[1])))
    halfmaxpt = [halfmaxpt[0] + idxmax, halfmaxpt[1]]

    limmag = round(bin_vals[halfmaxpt[0]], 1)

    print('Lim Mag = ', limmag)

    plt.figure(8, figsize=(24, 8))
    plt.bar(bin_vals, height=all_sources, width=0.1, align='edge', color='red', edgecolor='black')
    plt.axhline(halfmax, linestyle='--')
    plt.axvline(limmag, color='b', linewidth=2)
    # plt.bar(bin_vals,height=all_osaka_sources,width=0.1,align='edge',color='blue',edgecolor='black',alpha=0.5)
    # plt.bar(bin_vals-0.5,height=avginvsnr,width=0.5,align='edge',color='blue',edgecolor='black')
    # plt.axvline(x=bin_vals[imin+1],linestyle='--',linewidth=2)
    # plt.axvline(x=bin_vals[iimin+1],linestyle='--',linewidth=2,color='green')
    # plt.axhline(y=0.3,color='black', linestyle='-.', linewidth=1)
    # plt.axhline(y=0.2,color='black', linestyle='-.', linewidth=1)
    # bin_ticks = np.arange(12,21.25,0.5)
    xticks = np.arange(12, 25.5, 0.5)
    plt.xticks(xticks, fontsize=10)
    plt.grid()
    # if len(PSFsources) < 3500:
    #     plt.ylim(0, 250)
    # elif len(PSFsources) > 50000:
    #     plt.ylim(0, 5000)
    # else:
    #     plt.ylim(0, 2000)
    # plt.yscale('log')
    # plt.yticks([0,0.05,0.1,0.15,0.2,0.25])
    # plt.legend([r'5 $\sigma$ Limit'],loc='upper left',fontsize=12)
    # plt.legend([f'Bin of 0.2 Error = {bin_vals[imin+1]}',f'Bin of 0.3 Error = {bin_vals[iimin+1]}'])
    # plt.legend([f'Bin of 0.2 Error (>16 mag) = {bin_vals[imin+1]}'])
    # plt.legend(['Our Field','Osaka Field'],loc='upper right')
    plt.yscale('log')
    plt.title('PRIME Limiting Mag Plot')
    plt.ylabel('Number of Sources')
    plt.xlabel('%s Magnitude' % band)
    plt.legend(['Half Max = %s' % round(halfmax, 1), 'Limiting Mag = %s' % round(limmag, 1)], fontsize=15)
    plt.savefig('%s_C%s_lim_mag_plot_%s.png' % (survey, chip, num), dpi=300)
    print('Saved lim mag plot to dir!')
    plt.clf()

    # Crossmatch location check plot
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

    plt.savefig('%s_C%s_source_check_plot_%s.png' % (survey, chip, num), dpi=200)
    print('Saved source location check plot to dir!')
    plt.clf()

    plt.close('all')

    print('Writing relevant plot info to image header...')
    with fits.open(imageName, mode='update') as hdul:
        hdr = hdul[0].header
        hdr.set('fit_m', m_sig, 'WLS %s sig fit slope' % sigma, after='Survey')
        hdr.set('e_fit_m', m_sigerr, 'Error in WLS %s sig fit slope' % sigma, after='fit_m')
        hdr.set('fit_b', b_sig, 'WLS %s sig fit intercept' % sigma, after='e_fit_m')
        hdr.set('e_fit_b', b_sigerr, 'Error in WLS %s sig fit intercept' % sigma, after='fit_b')
        hdr.set('Lim_Mag', limmag, 'Source Histogram FWHM Limiting Mag', after='e_fit_b')
        hdul.close()

    return m_sig, b_sig

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
        sys.exit('Use supported units.')
    return arcconvert


#%% automated y int fit calibration


def int_calibration(
        name, directory, band, crop, sigma, given_catalog, survey,
        mag_low_lim, mag_high_lim, grb_ra, grb_dec,
        grb_coordlist, grb_radius, max_int
):
    print('3 sigma fit y-intercept > %s! Redoing photometry w/ sigma = %s, mag low cutoff = %s\n' % (max_int, sigma, mag_low_lim))
    data, header, w, raImage, decImage, bulge = img(directory, name, crop)
    Q, chosen_survey, mag_low_cutoff = query(raImage, decImage, band, given_catalog_path=given_catalog, mag_lower_lim=mag_low_lim,
                                             mag_upper_lim=mag_high_lim, bulge=bulge)
    psfcatalogName = []
    for f in os.listdir(directory):
        if f.endswith('.psf.cat'):
            psfcatalogName.append(f)
    psfcatalogName = ''.join(psfcatalogName)
    good_cat_stars, cleanPSFSources, PSFSources, idx_psfmass, idx_psfimage, massCatCoords = tables(Q, data, w, psfcatalogName,
                                                                                    crop, given_catalog)
    cleanPSFSources, PSFsources, psfweights_noclip, psf_clipped, ab_cat_stars = zeropt(good_cat_stars, cleanPSFSources, PSFSources,
                                                                         idx_psfmass, idx_psfimage,
                                                                         name, band, chosen_survey, sigma)
    if grb_ra:
        GRB(grb_ra, grb_dec, name, chosen_survey, band, grb_radius, massCatCoords, ab_cat_stars, directory)
    elif grb_coordlist:
        GRB(grb_ra, grb_dec, name, chosen_survey, band, grb_radius, massCatCoords, ab_cat_stars, directory, grb_coordlist)
    slope, intercept = photometry_plots(cleanPSFSources, PSFsources, data, name, chosen_survey, band, ab_cat_stars, idx_psfmass,
                                        idx_psfimage, psfweights_noclip, psf_clipped, sigma)

    return intercept

# %% optional removal of intermediate files


def removal(directory):
    fnames = ['.cat', '.psf']
    for f in os.listdir(directory):
        for name in fnames:
            if f.endswith(name):
                path = os.path.join(directory + f)
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
        grb_radius=defaults['grb_radius'], int_cal=defaults['int_cal']
):
    start_time = dt.now()

    det_thresh = 0.15   # percentage of median of weight image to cut sources under (ex. det_thresh = 0.15 -> cutoff = med * 0.15)
    max_int = 0.1   # max int value allowed for photometric fit

    directory = os.path.dirname(full_filename)
    if directory == '':
        directory = '.'
    directory = directory + '/'
    name = os.path.basename(full_filename)

    if grb_ra:
        grb_thresh = grb_rad_convert(grb_radius)

    if grb_only:
        if grb_dec is None:
            sys.exit('Only GRB RA is found, GRB Dec is None!  Make sure the -grb_dec flag is correctly formatted!')
        os.chdir(directory)
        data, header, w, raImage, decImage, bulge = img(directory, name, crop)
        Q, chosen_survey, mag_low_cutoff = query(raImage, decImage, band, survey, given_catalog, mag_low_lim,
                                                 mag_high_lim, bulge)
        psfcatalogName = name.replace('.fits', '.psf.cat')
        good_cat_stars, cleanPSFSources, PSFSources, idx_psfmass, idx_psfimage, massCatCoords = tables(Q, data, w, psfcatalogName,
                                                                                        crop, given_catalog)
        if grb_coordlist:
            GRB(grb_ra, grb_dec, name, chosen_survey, band, grb_thresh, massCatCoords, good_cat_stars, directory, grb_coordlist)
        else:
            if grb_thresh > 60:
                # newsourcesearch(grb_ra, grb_dec, w, name, chosen_survey, band, massCatCoords, grb_thresh)
                newsourcesearch(grb_ra, grb_dec, grb_thresh, w, name, chosen_survey, band, Q, data, crop)
            else:
                GRB(grb_ra, grb_dec, name, chosen_survey, band, grb_thresh, massCatCoords, good_cat_stars, directory)
    else:
        data, header, w, raImage, decImage, bulge = img(directory, name, crop)
        Q, chosen_survey, mag_low_cutoff = query(raImage, decImage, band, survey, given_catalog, mag_low_lim,
                                                 mag_high_lim, bulge)
        catalogName = sex1(name, det_cut=det_thresh)
        psfex(catalogName)
        psfcatalogName = sex2(name, det_cut=det_thresh)
        good_cat_stars, cleanPSFSources, PSFSources, idx_psfmass, idx_psfimage, massCatCoords = tables(Q, data, w, psfcatalogName,
                                                                                        crop, given_catalog)
        if len(idx_psfimage) == 0:
            print('No crossmatches found!  Cannot continue with photometry!  Is there something wrong with the image, '
                  'source catalogs, or psf model?  If those all seem normal, perhaps the image has had pixel values scaled'
                  'to uJy.  The current setup only applies the uJy/ADU conv factor as a header card, rerun the image '
                  'stacking and try again!')
        else:
            cleanPSFSources, PSFsources, psfweights_noclip, psf_clipped, ab_cat_stars = zeropt(good_cat_stars, cleanPSFSources, PSFSources, idx_psfmass, idx_psfimage,
                                                 name, band, chosen_survey, sigma)
            # if 'BUNIT' not in header:
            #     psfcatalogName = sex2(name, det_cut=det_thresh)
            #     good_cat_stars, cleanPSFSources, PSFSources, idx_psfmass, idx_psfimage, massCatCoords = tables(Q, data, w,
            #                                                                                                    psfcatalogName,
            #                                                                                                    crop,
            #                                                                                                    given_catalog)
            #     cleanPSFSources, PSFsources, psfweights_noclip, psf_clipped, ab_cat_stars = zeropt(good_cat_stars, cleanPSFSources, PSFSources, idx_psfmass, idx_psfimage,
            #                                          name, band, chosen_survey, sigma)
            if grb_ra:
                if grb_dec is None:
                    sys.exit('Only GRB RA is found, GRB Dec is None!  Make sure the -grb_dec flag is correctly formatted!')
                if grb_thresh > 60:
                    newsourcesearch(grb_ra, grb_dec, grb_thresh, w, name, chosen_survey, band, Q, data, crop)
                else:
                    GRB(grb_ra, grb_dec, name, chosen_survey, band, grb_thresh, massCatCoords, ab_cat_stars, directory)
            elif grb_coordlist:
                GRB(grb_ra, grb_dec, name, chosen_survey, band, grb_thresh, massCatCoords, ab_cat_stars, directory, grb_coordlist)
            if not no_plots:
                slope, intercept = photometry_plots(cleanPSFSources, PSFsources, data,  name, chosen_survey, band, ab_cat_stars, idx_psfmass,
                                 idx_psfimage, psfweights_noclip, psf_clipped, sigma)

            prev_intercept = intercept
            revert_flag = False

            while abs(intercept) > max_int:
                print('\nIntercept = %.4f\n' % intercept)
                # sigma -= 0.5
                mag_low_cutoff += 0.5
                new_intercept = int_calibration(name, directory, band, crop, sigma, given_catalog, chosen_survey,
                                                mag_low_cutoff, mag_high_lim,  grb_ra, grb_dec, grb_coordlist, grb_thresh,
                                                max_int=max_int)
                if abs(new_intercept) > abs(prev_intercept):
                    print("\nNew intercept: %.4f is higher than previous: %.4f! Reverting and "
                          "redoing...\n" % (new_intercept, prev_intercept))
                    intercept = prev_intercept
                    mag_low_cutoff -= 0.5
                    new_intercept = int_calibration(name, directory, band, crop, sigma, given_catalog, chosen_survey,
                                                    mag_low_cutoff, mag_high_lim, grb_ra,
                                                    grb_dec, grb_coordlist, grb_thresh, max_int=max_int)
                    revert_flag = True
                    break
                else:
                    intercept = new_intercept
                    prev_intercept = intercept
                    revert_flag = False
            if revert_flag:
                print("Loop stopped due to intercept reverting to the previous value: %.4f" % intercept)
            else:
                print(f"Final intercept below 0.15: %.4f" % intercept)

            prev_intercept = intercept
            revert_flag = False

            while abs(intercept) > max_int and sigma > 0:  # ensure sigma doesn't go negative
                print('\nIntercept = %.4f\n' % intercept)
                step = 0.25 if sigma <= 1.5 else 0.5
                sigma -= step
                new_intercept = int_calibration(name, directory, band, crop, sigma, given_catalog, chosen_survey,
                                                mag_low_cutoff, mag_high_lim, grb_ra, grb_dec, grb_coordlist,
                                                grb_thresh,
                                                max_int=max_int)
                if abs(new_intercept) > abs(prev_intercept):
                    print("\nNew intercept: %.4f is higher than previous: %.4f! Reverting and "
                          "redoing...\n" % (new_intercept, prev_intercept))
                    # revert
                    intercept = prev_intercept
                    sigma += step
                    new_intercept = int_calibration(name, directory, band, crop, sigma, given_catalog,
                                                    chosen_survey,
                                                    mag_low_cutoff, mag_high_lim, grb_ra,
                                                    grb_dec, grb_coordlist, grb_thresh, max_int=max_int)
                    revert_flag = True
                    break
                else:
                    intercept = new_intercept
                    prev_intercept = intercept
                    revert_flag = False

            if revert_flag:
                print("Sigma loop stopped due to intercept reverting to the previous value: %.4f" % intercept)
            else:
                print(f"Final intercept after sigma tuning: %.4f" % intercept)

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
    parser.add_argument('-int_cal', action='store_true',
                        help='optional flag, use to automatically improve 3 sigma fit y-int.  When y-int is >0.15, the '
                             'low mag cutoff value is increased by 0.5, only stopping when y-int < 0.15.',
                        default=defaults["int_cal"])
    args, unknown = parser.parse_known_args()
    # print(args)
    # print(unknown)

    photometry(args.filepath, args.band, args.crop, args.sigma, args.catalog, args.survey, args.mag_low,
               args.mag_high, args.no_plots, args.keep,
               args.grb_only, args.grb_ra, args.grb_dec, args.grb_coordlist, args.grb_radius, args.int_cal)


if __name__ == "__main__":
    main()
