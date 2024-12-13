"""
Calibrates photometry for stacked image
"""

import numpy as np
import numpy.ma as ma
import argparse
import sys
import astropy.units as u
from astroquery.vizier import Vizier
from astroquery.ipac.ned import Ned
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.wcs import utils
from astropy.stats import sigma_clip, sigma_clipped_stats
from astropy.io import fits
from astropy.io import ascii
import os
from astropy.table import Column
from astropy.table import Table
import matplotlib.pyplot as plt
import statsmodels.api as sm
import subprocess
from scipy.stats import skew

# sys.path.insert(0, 'C:\PycharmProjects\prime-photometry\photomitrus')
from photomitrus.settings import (gen_config_file_name, PHOTOMETRY_MAG_LOWER_LIMIT, PHOTOMETRY_MAG_UPPER_LIMIT,
                                  PHOTOMETRY_QUERY_WIDTH, PHOTOMETRY_QUERY_CATALOGS, PHOTOMETRY_LIM_MAGS)

# %%
defaults = dict(crop=500, RA=None, DEC=None, thresh='4.0', sigma=3)


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


# %%
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

    return data, header, w, raImage, decImage


# %%
# Use astroquery to get catalog search


def query(raImage, decImage, band, survey=None, given_catalog_path=None, mag_lower_lim=None, mag_upper_lim=None):
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

            coords = SkyCoord(ra=[raImage], dec=[decImage], unit=(u.deg, u.deg))

            # current columns
            v = Vizier(columns=['RAJ2000', 'DEJ2000', 'RAICRS', 'DEICRS', 'RA_ICRS', 'DE_ICRS', '%sap3' % band,
                                'e_%sap3' % band, '%smag' % band, 'e_%smag' % band, '%smag' % band.lower(),
                                'e_%smag' % band.lower(), '%sPSF' % band.lower(), 'e_%sPSF' % band.lower(),
                                '%spmag' % band.lower(), 'e_%spmag' % band.lower()])
            try:
                result = v.query_region(coords, width=str(width) + 'm', catalog=[f[1] for f in catalogs])
                test = result[0]
            except IndexError:
                print('Sadly, no current surveys available in current area in %s band' % band)
                sys.exit('No surveys available.')

            keys = result.format_table_list()

            match_keys = [(f[0], f[1]) for f in catalogs if f[1] in keys]

            print('RA: %.4f, DEC: %.4f, Box Width: %s arcmin... '
                  '\n%s band initial query resulting in: \n%s' % (raImage, decImage, width, band, keys))

            keycheck = result.keys()

            for f in catalogs:
                for k in keycheck:
                    if f[1] in k:
                        print('%s catalog found!' % k)
                        vhs_table = result[''.join(k)]
                        cols = vhs_table.colnames
                        vhs_band_col = vhs_table[cols[3]]
                        if not np.all(vhs_band_col.mask):
                            print('survey has coverage in %s band!' % band)
                            catNum = ''.join(k)
                            chosen_survey = f[0]
                            print('Survey = %s' % chosen_survey)
                            if chosen_survey == 'DES_Y' or chosen_survey == 'DES_Z' or chosen_survey == 'Skymapper':
                                print('\nMag upper limit active due to survey choice: %s' % mag_high_cutoff)
                                print(
                                    '\nQuerying Vizier %s around RA %.4f, Dec %.4f with a box of width %.2f arcmin' % (
                                        catNum, raImage, decImage, width))
                                try:
                                    v = Vizier(columns=['%s' % cols[1], '%s' % cols[2], '%s' % cols[3], '%s' % cols[4]],
                                               column_filters={"%s" % cols[3]: f"{mag_low_cutoff:f}..{mag_high_cutoff:f}",
                                                               "%sFlag" % band.lower(): "<4"
                                                               }, row_limit=-1)
                                    Q = v.query_region(SkyCoord(ra=raImage, dec=decImage, unit=(u.deg, u.deg)),
                                                       width=str(width) + 'm'
                                                       , catalog=catNum, cache=False)
                                    print('Queried source total = ', len(Q[0]))
                                except Exception as e:
                                    print(
                                        'Error in Vizier query.')
                            else:
                                print('\nQuerying Vizier %s around RA %.4f, Dec %.4f with a box of width %.2f arcmin' % (
                                    catNum, raImage, decImage, width))
                                try:
                                    v = Vizier(columns=['%s' % cols[1], '%s' % cols[2], '%s' % cols[3], '%s' % cols[4]],
                                               column_filters={"%s" % cols[3]: f">{mag_low_cutoff:f}"}, row_limit=-1)
                                    Q = v.query_region(SkyCoord(ra=raImage, dec=decImage, unit=(u.deg, u.deg)), width=str(width) + 'm'
                                                       , catalog=catNum, cache=False)
                                    print('Queried source total = ', len(Q[0]))
                                except Exception as e:
                                    print(
                                        'Error in Vizier query.')
                                    print(f"Error details: {e}")
                            break
                        else:
                            print('no %s mag sources found, defaulting to next survey...' % band)
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
                print('\nQuerying Vizier %s around RA %.4f, Dec %.4f with a box width %.3f arcmin' % (
                    catNum, raImage, decImage, width))
                try:
                    # You can set the bands for the individual columns (magnitude range, number of detections) inside the Vizier query
                    v = Vizier(columns=['RAJ2000', 'DEJ2000', '%smag' % band, 'e_%smag' % band],
                               column_filters={"%smag" % band: f">{mag_low_cutoff:f}", "Nd": ">6"},
                               row_limit=-1)
                    Q = v.query_region(SkyCoord(ra=raImage, dec=decImage, unit=(u.deg, u.deg)), width=str(width) + 'm'
                                       , catalog=catNum, cache=False)
                    # query vizier around (ra, dec) with a radius of boxsize
                    # print(Q[0])
                    print('Queried source total = ', len(Q[0]))
                except:
                    print('I cannnot reach the Vizier database. Is the internet working?')
            elif survey == 'VHS':
                catNum = 'II/367'  # changing to vista
                print('\nQuerying Vizier %s around RA %.4f, Dec %.4f with a box width %.3f arcmin' % (
                    catNum, raImage, decImage, width))
                try:
                    # You can set the bands for the individual columns (magnitude range, number of detections) inside the Vizier query
                    v = Vizier(columns=['RAJ2000', 'DEJ2000', '%sap3' % band, 'e_%sap3' % band],
                               column_filters={"%sap3" % band: f">{mag_low_cutoff:f}"}, row_limit=-1)
                    Q = v.query_region(SkyCoord(ra=raImage, dec=decImage, unit=(u.deg, u.deg)), width=str(width) + 'm'
                                       , catalog=catNum, cache=False)
                    # query vizier around (ra, dec) with a radius of boxsize
                    # print(Q[0])
                    print('Queried source total = ', len(Q[0]))
                except:
                    print(
                        'Error in Vizier query. Perhaps your image is not in the southern hemisphere sky?  H band is also not well covered!')
            elif survey == 'VIKING':
                catNum = 'II/343/viking2'
                print('\nQuerying Vizier %s around RA %.4f, Dec %.4f with a box width %.3f arcmin' % (
                    catNum, raImage, decImage, width))
                try:
                    v = Vizier(columns=['RAJ2000', 'DEJ2000', '%sap3' % band, 'e_%sap3' % band],
                               column_filters={"%sap3" % band: f">{mag_low_cutoff:f}"}, row_limit=-1)
                    Q = v.query_region(SkyCoord(ra=raImage, dec=decImage, unit=(u.deg, u.deg)), width=str(width) + 'm'
                                       , catalog=catNum, cache=False)
                    print('Queried source total = ', len(Q[0]))
                except:
                    print(
                        'Error in Vizier query. Perhaps your image is not in the southern hemisphere sky?  H band is also not well covered!'
                        ' If you are in S.H., VIKING is only in a relatively smaller strip!')
            elif survey == 'Skymapper':
                catNum = 'II/379/smssdr4'
                print('\nQuerying Vizier %s around RA %.4f, Dec %.4f with a box width %.3f arcmin' % (
                    catNum, raImage, decImage, width))
                try:
                    v = Vizier(columns=['RAICRS', 'DEICRS', '%sPSF' % band.lower(), 'e_%sPSF' % band.lower()],
                               column_filters={"%sPSF" % band.lower(): f">{mag_low_cutoff:f}"}, row_limit=-1)
                    Q = v.query_region(SkyCoord(ra=raImage, dec=decImage, unit=(u.deg, u.deg)), width=str(width) + 'm'
                                       , catalog=catNum, cache=False)
                    print('Queried source total = ', len(Q[0]))
                except:
                    print(
                        'Error in Vizier query. Perhaps your image is not in the southern hemisphere sky?'
                        '\n perhaps check skymapper coverage maps?')
            elif survey == 'SDSS':
                catNum = 'V/154/sdss16'
                print('\nQuerying Vizier %s around RA %.4f, Dec %.4f with a box width %.3f arcmin' % (
                    catNum, raImage, decImage, width))
                try:
                    v = Vizier(columns=['RA_ICRS', 'DE_ICRS', '%spmag' % band.lower(), 'e_%spmag' % band.lower()],
                               column_filters={"%spmag" % band.lower(): f">{mag_low_cutoff:f}"
                                   , "clean": "=1"}, row_limit=-1)
                    Q = v.query_region(SkyCoord(ra=raImage, dec=decImage, unit=(u.deg, u.deg)), width=str(width) + 'm'
                                       , catalog=catNum, cache=False)
                    print('Queried source total = ', len(Q[0]))
                except:
                    print(
                        'Error in Vizier query. Perhaps your image is not in the southern hemisphere sky?'
                        '\n perhaps check SDSS coverage maps?')
            elif survey == 'UKIDSS':
                catNum = 'II/319/las9'
                print('\nQuerying Vizier %s around RA %.4f, Dec %.4f with a box width %.3f arcmin' % (
                    catNum, raImage, decImage, width))
                try:
                    v = Vizier(columns=['RAJ2000', 'DEJ2000', '%smag' % band, 'e_%smag' % band],
                               column_filters={"%smag" % band: f">{mag_low_cutoff:f}"}, row_limit=-1)
                    Q = v.query_region(SkyCoord(ra=raImage, dec=decImage, unit=(u.deg, u.deg)), width=str(width) + 'm',
                                       catalog=catNum, cache=False)
                    print('Queried source total = ', len(Q[0]))
                except:
                    print(
                        'Error in Vizier query. Perhaps your image is not in the southern hemisphere sky?'
                        '\n perhaps check UKIDSS coverage maps?')
            elif survey == 'DES':
                if band == 'Y':
                    catNum = 'II/371/des_dr2'
                    print('\nQuerying Vizier %s around RA %.4f, Dec %.4f with a box width %.3f arcmin' % (
                        catNum, raImage, decImage, width))
                    try:
                        v = Vizier(columns=['RA_ICRS', 'DE_ICRS', '%smag' % band, 'e_%smag' % band],
                                   column_filters={"%smag" % band: f"{mag_low_cutoff:f}..{mag_high_cutoff:f}"}, row_limit=-1)
                        Q = v.query_region(SkyCoord(ra=raImage, dec=decImage, unit=(u.deg, u.deg)), width=str(width) + 'm',
                                           catalog=catNum, cache=False)
                        print('Queried source total = ', len(Q[0]))
                    except:
                        print(
                            'Error in Vizier query. Perhaps your image is not in the southern hemisphere sky?'
                            '\n perhaps check DES coverage maps?')
                elif band == 'Z':
                    catNum = 'II/371/des_dr2'
                    print('\nQuerying Vizier %s around RA %.4f, Dec %.4f with a box width %.3f arcmin' % (
                        catNum, raImage, decImage, width))
                    try:
                        v = Vizier(columns=['RA_ICRS', 'DE_ICRS', '%smag' % band.lower(), 'e_%smag' % band.lower()],
                                   column_filters={"%smag" % band: f"{mag_low_cutoff:f}..{mag_high_cutoff:f}"}, row_limit=-1)
                        Q = v.query_region(SkyCoord(ra=raImage, dec=decImage, unit=(u.deg, u.deg)), width=str(width) + 'm',
                                           catalog=catNum, cache=False)
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


def sex1(imageName):
    print('Running sextractor on img to initially find sources...')
    configFile = gen_config_file_name('sex2.config')
    paramName = gen_config_file_name('tempsource.param')
    catalogName = imageName + '.cat'
    weightName = 'weight'+imageName[5:]
    if os.path.isfile(weightName):
        try:
            print('Including weight map!')
            command = ('sex %s -c %s -CATALOG_NAME %s -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE %s -PARAMETERS_NAME %s' %
                       (imageName, configFile, catalogName, weightName, paramName))
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
    try:
        command = 'psfex %s -c %s' % (catalogName, psfConfigFile)
        # print('Executing command: %s' % command)
        subprocess.run(command.split(), check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError as err:
        print('Could not run psfex with exit error %s' % err)


# %%
# feed generated psf model back into sextractor w/ diff param (or could use that param from the start but its slower)


def sex2(imageName):
    print('Feeding psf model back into sextractor for fitting and flux calculation...')
    psfName = imageName + '.psf'
    psfcatalogName = imageName.replace('.fits', '.psf.cat')
    configFile = gen_config_file_name('sex2.config')
    psfparamName = gen_config_file_name('photomPSF.param')
    weightName = 'weight' + imageName[5:]
    if os.path.isfile(weightName):
        try:
            # We are supplying SExtactor with the PSF model with the PSF_NAME option
            command = 'sex %s -c %s -CATALOG_NAME %s -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE %s -PSF_NAME %s -PARAMETERS_NAME %s' % (
            imageName, configFile, psfcatalogName, weightName, psfName, psfparamName)
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

        PSFSources = psfsourceTable[
            (psfsourceTable['XMODEL_IMAGE'] < (max_x - crop)) & (psfsourceTable['XMODEL_IMAGE'] > crop)
            & (psfsourceTable['YMODEL_IMAGE'] < (max_y) - crop) & (psfsourceTable['YMODEL_IMAGE'] > crop) &
            (psfsourceTable['FLUX_RADIUS'] >= 1.0)]

        cleanPSFSources = psfsourceTable[
            (psfsourceTable['FLAGS'] == 0) & (psfsourceTable['FLAGS_MODEL'] == 0) & (psfsourceTable['XMODEL_IMAGE']
            < (max_x - crop)) & (psfsourceTable['XMODEL_IMAGE'] > crop) & (psfsourceTable['YMODEL_IMAGE']
            < (max_y) - crop) & (psfsourceTable['YMODEL_IMAGE'] > crop) & (psfsourceTable['FLUX_RADIUS'] >= 1.0)]

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

        PSFSources = psfsourceTable[
            (psfsourceTable['XMODEL_IMAGE'] < (max_x - crop)) & (psfsourceTable['XMODEL_IMAGE'] > crop)
            & (psfsourceTable['YMODEL_IMAGE'] < (max_y) - crop) & (psfsourceTable['YMODEL_IMAGE'] > crop)
            & (psfsourceTable['FLUX_RADIUS'] >= 1.0)]

        cleanPSFSources = psfsourceTable[
            (psfsourceTable['FLAGS'] == 0) & (psfsourceTable['FLAGS_MODEL'] == 0) & (psfsourceTable['XMODEL_IMAGE']
                                                                                     < (max_x - crop)) & (
                        psfsourceTable['XMODEL_IMAGE'] > crop) & (psfsourceTable['YMODEL_IMAGE'] < (max_y) - crop)
            & (psfsourceTable['YMODEL_IMAGE'] > crop) & (psfsourceTable['FLUX_RADIUS'] >= 1.0)]

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


def ab_convert(mag, band):
    # jy zero points
    # TODO add more survey-specific zps once query is overhauled (ex. VISTA, etc)
    # from 2MASS
    zp_dict = {'zp_J': 1594, 'zp_H': 1024}
    pick_zp = 'zp_' + band

    for k, v in zp_dict.items():
        if k == pick_zp:
            zp = v

    flx = zp * 10 ** (-mag / 2.5)
    ab_mag = -2.5 * np.log10(flx / 3631)
    ab_mag = round(ab_mag, 3)
    return ab_mag


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
    # catalog for just clean sources (no flags)
    psfmag_clean = zero_psfmean + cleanPSFSources['MAG_POINTSOURCE']
    psfmagerr_clean = np.sqrt(cleanPSFSources['MAGERR_POINTSOURCE'] ** 2 + zero_psfstd ** 2)

    psfmagcol_clean = Column(psfmag_clean, name='%sMAG_PSF' % band, unit='mag')
    psfmagerrcol_clean = Column(psfmagerr_clean, name='e_%sMAG_PSF' % band, unit='mag')
    cleanPSFSources.add_column(psfmagcol_clean)
    cleanPSFSources.add_column(psfmagerrcol_clean)
    cleanPSFSources.remove_column('VIGNET')

    # catalog for all detected sources
    psfmag = zero_psfmean + PSFSources['MAG_POINTSOURCE']

    # caterr_all = good_cat_stars[magerrcolname]
    # primeerr_all = cleanPSFSources['MAGERR_POINTSOURCE']

    # np.sqrt(caterr_all**2 + primeerr_all**2)
    psfmagerr = np.sqrt(PSFSources['MAGERR_POINTSOURCE'] ** 2 + zero_psfstd ** 2)

    psfmagcol = Column(psfmag, name='%sMAG_PSF' % band, unit='mag')
    psfmagerrcol = Column(psfmagerr, name='e_%sMAG_PSF' % band, unit='mag')
    PSFSources.add_column(psfmagcol)
    PSFSources.add_column(psfmagerrcol)
    PSFSources.remove_column('VIGNET')
    print('Total PRIME source # = ', len(PSFSources))

    PSFSources.write('%s.%s.ecsv' % (imageName, survey), overwrite=True)
    print('%s.%s.ecsv written, CSV w/ corrected mags' % (imageName, survey))

    return cleanPSFSources, PSFSources, psfweights_noclip, psf_clipped


#%%


# New Source Search
def newsourcesearch(source_ra, source_dec, w, imageName, survey, band, massCatCoords, thresh):
    # crossmatch for all detected sources
    print('Large grb radius inputted, using new source search to find'
          ' detected sources brighter than survey lim mag w/ no crossmatch')

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


def GRB(ra, dec, imageName, survey, band, thresh, coordlist=None):
    mag_ecsvname = '%s.%s.ecsv' % (imageName, survey)
    mag_ecsvtable = ascii.read(mag_ecsvname)
    mag_ecsvcleanSources = mag_ecsvtable  # [(mag_ecsvtable['FLAGS'] == 0) & (mag_ecsvtable['FLAGS_MODEL'] == 0)]
    mag_ecsvsourceCatCoords = SkyCoord(ra=mag_ecsvcleanSources['ALPHA_J2000'], dec=mag_ecsvcleanSources['DELTA_J2000'],
                                       frame='icrs',
                                       unit='degree')

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
                idx_GRBpsfdict[keys[i]] = idx_GRBcleanpsf, coordlist[i]
        except NameError:
            print('No Sources found!')
    else:
        GRBcoords = SkyCoord(ra=[ra], dec=[dec], frame='icrs', unit='degree')
        idx_GRB, idx_GRBcleanpsf, d2d, d3d = mag_ecsvsourceCatCoords.search_around_sky(GRBcoords,
                                                                                       photoDistThresh * u.arcsec)

    if coordlist:
        for key in idx_GRBpsfdict:
            values = idx_GRBpsfdict[key]
            idx_GRBcleanpsf = values[0]
            ra = values[1][0]
            dec = values[1][1]
            print('idx size = %d for location %d' % (len(idx_GRBcleanpsf), key))
            if len(idx_GRBcleanpsf) == 1:
                print('GRB source at inputted coords %s and %s, rad = %s arcsec found!' % (ra, dec, photoDistThresh))

                grb_mag = mag_ecsvcleanSources[idx_GRBcleanpsf]['%sMAG_PSF' % band][0]
                grb_magerr = mag_ecsvcleanSources[idx_GRBcleanpsf]['e_%sMAG_PSF' % band][0]

                grb_ra = mag_ecsvcleanSources[idx_GRBcleanpsf]['ALPHA_J2000'][0]
                grb_dec = mag_ecsvcleanSources[idx_GRBcleanpsf]['DELTA_J2000'][0]
                grb_rad = mag_ecsvcleanSources[idx_GRBcleanpsf]['FLUX_RADIUS'][0]
                grb_snr = mag_ecsvcleanSources[idx_GRBcleanpsf]['SNR_WIN'][0]
                grb_dist = d2d[0].to(u.arcsec)
                grb_dist = grb_dist / u.arcsec

                print('Detected GRB ra = %.6f, dec = %.6f, with 50 percent flux radius = %.3f arcsec and SNR = %.3f' % (
                    grb_ra, grb_dec, grb_rad, grb_snr))
                print('%s magnitude of GRB is %.2f +/- %.2f' % (band, grb_mag, grb_magerr))

                grbdata = Table()
                grbdata['RA (deg)'] = np.array([grb_ra])
                grbdata['DEC (deg)'] = np.array([grb_dec])
                grbdata['%sMag' % band] = np.array([grb_mag])
                grbdata['%sMag_Err' % band] = np.array([grb_magerr])
                grbdata['Radius (arcsec)'] = np.array([grb_rad])
                grbdata['SNR'] = np.array([grb_snr])
                grbdata['Distance (arcsec)'] = np.array([grb_dist])

                grbdata.write('GRB_%s_Data_%s_loc_%d.ecsv' % (band, survey, key), overwrite=True)
                print('Generated GRB data table!')
            elif len(idx_GRBcleanpsf) > 1:
                print('Multiple sources detected in search radius (ra = %.6f, dec = %.6f, rad = %s arcsec)'
                      ', refer to .ecsv file for source info!' % (ra, dec, photoDistThresh))
                mag_ar = []
                mag_err_ar = []
                ra_ar = []
                dec_ar = []
                rad_ar = []
                snr_ar = []
                dist_ar = []
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
                grbdata = Table()
                grbdata['RA (deg)'] = np.array(ra_ar)
                grbdata['DEC (deg)'] = np.array(dec_ar)
                grbdata['%sMag' % band] = np.array(mag_ar)
                grbdata['%sMag_Err' % band] = np.array(mag_err_ar)
                grbdata['Radius (arcsec)'] = np.array(rad_ar)
                grbdata['SNR'] = np.array(snr_ar)
                grbdata['Distance (arcsec)'] = np.array(dist_ar)
                grbdata.write('GRB_Multisource_%s_Data_%s_loc_%d.ecsv' % (band, survey, key), overwrite=True)
                print('Generated GRB data table!')
            else:
                print(
                    'GRB source at inputted coords %s and %s not found, perhaps increase photoDistThresh?' % (ra, dec))
    else:
        print('idx size = %d' % len(idx_GRBcleanpsf))
        if len(idx_GRBcleanpsf) == 1:
            print('GRB source at inputted coords %s and %s, rad = %s arcsec found!' % (ra, dec, photoDistThresh))

            grb_mag = mag_ecsvcleanSources[idx_GRBcleanpsf]['%sMAG_PSF' % band][0]
            grb_magerr = mag_ecsvcleanSources[idx_GRBcleanpsf]['e_%sMAG_PSF' % band][0]

            grb_ra = mag_ecsvcleanSources[idx_GRBcleanpsf]['ALPHA_J2000'][0]
            grb_dec = mag_ecsvcleanSources[idx_GRBcleanpsf]['DELTA_J2000'][0]
            grb_rad = mag_ecsvcleanSources[idx_GRBcleanpsf]['FLUX_RADIUS'][0]
            grb_snr = mag_ecsvcleanSources[idx_GRBcleanpsf]['SNR_WIN'][0]
            grb_dist = d2d[0].to(u.arcsec)
            grb_dist = grb_dist / u.arcsec

            print('Detected GRB ra = %.6f, dec = %.6f, with 50 percent flux radius = %.3f arcsec and SNR = %.3f' % (
                grb_ra, grb_dec, grb_rad, grb_snr))
            print('%s magnitude of GRB is %.2f +/- %.2f' % (band, grb_mag, grb_magerr))

            grbdata = Table()
            grbdata['RA (deg)'] = np.array([grb_ra])
            grbdata['DEC (deg)'] = np.array([grb_dec])
            grbdata['%sMag' % band] = np.array([grb_mag])
            grbdata['%sMag_Err' % band] = np.array([grb_magerr])
            grbdata['Radius (arcsec)'] = np.array([grb_rad])
            grbdata['SNR'] = np.array([grb_snr])
            grbdata['Distance (arcsec)'] = np.array([grb_dist])

            grbdata.write('GRB_%s_Data_%s.ecsv' % (band, survey), overwrite=True)
            print('Generated GRB data table!')
        elif len(idx_GRBcleanpsf) > 1:
            print('Multiple sources detected in search radius (ra = %.6f, dec = %.6f, rad = %s arcsec)'
                  ', refer to .ecsv file for source info!' % (ra, dec, photoDistThresh))
            mag_ar = []
            mag_err_ar = []
            ra_ar = []
            dec_ar = []
            rad_ar = []
            snr_ar = []
            dist_ar = []
            idx_GRBcleanpsflist = idx_GRBcleanpsf.tolist()
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
            grbdata = Table()
            grbdata['RA (deg)'] = np.array(ra_ar)
            grbdata['DEC (deg)'] = np.array(dec_ar)
            grbdata['%sMag' % band] = np.array(mag_ar)
            grbdata['%sMag_Err' % band] = np.array(mag_err_ar)
            grbdata['Radius (arcsec)'] = np.array(rad_ar)
            grbdata['SNR'] = np.array(snr_ar)
            grbdata['Distance (arcsec)'] = np.array(dist_ar)
            grbdata.write('GRB_Multisource_%s_Data_%s.ecsv' % (band, survey), overwrite=True)
            print('Generated GRB data table!')
        else:
            print('GRB source at inputted coords %s and %s not found, perhaps increase photoDistThresh?' % (ra, dec))


# %% optional plots


def photometry_plots(cleanPSFsources, PSFsources, data, imageName, survey, band, good_cat_stars, idx_psfmass, idx_psfimage,
                     psfweights_noclip, psf_clipped, sigma):
    # appropriate mag column
    colnames = good_cat_stars.colnames
    magcol = colnames[2]

    chip = imageName[-6]
    if len(imageName) <= 16:
        num = 'img'
    else:
        num = imageName[-16:-8]

    # mag comparison plot
    plt.figure(1, figsize=(8, 8))
    plt.clf()
    plt.plot(cleanPSFsources['%sMAG_PSF' % band][idx_psfimage], good_cat_stars['%s' % magcol][idx_psfmass],
             'r.', markersize=14, markeredgecolor='black')
    plt.xlim(10, 22)
    plt.ylim(10, 22)
    plt.title('PRIME Mags vs %s Mags' % survey)
    plt.xlabel('PRIME %s Mags' % band, fontsize=15)
    plt.ylabel('%s %s Mags' % (survey, band), fontsize=15)
    plt.grid()
    plt.savefig('%s_C%s_mag_comp_plot_%s.png' % (survey, chip, num))
    print('Saved mag comparison plot to dir!')

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
    x_sig = cleanPSFsources['%sMAG_PSF' % band][idx_psfimage][~psf_clipped.mask]
    y_sig = good_cat_stars['%s' % magcol][idx_psfmass][~psf_clipped.mask]
    x_const_sig = sm.add_constant(x_sig)
    model_sig = sm.WLS(y_sig, x_const_sig, weights=psfweights_noclip[~psf_clipped.mask]).fit()
    m_sig = model_sig.params[1]
    m_sigerr = model_sig.bse[1]
    b_sig = model_sig.params[0]
    b_sigerr = model_sig.bse[0]
    rsquare_sig = model_sig.rsquared
    rss_sig = model_sig.ssr

    print('# of crossmatched sources used in %s sig fit: %i'
          % (sigma, len(cleanPSFsources['%sMAG_PSF' % band][idx_psfimage][~psf_clipped.mask])))
    print('%s sig fit: slope = %.4f +/- %.4f' % (sigma, m_sig, m_sigerr))
    print('%s sig fit: y-int = %.4f +/- %.4f' % (sigma, b_sig, b_sigerr))

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
        res_mask = model_sig.resid[mask]
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

    # res plot, y int included
    plt.figure(2, figsize=(8, 6))
    plt.clf()
    plt.scatter(cleanPSFsources['%sMAG_PSF' % band][idx_psfimage][~psf_clipped.mask], model_sig.resid, color='red')
    plt.ylim(-1, 1)
    plt.xlim(10, 21)
    plt.title('PRIME vs %s Residuals - %s Sigma Clip' % (survey, sigma))
    plt.ylabel('Residuals')
    plt.xlabel('%s Mags' % band)
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
    plt.savefig('%s_C%s_residual_plot_int_%s.png' % (survey, chip, num), dpi=300)

    # res plot y int, histogram
    if len(idx_psfimage) >= 5000:
        bin_num_int = round(len(idx_psfimage) / 75)
    elif 1000 <= len(idx_psfimage) <= 5000:
        bin_num_int = round(len(idx_psfimage) / 50)
    elif len(idx_psfimage) <= 1000:
        bin_num_int = 50

    plt.figure(3, figsize=(10, 6))
    plt.clf()
    plt.hist2d(x=cleanPSFsources['%sMAG_PSF' % band][idx_psfimage][~psf_clipped.mask], y=model_sig.resid,
               bins=[bin_num_int, bin_num_int], range=[[10, 21],[-1, 1]], cmap='gist_heat_r')
    plt.colorbar(label='Density')
    plt.ylim(-1, 1)
    plt.xlim(10, 21)
    plt.title('PRIME vs %s Residuals w/ %s Sigma Clip- Density Histogram' % (survey, sigma))
    plt.ylabel('Residuals')
    plt.xlabel('%s Mags' % band)
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

    print('Saved y-int residual plots to dir!')

    # WLS fit line over data plot
    def predict_y_for(x, m, b):
        return m * x + b

    txt = ('slope = %.4f' % m2 + '\nslope err = %.4f' % m2err + '\nint = %.4f' % b2 + '\nint err = %.4f' % b2err)

    plt.figure(4, figsize=(8, 8))
    plt.clf()
    plt.xlim(10, 22)
    plt.ylim(10, 22)
    plt.title('PRIME vs %s w/ Weighted Fit' % survey)
    plt.grid()
    plt.xlabel('PRIME %s Mags' % band, fontsize=15)
    plt.ylabel('%s %s Mags' % (survey, band), fontsize=15)
    plt.scatter(cleanPSFsources['%sMAG_PSF' % band][idx_psfimage], good_cat_stars['%s' % magcol][idx_psfmass])
    plt.plot(cleanPSFsources['%sMAG_PSF' % band][idx_psfimage],
             predict_y_for(cleanPSFsources['%sMAG_PSF' % band][idx_psfimage], m2, b2), c='r')
    box = dict(facecolor='white')
    plt.text(11, 18, txt, fontsize=12, bbox=box)
    plt.savefig('%s_C%s_WLS_fit_plot_%s.png' % (survey, chip, num), dpi=300)

    # WLS hist density plot
    if len(idx_psfimage) >= 5000:
        bin_num = round(len(idx_psfimage) / 50)
    elif 1000 <= len(idx_psfimage) <= 5000:
        bin_num = round(len(idx_psfimage) / 20)
    elif len(idx_psfimage) <= 1000:
        bin_num = 100

    plt.figure(5, figsize=(10, 8))
    plt.clf()
    plt.xlim(10, 22)
    plt.ylim(10, 22)
    plt.title('PRIME vs %s w/ Weighted Fit - Density Histogram' % survey)
    plt.grid()
    plt.xlabel('PRIME %s Mags' % band, fontsize=15)
    plt.ylabel('%s %s Mags' % (survey, band), fontsize=15)
    plt.hist2d(x=cleanPSFsources['%sMAG_PSF' % band][idx_psfimage], y=good_cat_stars['%s' % magcol][idx_psfmass],
               bins=[bin_num, bin_num], range=[[10, 22],[10, 22]], cmap='gist_heat_r')
    plt.plot(cleanPSFsources['%sMAG_PSF' % band][idx_psfimage],
             predict_y_for(cleanPSFsources['%sMAG_PSF' % band][idx_psfimage], m2, b2), c='b')
    plt.colorbar(label='Density')
    box = dict(facecolor='white')
    plt.text(11, 18, txt, fontsize=12, bbox=box)
    plt.savefig('%s_C%s_WLS_fit_hist_plot_%s.png' % (survey, chip, num), dpi=300)

    print('Saved WLS fit plots to dir!')

    # WLS fit for 3 sig clip of data

    sigtxt = ('slope = %.4f' % m_sig + '\nslope err = %.4f' % m_sigerr + '\nint = %.4f' % b_sig +
              '\nint err = %.4f' % b_sigerr)

    plt.figure(6, figsize=(8, 8))
    plt.clf()
    plt.xlim(10, 22)
    plt.ylim(10, 22)
    plt.title('PRIME vs %s w/ Weighted Fit - %s Sigma Clip' % (survey, sigma))
    plt.grid()
    plt.xlabel('PRIME %s Mags' % band, fontsize=15)
    plt.ylabel('%s %s Mags' % (survey, band), fontsize=15)
    plt.scatter(cleanPSFsources['%sMAG_PSF' % band][idx_psfimage][~psf_clipped.mask],
                good_cat_stars['%s' % magcol][idx_psfmass][~psf_clipped.mask])
    plt.plot(cleanPSFsources['%sMAG_PSF' % band][idx_psfimage][~psf_clipped.mask],
             predict_y_for(cleanPSFsources['%sMAG_PSF' % band][idx_psfimage][~psf_clipped.mask], m_sig, b_sig), c='r')
    box = dict(facecolor='white')
    plt.text(11, 18, sigtxt, fontsize=12, bbox=box)
    plt.savefig('%s_C%s_WLS_fit_3sig_plot_%s.png' % (survey, chip, num), dpi=300)

    # WLS 3 sig hist density plot
    plt.figure(7, figsize=(10, 8))
    plt.clf()
    plt.xlim(10, 22)
    plt.ylim(10, 22)
    plt.title('PRIME vs %s w/ Weighted Fit - %s Sigma Clip - Density Histogram' % (survey, sigma))
    plt.grid()
    plt.xlabel('PRIME %s Mags' % band, fontsize=15)
    plt.ylabel('%s %s Mags' % (survey, band), fontsize=15)
    plt.hist2d(x=cleanPSFsources['%sMAG_PSF' % band][idx_psfimage][~psf_clipped.mask],
               y=good_cat_stars['%s' % magcol][idx_psfmass][~psf_clipped.mask],
               bins=[bin_num, bin_num], range=[[10, 22],[10, 22]], cmap='gist_heat_r')
    plt.plot(cleanPSFsources['%sMAG_PSF' % band][idx_psfimage][~psf_clipped.mask],
             predict_y_for(cleanPSFsources['%sMAG_PSF' % band][idx_psfimage][~psf_clipped.mask], m_sig, b_sig), c='b')
    plt.colorbar(label='Density')
    box = dict(facecolor='white')
    plt.text(11, 18, sigtxt, fontsize=12, bbox=box)
    plt.savefig('%s_C%s_WLS_fit_3sig_hist_plot_%s.png' % (survey, chip, num), dpi=300)

    print('Saved WLS 3 sig fit plots to dir!')

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
    plt.clf()
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
    if len(PSFsources) < 1500:
        plt.ylim(0, 250)
    else:
        plt.ylim(0, 2000)
    # plt.yscale('log')
    # plt.yticks([0,0.05,0.1,0.15,0.2,0.25])
    # plt.legend([r'5 $\sigma$ Limit'],loc='upper left',fontsize=12)
    # plt.legend([f'Bin of 0.2 Error = {bin_vals[imin+1]}',f'Bin of 0.3 Error = {bin_vals[iimin+1]}'])
    # plt.legend([f'Bin of 0.2 Error (>16 mag) = {bin_vals[imin+1]}'])
    # plt.legend(['Our Field','Osaka Field'],loc='upper right')
    plt.title('PRIME Limiting Mag Plot')
    plt.ylabel('Number of Sources')
    plt.xlabel('%s Magnitude' % band)
    plt.legend(['Half Max = %s' % round(halfmax, 1), 'Limiting Mag = %s' % round(limmag, 1)], fontsize=15)
    plt.savefig('%s_C%s_lim_mag_plot_%s.png' % (survey, chip, num), dpi=300)
    print('Saved lim mag plot to dir!')

    # Crossmatch location check plot

    mean, median, sigma = sigma_clipped_stats(data)

    fig = plt.figure(figsize=(10, 10))
    ax = fig.gca()
    plt.imshow(data, vmin=median - 1.5 * sigma, vmax=median + 1.5 * sigma)
    circles = [
        plt.Circle((cleanPSFsources['X_IMAGE'][idx_psfimage][i], cleanPSFsources['Y_IMAGE'][idx_psfimage][i]), radius=5,
                   edgecolor='r', facecolor='None') for i in range(len(cleanPSFsources['X_IMAGE'][idx_psfimage]))]
    for c in circles:
        ax.add_artist(c)
    plt.savefig('%s_C%s_source_check_plot_%s.png' % (survey, chip, num), dpi=300)
    print('Saved source location check plot to dir!')

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
        name, directory, band, crop=defaults['crop'], sigma=defaults['sigma'], given_catalog=None, survey=None,
        mag_low_lim=None, mag_high_lim=None, grb_ra=None, grb_dec=None,
        grb_coordlist=None, grb_radius=4.0
):
    print('3 sigma fit y-intercept > 0.5! Redoing photometry w/ mag low cutoff = %s\n' % mag_low_lim)
    data, header, w, raImage, decImage = img(directory, name, crop)
    Q, chosen_survey, mag_low_cutoff = query(raImage, decImage, band, survey, given_catalog, mag_lower_lim=mag_low_lim,
                                             mag_upper_lim=mag_high_lim)
    psfcatalogName = []
    for f in os.listdir(directory):
        if f.endswith('.psf.cat'):
            psfcatalogName.append(f)
    psfcatalogName = ''.join(psfcatalogName)
    good_cat_stars, cleanPSFSources, PSFSources, idx_psfmass, idx_psfimage, massCatCoords = tables(Q, data, w, psfcatalogName,
                                                                                    crop, given_catalog)
    cleanPSFSources, PSFsources, psfweights_noclip, psf_clipped = zeropt(good_cat_stars, cleanPSFSources, PSFSources,
                                                                         idx_psfmass, idx_psfimage,
                                                                         name, band, chosen_survey, sigma)
    if grb_ra:
        GRB(grb_ra, grb_dec, name, survey, band, grb_radius)
    elif grb_coordlist:
        GRB(grb_ra, grb_dec, name, survey, band, grb_radius, grb_coordlist)
    slope, intercept = photometry_plots(cleanPSFSources, PSFsources, data, name, chosen_survey, band, good_cat_stars, idx_psfmass,
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
        full_filename, band, crop=defaults['crop'], sigma=defaults['sigma'], given_catalog=None, survey=None,
        mag_low_lim=None, mag_high_lim=None, no_plots=False,
        keep=False, grb_only=False, grb_ra=None, grb_dec=None, grb_coordlist=None, grb_radius=defaults['thresh'],
        int_cal=False
):
    directory = os.path.dirname(full_filename)
    if directory == '':
        directory = '.'
    directory = directory + '/'
    name = os.path.basename(full_filename)

    grb_thresh = grb_rad_convert(grb_radius)

    if grb_only:
        os.chdir(directory)
        data, header, w, raImage, decImage = img(directory, name, crop)
        Q, chosen_survey, mag_low_cutoff = query(raImage, decImage, band, survey, given_catalog, mag_low_lim,
                                                 mag_high_lim)
        if grb_coordlist:
            GRB(grb_ra, grb_dec, name, chosen_survey, band, grb_thresh, grb_coordlist)
        else:
            GRB(grb_ra, grb_dec, name, survey, band, grb_thresh)
    else:
        data, header, w, raImage, decImage = img(directory, name, crop)
        Q, chosen_survey, mag_low_cutoff = query(raImage, decImage, band, survey, given_catalog, mag_low_lim, mag_high_lim)
        catalogName = sex1(name)
        psfex(catalogName)
        psfcatalogName = sex2(name)
        good_cat_stars, cleanPSFSources, PSFSources, idx_psfmass, idx_psfimage, massCatCoords = tables(Q, data, w, psfcatalogName,
                                                                                        crop, given_catalog)
        cleanPSFSources, PSFsources, psfweights_noclip, psf_clipped = zeropt(good_cat_stars, cleanPSFSources, PSFSources, idx_psfmass, idx_psfimage,
                                             name, band, chosen_survey, sigma)
        if grb_ra:
            if grb_thresh > 60:
                newsourcesearch(grb_ra, grb_dec, w, name, chosen_survey, band, massCatCoords, grb_thresh)
            else:
                GRB(grb_ra, grb_dec, name, chosen_survey, band, grb_thresh)
        elif grb_coordlist:
            GRB(grb_ra, grb_dec, name, chosen_survey, band, grb_thresh, grb_coordlist)
        if not no_plots:
            slope, intercept = photometry_plots(cleanPSFSources, PSFsources, data,  name, chosen_survey, band, good_cat_stars, idx_psfmass,
                             idx_psfimage, psfweights_noclip, psf_clipped, sigma)
        if not int_cal:
            if not keep:
                removal(directory)
        else:
            prev_intercept = intercept
            revert_flag = False
            while intercept > 0.5:
                print('\nIntercept = %.4f\n' % intercept)
                mag_low_cutoff += 0.5
                new_intercept = int_calibration(name, directory, band, crop, sigma, given_catalog, chosen_survey,
                                                mag_low_cutoff, mag_high_lim,  grb_ra, grb_dec, grb_coordlist, grb_thresh)
                if new_intercept > prev_intercept:
                    print("\nNew intercept: %.4f is higher than previous: %.4f! Reverting and "
                          "redoing...\n" % (new_intercept, prev_intercept))
                    intercept = prev_intercept
                    new_intercept = int_calibration(name, directory, band, crop, sigma, given_catalog, chosen_survey,
                                                    mag_low_cutoff-0.5, mag_high_lim, grb_ra,
                                                    grb_dec, grb_coordlist, grb_thresh)
                    revert_flag = True
                    break
                else:
                    intercept = new_intercept
                    prev_intercept = intercept
                    revert_flag = False

            if revert_flag:
                print("Loop stopped due to intercept reverting to the previous value: %.4f" % intercept)
            else:
                print(f"Final intercept below 0.5: %.4f" % intercept)


def main():

    parser = argparse.ArgumentParser(
        description='runs sextractor and psfex on swarped img to get psf fit photometry, then '
                    'outputs ecsv w/ corrected mags')
    parser.add_argument('-exp_query', action='store_true', help='optional flag, exports ecsv of astroquery results '
                                                                'along with photometry')
    parser.add_argument('-exp_query_only', action='store_true',
                        help='optional flag, use if photom is run already to only generate query results')
    parser.add_argument('-no_plots', action='store_true',
                        help='optional flag, stops creation of mag comparison plot betw. PRIME and survey, '
                             'along with residual plot w/ statistics, lim mag plot')
    parser.add_argument('-keep', action='store_true',
                        help='optional flag, use if you DONT want to remove intermediate products after getting photom,'
                             ' i.e. the ".cat" and ".psf" files')
    parser.add_argument('-grb_only', action='store_true',
                        help='optional flag, use if running -grb again on already created catalog')
    parser.add_argument('-filepath', type=str, help='[str], full file path of stacked image, can also place just'
                                                    'filename and it will default to current directory')
    parser.add_argument('-band', type=str, help='[str], band, ex. "J"')
    parser.add_argument('-survey', type=str,
                        help='[str], *NOW OPTIONAL* manually specify which survey to query, choose from VHS, 2MASS'
                             ', VIKING, Skymapper, SDSS, UKIDSS, & DES.  If you leave out this arg, it will automatically'
                             'pick a survey from the above list depending on the area and coverage.',
                        default=None)
    parser.add_argument('-crop', type=int, help='[int], # of pixels from edge of image to crop, default = 300',
                        default=defaults["crop"])
    parser.add_argument('-sigma', type=float, help='[float], # of sigma w/ which to sigma clip for '
                                                  'zero point calculation, default = 3',
                        default=defaults["sigma"])
    parser.add_argument('-catalog', type=str, help='[str], optional field to supply an already generated'
                                                   'catalog for photometry INSTEAD of querying, put in full file path.',
                        default=None)
    parser.add_argument('-mag_low', type=float, help='[float], Lower mag cutoff for survey query & crossmatch'
                                                        ' settings default = 12.5',
                        default=None)
    parser.add_argument('-mag_high', type=float, help='[float], Higher mag cutoff for survey query & crossmatch'
                                                        ' settings default = 21, currently only applies to DES & Skymapper',
                        default=None)
    parser.add_argument('-grb_ra', type=str, help='[str], RA for GRB source, either in hh:mm:ss or decimal'
                                                  '*NOTE* When using sexagesimal, use "-grb_ra=value_here" NOT "-grb_ra '
                                                  'value_here", as argparse doesnt like negative sexagesimals',
                        default=defaults["RA"])
    parser.add_argument('-grb_dec', type=str, help='[str], DEC for GRB source, either in dd:mm:ss or decimal'
                                                   '*NOTE* When using sexagesimal, use "-grb_dec=value_here" NOT "-grb_dec '
                                                  'value_here", as argparse doesnt like negative sexagesimals',
                        default=defaults["DEC"])
    parser.add_argument('-grb_coordlist', type=str, nargs='+',
                        help='[float] Used to check multiple GRB locations.  Input RA and DECs of locations '
                             'with the format: -coordlist 123,45 -123,-45 etc..  *DONT USE -RA '
                             '& -DEC BUT INCLUDE -grb_radius*', default=None)
    parser.add_argument('-grb_radius', type=str,
                        help='[str], # of arcsec diameter to search for GRB, default = 4.0".  You can specify arcsec,'
                             ' arcmin, or deg w/ an underscore.  Ex. "-grb_radius 3_arcmin" will specify an area of 3 '
                             'arcminutes.  If just a number is applied, it defaults to arcsec.',
                        default=defaults["thresh"])
    parser.add_argument('-int_cal', action='store_true',
                        help='optional flag, use to automatically improve 3 sigma fit y-int.  When y-int is >0.5, the '
                             'low mag cutoff value is increased by 0.5, only stopping when y-int < 0.5.')
    args, unknown = parser.parse_known_args()
    # print(args)
    # print(unknown)

    photometry(args.filepath, args.band, args.crop, args.sigma, args.catalog, args.survey, args.mag_low,
               args.mag_high, args.no_plots, args.keep,
               args.grb_only, args.grb_ra, args.grb_dec, args.grb_coordlist, args.grb_radius, args.int_cal)


if __name__ == "__main__":
    main()
