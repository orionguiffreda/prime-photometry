"""
Calibrates photometry for stacked image
"""

import numpy as np
import numpy.ma as ma
import argparse
import sys
import astropy.units as u
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.wcs import utils
# from astropy.stats import sigma_clipped_stats
from astropy.stats import sigma_clip
from astropy.io import fits
from astropy.io import ascii
import os
from astropy.table import Column
from astropy.table import Table
import matplotlib.pyplot as plt
import statsmodels.api as sm
import subprocess

# from scipy import stats

sys.path.insert(0, 'C:\PycharmProjects\prime-photometry\photomitrus')
from photomitrus.settings import gen_config_file_name

# %%
defaults = dict(crop=300, survey='2MASS', RA=None, DEC=None, thresh=4.0)


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

    return data, header, w, raImage, decImage, width, height


# %%
# Use astroquery to get catalog search


def query(raImage, decImage, band, width, height, survey):
    # Use astroquery to get catalog search
    # Vizier.VIZIER_SERVER = 'vizier.ast.cam.ac.uk'
    if survey == '2MASS':
        catNum = 'II/246'  # changing to 2mass
        print('\nQuerying Vizier %s around RA %.4f, Dec %.4f with a box of width %.3f and height %.3f arcmin' % (
            catNum, raImage, decImage, (width), (height)))
        try:
            # You can set the bands for the individual columns (magnitude range, number of detections) inside the Vizier query
            v = Vizier(columns=['RAJ2000', 'DEJ2000', '%smag' % band, 'e_%smag' % band],
                       column_filters={"%smag" % band: ">12.5", "Nd": ">6"}, row_limit=-1)
            Q = v.query_region(SkyCoord(ra=raImage, dec=decImage, unit=(u.deg, u.deg)), width=str(width) + 'm',
                               height=str(height) + 'm', catalog=catNum, cache=False)
            # query vizier around (ra, dec) with a radius of boxsize
            # print(Q[0])
            print('Queried source total = ', len(Q[0]))
        except:
            print('I cannnot reach the Vizier database. Is the internet working?')
    elif survey == 'VHS':
        catNum = 'II/367'  # changing to vista
        print('\nQuerying Vizier %s around RA %.4f, Dec %.4f with a box of width %.3f and height %.3f arcmin' % (
            catNum, raImage, decImage, width, height))
        try:
            # You can set the bands for the individual columns (magnitude range, number of detections) inside the Vizier query
            v = Vizier(columns=['RAJ2000', 'DEJ2000', '%sap3' % band, 'e_%sap3' % band],
                       column_filters={"%sap3" % band: ">12.5"}, row_limit=-1)
            Q = v.query_region(SkyCoord(ra=raImage, dec=decImage, unit=(u.deg, u.deg)), width=str(width) + 'm',
                               height=str(height) + 'm', catalog=catNum, cache=False)
            # query vizier around (ra, dec) with a radius of boxsize
            # print(Q[0])
            print('Queried source total = ', len(Q[0]))
        except:
            print(
                'Error in Vizier query. Perhaps your image is not in the southern hemisphere sky?  H band is also not well covered!')
    elif survey == 'VIKING':
        # catNum = 'II/348/vvv2'  # changing to vista
        catNum = 'II/343/viking2'
        print('\nQuerying Vizier %s around RA %.4f, Dec %.4f with a box of width %.3f and height %.3f arcmin' % (
            catNum, raImage, decImage, width, height))
        try:
            # You can set the bands for the individual columns (magnitude range, number of detections) inside the Vizier query
            v = Vizier(columns=['RAJ2000', 'DEJ2000', '%sap3' % band, 'e_%sap3' % band],
                       column_filters={"%sap3" % band: ">12.5"}, row_limit=-1)
            print(v)
            Q = v.query_region(SkyCoord(ra=raImage, dec=decImage, unit=(u.deg, u.deg)), width=str(width) + 'm',
                               height=str(height) + 'm', catalog=catNum, cache=False)
            # query vizier around (ra, dec) with a radius of boxsize
            # print(Q[0])
            print('Queried source total = ', len(Q[0]))
        except:
            print(
                'Error in Vizier query. Perhaps your image is not in the southern hemisphere sky?  H band is also not well covered!'
                ' If you are in S.H., VIKING is only in a relatively smaller strip!')
    elif survey == 'Skymapper':
        # catNum = 'II/348/vvv2'  # changing to vista
        catNum = 'II/379/smssdr4'
        print('\nQuerying Vizier %s around RA %.4f, Dec %.4f with a box of width %.3f and height %.3f arcmin' % (
            catNum, raImage, decImage, width, height))
        try:
            # You can set the bands for the individual columns (magnitude range, number of detections) inside the Vizier query
            v = Vizier(columns=['RAICRS', 'DEICRS', '%sPSF' % band.lower(), 'e_%sPSF' % band.lower()],
                       column_filters={"%sPSF" % band.lower(): ">12.5"}, row_limit=-1)
            print(v)
            Q = v.query_region(SkyCoord(ra=raImage, dec=decImage, unit=(u.deg, u.deg)), width=str(width) + 'm',
                               height=str(height) + 'm', catalog=catNum, cache=False)
            # query vizier around (ra, dec) with a radius of boxsize
            # print(Q[0])
            print('Queried source total = ', len(Q[0]))
        except:
            print(
                'Error in Vizier query. Perhaps your image is not in the southern hemisphere sky?'
                '\n perhaps check skymapper coverage maps?')
    elif survey == 'SDSS':
        # catNum = 'II/348/vvv2'  # changing to vista
        catNum = 'V/154/sdss16'
        print('\nQuerying Vizier %s around RA %.4f, Dec %.4f with a box of width %.3f and height %.3f arcmin' % (
            catNum, raImage, decImage, width, height))
        try:
            # You can set the bands for the individual columns (magnitude range, number of detections) inside the Vizier query
            v = Vizier(columns=['RA_ICRS', 'DE_ICRS', '%spmag' % band.lower(), 'e_%spmag' % band.lower()],
                       column_filters={"%spmag" % band.lower(): ">12.5"
                           , "clean": "=1"}, row_limit=-1)
            print(v)
            Q = v.query_region(SkyCoord(ra=raImage, dec=decImage, unit=(u.deg, u.deg)), width=str(width) + 'm',
                               height=str(height) + 'm', catalog=catNum, cache=False)
            # query vizier around (ra, dec) with a radius of boxsize
            # print(Q[0])
            print('Queried source total = ', len(Q[0]))
        except:
            print(
                'Error in Vizier query. Perhaps your image is not in the southern hemisphere sky?'
                '\n perhaps check SDSS coverage maps?')
    elif survey == 'UKIDSS':
        # catNum = 'II/348/vvv2'  # changing to vista
        catNum = 'II/319/las9'
        print('\nQuerying Vizier %s around RA %.4f, Dec %.4f with a box of width %.3f and height %.3f arcmin' % (
            catNum, raImage, decImage, width, height))
        try:
            # You can set the bands for the individual columns (magnitude range, number of detections) inside the Vizier query
            v = Vizier(columns=['RAJ2000', 'DEJ2000', '%smag' % band, 'e_%smag' % band],
                       column_filters={"%smag" % band: ">12.5"}, row_limit=-1)
            print(v)
            Q = v.query_region(SkyCoord(ra=raImage, dec=decImage, unit=(u.deg, u.deg)), width=str(width) + 'm',
                               height=str(height) + 'm', catalog=catNum, cache=False)
            # query vizier around (ra, dec) with a radius of boxsize
            # print(Q[0])
            print('Queried source total = ', len(Q[0]))
        except:
            print(
                'Error in Vizier query. Perhaps your image is not in the southern hemisphere sky?'
                '\n perhaps check UKIDSS coverage maps?')
    else:
        print('No supported survey found, currently use either 2MASS, VHS, VIKING, Skymapper, SDSS, or UKIDSS')
    return Q


# %%
# run sextractor on swarped img to find sources


def sex1(imageName):
    print('Running sextractor on img to initially find sources...')
    configFile = gen_config_file_name('sex2.config')
    paramName = gen_config_file_name('tempsource.param')
    catalogName = imageName + '.cat'
    try:
        command = 'sex %s -c %s -CATALOG_NAME %s -PARAMETERS_NAME %s' % (imageName, configFile, catalogName, paramName)
        # print('Executing command: %s' % command)
        rval = subprocess.run(command.split(), check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
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
        rval = subprocess.run(command.split(), check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
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
    try:
        # We are supplying SExtactor with the PSF model with the PSF_NAME option
        command = 'sex %s -c %s -CATALOG_NAME %s -PSF_NAME %s -PARAMETERS_NAME %s' % (
        imageName, configFile, psfcatalogName, psfName, psfparamName)
        # print("Executing command: %s" % command)
        rval = subprocess.run(command.split(), check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError as err:
        print('Could not run sextractor with exit error %s' % err)
    return psfcatalogName


# %%
# read in tables, crossmatch
def tables(Q, data, w, psfcatalogName, crop):
    print('Cross-matching sextracted and catalog sources...')
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
    print('Catalogue cropped, source total = ', len(good_cat_stars))

    psfsourceTable = get_table_from_ldac(psfcatalogName)

    PSFSources = psfsourceTable[
        (psfsourceTable['XMODEL_IMAGE'] < (max_x - crop)) & (psfsourceTable['XMODEL_IMAGE'] > crop)
        & (psfsourceTable['YMODEL_IMAGE'] < (max_y) - crop) & (psfsourceTable['YMODEL_IMAGE'] > crop)]

    cleanPSFSources = psfsourceTable[
        (psfsourceTable['FLAGS'] == 0) & (psfsourceTable['FLAGS_MODEL'] == 0) & (psfsourceTable['XMODEL_IMAGE']
                                                                                 < (max_x - crop)) & (
                    psfsourceTable['XMODEL_IMAGE'] > crop) & (psfsourceTable['YMODEL_IMAGE'] < (max_y) - crop)
        & (psfsourceTable['YMODEL_IMAGE'] > crop)]

    # magmean = cleanPSFSources['MAG_POINTSOURCE'].mean()  #trimming sources more than 4 sig away from mean (outliers)
    # magerr = cleanPSFSources['MAG_POINTSOURCE'].std()
    # maglow = magmean - 4 * magerr
    # maghigh = magmean + 4 * magerr

    # cleanPSFSources = cleanPSFSources[(cleanPSFSources['MAG_POINTSOURCE'] >= maglow) & (cleanPSFSources['MAG_POINTSOURCE'] <= maghigh)]

    # psfsourceCatCoords = SkyCoord(ra=cleanPSFSources['ALPHA_J2000'], dec=cleanPSFSources['DELTA_J2000'], frame='icrs', unit='degree')

    psfsourceCatCoords = utils.pixel_to_skycoord(cleanPSFSources['X_IMAGE'], cleanPSFSources['Y_IMAGE'], w, origin=1)

    massCatCoords = SkyCoord(ra=good_cat_stars[RA], dec=good_cat_stars[DEC], frame='icrs', unit='degree')
    # Now cross match sources
    # Set the cross-match distance threshold to 0.6 arcsec, or just about one pixel
    photoDistThresh = 0.6
    idx_psfimage, idx_psfmass, d2d, d3d = massCatCoords.search_around_sky(psfsourceCatCoords,
                                                                          photoDistThresh * u.arcsec)
    # idx_psfimage are indexes into psfsourceCatCoords for the matched sources, while idx_psfmass are indexes into massCatCoords for the matched sources

    print('Found %d good cross-matches' % len(idx_psfmass))
    return good_cat_stars, cleanPSFSources, PSFSources, idx_psfmass, idx_psfimage


def queryexport(good_cat_stars, imageName, survey):
    table = good_cat_stars
    chip = imageName[-6]
    num = imageName[-16:-8]
    table.write('%s_C%s_query_%s.ecsv' % (survey, chip, num), overwrite=True)
    print('%s_C%s_query_%s.ecsv written!' % (survey, chip, num))


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


def zeropt(good_cat_stars, cleanPSFSources, PSFSources, idx_psfmass, idx_psfimage, imageName, band, survey):
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
    psf_clipped = sigma_clip(psfoffsets)
    psfoffsets = psfoffsets[~psf_clipped.mask]
    psfweights = psfweights_noclip[~psf_clipped.mask]

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
    print('Total source # = ', len(PSFSources))

    PSFSources.write('%s.%s.ecsv' % (imageName, survey), overwrite=True)
    print('%s.%s.ecsv written, CSV w/ corrected mags' % (imageName, survey))

    return cleanPSFSources, PSFSources, psfweights_noclip, psf_clipped


# %% optional GRB-specific photom


def GRB(ra, dec, imageName, survey, band, thresh, coordlist=None):
    mag_ecsvname = '%s.%s.ecsv' % (imageName, survey)
    mag_ecsvtable = ascii.read(mag_ecsvname)
    mag_ecsvcleanSources = mag_ecsvtable  # [(mag_ecsvtable['FLAGS'] == 0) & (mag_ecsvtable['FLAGS_MODEL'] == 0)]
    mag_ecsvsourceCatCoords = SkyCoord(ra=mag_ecsvcleanSources['ALPHA_J2000'], dec=mag_ecsvcleanSources['DELTA_J2000'],
                                       frame='icrs',
                                       unit='degree')

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
                print('GRB source at inputted coords %s and %s found!' % (ra, dec))

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
                print('Multiple sources detected in search radius, refer to .ecsv file for source info!')
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
            print('GRB source at inputted coords %s and %s found!' % (ra, dec))

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
            print('Multiple sources detected in search radius, refer to .ecsv file for source info!')
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


def photometry_plots(cleanPSFsources, PSFsources, imageName, survey, band, good_cat_stars, idx_psfmass, idx_psfimage,
                     psfweights_noclip, psf_clipped):
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

    avgsig = np.average(model_sig.resid, weights=psfweights_noclip[~psf_clipped.mask])
    varsig = np.average((model_sig.resid - avgsig)**2, weights=psfweights_noclip[~psf_clipped.mask])
    ressig_err = np.sqrt(varsig)

    # t = Table()
    # t['Residuals'] = model_sig.resid
    # t['Jmags'] = x_sig
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
    plt.figure(3, figsize=(8, 6))
    plt.scatter(cleanPSFsources['%sMAG_PSF' % band][idx_psfimage][~psf_clipped.mask], model_sig.resid, color='red')
    plt.ylim(-1.5, 1.5)
    plt.xlim(10, 21)
    plt.title('PRIME vs %s Residuals' % survey)
    plt.ylabel('Residuals')
    plt.xlabel('%s Mags' % band)
    plt.axhline(y=ressig_err, color='tab:orange', linestyle='--', linewidth=1)
    plt.axhline(y=ressig_err * 2, color='green', linestyle='--', linewidth=1)
    plt.axhline(y=-ressig_err, color='tab:orange', linestyle='--', linewidth=1)
    plt.axhline(y=-ressig_err * 2, color='green', linestyle='--', linewidth=1)
    plt.axhline(y=0, color='black', linestyle='--', linewidth=1)
    plt.legend([r'1 $\sigma$ = %.3f' % ressig_err, r'2 $\sigma$ = %.3f' % (2 * ressig_err)], loc='lower left')
    info2 = ('eqn: y = mx+b' + '\nslope = %.4f +/- %.4f' % (m_sig, m_sigerr)) + (
                '\nintercept = %.3f +/- %.3f' % (b_sig, b_sigerr)) + ('\nR$^{2}$ = %.3f' % rsquare_sig) + ('\nRSS = %d' % rss_sig)
    plt.text(15, -1.25, info2, bbox=dict(facecolor='white', edgecolor='black', pad=5.0))
    plt.savefig('%s_C%s_residual_plot_int_%s.png' % (survey, chip, num), dpi=300)

    # res plot y int, histogram
    bin_num_int = round(len(idx_psfimage) / 50)

    plt.figure(8, figsize=(10, 6))
    plt.hist2d(x=cleanPSFsources['%sMAG_PSF' % band][idx_psfimage][~psf_clipped.mask], y=model_sig.resid,
               bins=[bin_num_int, bin_num_int], range=[[10, 21],[-1, 1]], cmap='gist_heat_r')
    plt.colorbar(label='Density')
    plt.ylim(-1, 1)
    plt.xlim(10, 21)
    plt.title('PRIME vs %s Residuals w/ 3 Sigma Clip- Density Histogram' % survey)
    plt.ylabel('Residuals')
    plt.xlabel('%s Mags' % band)
    plt.axhline(y=ressig_err, color='tab:orange', linestyle='--', linewidth=1)
    plt.axhline(y=ressig_err * 2, color='green', linestyle='--', linewidth=1)
    plt.axhline(y=-ressig_err, color='tab:orange', linestyle='--', linewidth=1)
    plt.axhline(y=-ressig_err * 2, color='green', linestyle='--', linewidth=1)
    plt.axhline(y=0, color='black', linestyle='--', linewidth=1)
    plt.legend(['Residuals', r'1 $\sigma$ = %.3f' % ressig_err, r'2 $\sigma$ = %.3f' % (2 * ressig_err)], loc='lower left')
    infohist = ('eqn: y = mx+b' + '\nslope = %.4f +/- %.4f' % (m_sig, m_sigerr)) + (
                '\nintercept = %.3f +/- %.3f' % (b_sig, b_sigerr)) + ('\nR$^{2}$ = %.3f' % rsquare_sig) + ('\nRSS = %d' % rss_sig)
    plt.text(15, -0.9, infohist, fontsize=9, bbox=dict(facecolor='white', edgecolor='black', pad=5.0))
    plt.savefig('%s_C%s_residual_plot_int_hist_%s.png' % (survey, chip, num), dpi=300)

    print('Saved y-int residual plots to dir!')

    # WLS fit line over data plot
    def predict_y_for(x):
        return m2 * x + b2

    txt = ('slope = %.4f' % m2 + '\nslope err = %.4f' % m2err + '\nint = %.4f' % b2 + '\nint err = %.4f' % b2err)

    plt.figure(4, figsize=(8, 8))
    plt.xlim(10, 22)
    plt.ylim(10, 22)
    plt.title('PRIME vs %s w/ Weighted Fit' % survey)
    plt.grid()
    plt.xlabel('PRIME %s Mags' % band, fontsize=15)
    plt.ylabel('%s %s Mags' % (survey, band), fontsize=15)
    plt.scatter(cleanPSFsources['%sMAG_PSF' % band][idx_psfimage], good_cat_stars['%s' % magcol][idx_psfmass])
    plt.plot(cleanPSFsources['%sMAG_PSF' % band][idx_psfimage],
             predict_y_for(cleanPSFsources['%sMAG_PSF' % band][idx_psfimage]), c='r')
    box = dict(facecolor='white')
    plt.text(11, 18, txt, fontsize=12, bbox=box)
    plt.savefig('%s_C%s_WLS_fit_plot_%s.png' % (survey, chip, num), dpi=300)

    # WLS hist density plot
    bin_num = round(len(idx_psfimage) / 20)

    plt.figure(5, figsize=(10, 8))
    plt.xlim(10, 22)
    plt.ylim(10, 22)
    plt.title('PRIME vs %s w/ Weighted Fit - Density Histogram' % survey)
    plt.grid()
    plt.xlabel('PRIME %s Mags' % band, fontsize=15)
    plt.ylabel('%s %s Mags' % (survey, band), fontsize=15)
    plt.hist2d(x=cleanPSFsources['%sMAG_PSF' % band][idx_psfimage], y=good_cat_stars['%s' % magcol][idx_psfmass],
               bins=[bin_num, bin_num], range=[[10, 22],[10, 22]], cmap='gist_heat_r')
    plt.plot(cleanPSFsources['%sMAG_PSF' % band][idx_psfimage],
             predict_y_for(cleanPSFsources['%sMAG_PSF' % band][idx_psfimage]), c='r')
    plt.colorbar(label='Density')
    box = dict(facecolor='white')
    plt.text(11, 18, txt, fontsize=12, bbox=box)
    plt.savefig('%s_C%s_WLS_fit_hist_plot_%s.png' % (survey, chip, num), dpi=300)

    print('Saved WLS fit plots to dir!')

    # WLS fit for 3 sig clip of data

    sigtxt = ('slope = %.4f' % m_sig + '\nslope err = %.4f' % m_sigerr + '\nint = %.4f' % b_sig +
              '\nint err = %.4f' % b_sigerr)

    plt.figure(6, figsize=(8, 8))
    plt.xlim(10, 22)
    plt.ylim(10, 22)
    plt.title('PRIME vs %s w/ Weighted Fit - 3 Sigma Clip' % survey)
    plt.grid()
    plt.xlabel('PRIME %s Mags' % band, fontsize=15)
    plt.ylabel('%s %s Mags' % (survey, band), fontsize=15)
    plt.scatter(cleanPSFsources['%sMAG_PSF' % band][idx_psfimage][~psf_clipped.mask],
                good_cat_stars['%s' % magcol][idx_psfmass][~psf_clipped.mask])
    plt.plot(cleanPSFsources['%sMAG_PSF' % band][idx_psfimage][~psf_clipped.mask],
             predict_y_for(cleanPSFsources['%sMAG_PSF' % band][idx_psfimage][~psf_clipped.mask]), c='r')
    box = dict(facecolor='white')
    plt.text(11, 18, sigtxt, fontsize=12, bbox=box)
    plt.savefig('%s_C%s_WLS_fit_3sig_plot%s.png' % (survey, chip, num), dpi=300)

    # WLS 3 sig hist density plot
    plt.figure(7, figsize=(10, 8))
    plt.xlim(10, 22)
    plt.ylim(10, 22)
    plt.title('PRIME vs %s w/ Weighted Fit - 3 Sigma Clip - Density Histogram' % survey)
    plt.grid()
    plt.xlabel('PRIME %s Mags' % band, fontsize=15)
    plt.ylabel('%s %s Mags' % (survey, band), fontsize=15)
    plt.hist2d(x=cleanPSFsources['%sMAG_PSF' % band][idx_psfimage][~psf_clipped.mask],
               y=good_cat_stars['%s' % magcol][idx_psfmass][~psf_clipped.mask],
               bins=[bin_num, bin_num], range=[[10, 22],[10, 22]], cmap='gist_heat_r')
    plt.plot(cleanPSFsources['%sMAG_PSF' % band][idx_psfimage][~psf_clipped.mask],
             predict_y_for(cleanPSFsources['%sMAG_PSF' % band][idx_psfimage][~psf_clipped.mask]), c='r')
    plt.colorbar(label='Density')
    box = dict(facecolor='white')
    plt.text(11, 18, sigtxt, fontsize=12, bbox=box)
    plt.savefig('%s_C%s_WLS_fit_3sig_hist_plot%s.png' % (survey, chip, num), dpi=300)

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

    plt.figure(figsize=(24, 8))
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


def photometry(
        full_filename, band, survey=defaults['survey'], crop=defaults['crop'], no_plots=False, plots_only=False,
        keep=False, grb_only=False, grb_ra=None, grb_dec=None, grb_coordlist=None,
        grb_radius=defaults['thresh']
):
    directory = os.path.dirname(full_filename)
    if directory == '':
        directory = '.'
    directory = directory + '/'
    name = os.path.basename(full_filename)

    if plots_only:
        psfcatalogName = []
        for f in os.listdir(directory):
            if f.endswith('.psf.cat'):
                psfcatalogName.append(f)
        psfcatalogName = ' '.join(psfcatalogName)
        data, header, w, raImage, decImage, width, height = img(directory, name, crop)
        Q = query(raImage, decImage, band, width, height, survey)
        good_cat_stars, cleanPSFSources, PSFsources, idx_psfmass, idx_psfimage = tables(Q, data, w, psfcatalogName,
                                                                                        crop)
        cleanPSFSources, PSFsources, psfweights_noclip, psf_clipped = zeropt(good_cat_stars, cleanPSFSources, PSFsources, idx_psfmass, idx_psfimage,
                                             name, band, survey)
        photometry_plots(cleanPSFSources, PSFsources, name, survey, band, good_cat_stars, idx_psfmass,
                         idx_psfimage, psfweights_noclip, psf_clipped)
    elif grb_only:
        os.chdir(directory)
        if grb_coordlist:
            GRB(grb_ra, grb_dec, name, survey, band, grb_radius, grb_coordlist)
        else:
            GRB(grb_ra, grb_dec, name, survey, band, grb_radius)
    else:
        data, header, w, raImage, decImage, width, height = img(directory, name, crop)
        Q = query(raImage, decImage, band, width, height, survey)
        catalogName = sex1(name)
        psfex(catalogName)
        psfcatalogName = sex2(name)
        good_cat_stars, cleanPSFSources, PSFSources, idx_psfmass, idx_psfimage = tables(Q, data, w, psfcatalogName,
                                                                                        crop)
        cleanPSFSources, PSFsources, psfweights_noclip, psf_clipped = zeropt(good_cat_stars, cleanPSFSources, PSFSources, idx_psfmass, idx_psfimage,
                                             name, band, survey)
        if grb_ra:
            GRB(grb_ra, grb_dec, name, survey, band, grb_radius)
        elif grb_coordlist:
            GRB(grb_ra, grb_dec, name, survey, band, grb_radius, grb_coordlist)
        if not no_plots:
            photometry_plots(cleanPSFSources, PSFsources, name, survey, band, good_cat_stars, idx_psfmass,
                             idx_psfimage, psfweights_noclip, psf_clipped)
        if not keep:
            removal(directory)


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
    parser.add_argument('-plots_only', action='store_true',
                        help='optional flag, use if photom is run already to only generate plots quicker *CURRENTLY '
                             'BROKEN DO NOT USE*')
    parser.add_argument('-keep', action='store_true',
                        help='optional flag, use if you DONT want to remove intermediate products after getting photom,'
                             ' i.e. the ".cat" and ".psf" files')
    parser.add_argument('-grb_only', action='store_true',
                        help='optional flag, use if running -grb again on already created catalog')
    parser.add_argument('-filepath', type=str, help='[str], full file path of stacked image, can also place just'
                                                    'filename and it will default to current directory')
    parser.add_argument('-band', type=str, help='[str], band, ex. "J"')
    parser.add_argument('-survey', type=str,
                        help='[str], survey to query, choose from VHS (best / reliable for J), 2MASS (reliable for H)'
                             ', VIKING, Skymapper (reliable for Z), SDSS, and UKIDSS, default = 2MASS',
                        default=defaults["survey"])
    parser.add_argument('-crop', type=int, help='[int], # of pixels from edge of image to crop, default = 300',
                        default=defaults["crop"])
    parser.add_argument('-grb_ra', type=float, help='[float], RA for GRB source',
                        default=defaults["RA"])
    parser.add_argument('-grb_dec', type=float, help='[float], DEC for GRB source',
                        default=defaults["DEC"])
    parser.add_argument('-grb_coordlist', type=str, nargs='+',
                        help='[float] Used to check multiple GRB locations.  Input RA and DECs of locations '
                             'with the format: -coordlist 123,45 -123,-45 etc..  *DONT USE -RA '
                             '& -DEC BUT INCLUDE -grb_radius*', default=None)
    parser.add_argument('-grb_radius', type=float,
                        help='[float], # of arcsec diameter to search for GRB, default = 4.0"',
                        default=defaults["thresh"])
    args, unknown = parser.parse_known_args()

    photometry(args.filepath, args.band, args.survey, args.crop, args.no_plots, args.plots_only, args.keep,
               args.grb_only, args.grb_ra, args.grb_dec, args.grb_coordlist, args.grb_radius)


# %%


if __name__ == "__main__":
    main()
