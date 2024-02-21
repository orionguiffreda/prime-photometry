"""
Calibrates photometry for stacked image
"""

import numpy as np
import numpy.ma as ma
import warnings
import argparse
import sys
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.stats import sigma_clipped_stats, sigma_clip
import subprocess
from astropy.io import fits
from astropy.io import ascii
import os
from astropy.table import Column
import matplotlib.pyplot as plt
import statsmodels.api as sm
from scipy import stats

sys.path.insert(0,'C:\PycharmProjects\prime-photometry\photomitrus')
from settings import gen_config_file_name
#%%
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

#%%
#import img and get wcs
"""directory = '/Users/orion/Desktop/PRIME/GRB/preproc_swarp/'
imageName = 'coadd_cds.fits'"""

def img(directory,imageName,crop):
    os.chdir(directory)
    f = fits.open(imageName)
    data = f[0].data  #This is the image array
    header = f[0].header

    #strong the image WCS into an object
    w = WCS(header)

    #Get the RA and Dec of the center of the image
    [raImage, decImage] = w.all_pix2world(data.shape[0] / 2, data.shape[1] / 2, 1)

    #get ra and dec of right, left, top and bottom of img for box size
    [raREdge, decREdge] = w.all_pix2world(data.shape[0] - int(crop), data.shape[1] / 2, 1)
    [raLEdge, decLEdge] = w.all_pix2world(int(crop), data.shape[1] / 2, 1)
    [raTop, decTop] = w.all_pix2world(data.shape[0] / 2, data.shape[1] - int(crop), 1)
    [raBot, decBot] = w.all_pix2world(data.shape[0] / 2, int(crop), 1)

    height = (decTop - decBot) * 60
    width = (raLEdge - raREdge) * 60

    return data,header,w,raImage, decImage, width, height

#data,header,w,raImage, decImage,boxsize = img('/mnt/d/PRIME_photometry_test_files/GRB_Followup/GRB_Followup_flat_tests/flatstack/','coadd.Open-J.00888557-00888869.C4.fits',300)
#%%
# Use astroquery to get catalog search
from astroquery.vizier import Vizier
#Vizier.VIZIER_SERVER = 'vizier.ast.cam.ac.uk'

def query(raImage, decImage,filter, width, height, survey):
    # Use astroquery to get catalog search
    from astroquery.vizier import Vizier
    # Vizier.VIZIER_SERVER = 'vizier.ast.cam.ac.uk'
    if survey == '2MASS':
        catNum = 'II/246'  # changing to 2mass
        print('\nQuerying Vizier %s around RA %.4f, Dec %.4f with a box of width %.3f and height %.3f arcmin' % (
        catNum, raImage, decImage, width, height))
        try:
            # You can set the filters for the individual columns (magnitude range, number of detections) inside the Vizier query
            v = Vizier(columns=['*'], column_filters={"%smag" % filter: ">12", "Nd": ">6"}, row_limit=-1)
            Q = v.query_region(SkyCoord(ra=raImage, dec=decImage, unit=(u.deg, u.deg)), width=str(width) + 'm',
                               height=str(height) + 'm', catalog=catNum, cache=False)
            # query vizier around (ra, dec) with a radius of boxsize
            #print(Q[0])
            print('Queried source total = ', len(Q[0]))
        except:
            print('I cannnot reach the Vizier database. Is the internet working?')
    elif survey == 'VHS':
        catNum = 'II/367'  # changing to vista
        print('\nQuerying Vizier %s around RA %.4f, Dec %.4f with a box of width %.3f and height %.3f arcmin' % (
        catNum, raImage, decImage, width, height))
        try:
            # You can set the filters for the individual columns (magnitude range, number of detections) inside the Vizier query
            v = Vizier(columns=['RAJ2000','DEJ2000','%sap3' % filter,'e_%sap3' % filter], column_filters={"%sap3" % filter: ">12"}, row_limit=-1)
            Q = v.query_region(SkyCoord(ra=raImage, dec=decImage, unit=(u.deg, u.deg)), width=str(width) + 'm',
                               height=str(height) + 'm', catalog=catNum, cache=False)
            # query vizier around (ra, dec) with a radius of boxsize
            #print(Q[0])
            print('Queried source total = ', len(Q[0]))
        except:
            print(
                'Error in Vizier query. Perhaps your image is not in the southern hemisphere sky?')
    else:
        print('No supported survey found, currently use either 2MASS or VHS')
    return Q

def queryexport(Q, survey):
    table = Q[0]
    table.write('%s_query.ecsv' % survey, overwrite=True)
    print('%s_query.ecsv written!' % survey)

#Q = query(raImage, decImage, boxsize)
#table = Q[0]
#table.write('2MASS_query.ecsv', overwrite=True)
#%%
#run sextractor on swarped img to find sources
"""configFile = '/Users/orion/Desktop/PRIME/sex.config'
catalogName = imageName+'.cat'
paramName = '/Users/orion/Desktop/PRIME/tempsource.param'"""
def sex1(imageName):
    print('Running sextractor on img to initially find sources...')
    configFile = gen_config_file_name('sex2.config')
    paramName = gen_config_file_name('tempsource.param')
    catalogName = imageName + '.cat'
    try:
        command = 'sex %s -c %s -CATALOG_NAME %s -PARAMETERS_NAME %s' % (imageName, configFile, catalogName, paramName)
        #print('Executing command: %s' % command)
        rval = subprocess.run(command.split(), check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError as err:
        print('Could not run sextractor with exit error %s'%err)
    return catalogName

#%%
#run psfex on sextractor LDAC from previous step
def psfex(catalogName):
    print('Running PSFex on sextrctr catalogue to generate psf for stars in the img...')
    psfConfigFile = gen_config_file_name('default.psfex')
    try:
        command = 'psfex %s -c %s' % (catalogName,psfConfigFile)
        #print('Executing command: %s' % command)
        rval = subprocess.run(command.split(), check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError as err:
        print('Could not run psfex with exit error %s'%err)

#%%
#feed generated psf model back into sextractor w/ diff param (or you could use that param from the start but its slower)
"""psfName = imageName + '.psf'
psfcatalogName = imageName.replace('.fits','.psf.cat')
psfparamName = '/Users/orion/Desktop/PRIME/photomPSF.param'""" #This is a new set of parameters to be obtained from SExtractor, including PSF-fit magnitudes
def sex2(imageName):
    print('Feeding psf model back into sextractor for fitting and flux calculation...')
    psfName = imageName + '.psf'
    psfcatalogName = imageName.replace('.fits','.psf.cat')
    configFile = gen_config_file_name('sex2.config')
    psfparamName = gen_config_file_name('photomPSF.param')
    try:
        #We are supplying SExtactor with the PSF model with the PSF_NAME option
        command = 'sex %s -c %s -CATALOG_NAME %s -PSF_NAME %s -PARAMETERS_NAME %s' % (imageName, configFile, psfcatalogName, psfName, psfparamName)
        #print("Executing command: %s" % command)
        rval = subprocess.run(command.split(), check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError as err:
        print('Could not run sextractor with exit error %s'%err)
    return psfcatalogName

#%%
#read in tables, clean edges (potentially unnecessary), crossmatch
def tables(Q,psfcatalogName,crop):
    print('Cross-matching sextracted and catalog sources...')
    crop = int(crop)
    mass_imCoords = w.all_world2pix(Q[0]['RAJ2000'], Q[0]['DEJ2000'], 1)
    good_cat_stars = Q[0][np.where((mass_imCoords[0] > crop) & (mass_imCoords[0] < (max(mass_imCoords[0])-crop)) & (mass_imCoords[1] > crop) & (mass_imCoords[1] < (max(mass_imCoords[0])-crop)))]

    psfsourceTable = get_table_from_ldac(psfcatalogName)
    cleanPSFSources = psfsourceTable[(psfsourceTable['FLAGS']==0) & (psfsourceTable['FLAGS_MODEL']==0) & (psfsourceTable['XMODEL_IMAGE']
                < (max(psfsourceTable['XMODEL_IMAGE'])-crop)) & (psfsourceTable['XMODEL_IMAGE']>crop) &(psfsourceTable['YMODEL_IMAGE']<(max(psfsourceTable['YMODEL_IMAGE'])-crop))
                & (psfsourceTable['YMODEL_IMAGE']>crop)] #& (psfsourceTable['SNR_WIN']>3.0)]

    psfsourceCatCoords = SkyCoord(ra=cleanPSFSources['ALPHA_J2000'], dec=cleanPSFSources['DELTA_J2000'], frame='icrs', unit='degree')

    massCatCoords = SkyCoord(ra=good_cat_stars['RAJ2000'], dec=good_cat_stars['DEJ2000'], frame='icrs', unit='degree')
    #Now cross match sources
    #Set the cross-match distance threshold to 0.6 arcsec, or just about one pixel
    photoDistThresh = 1.0
    idx_psfimage, idx_psfmass, d2d, d3d = massCatCoords.search_around_sky(psfsourceCatCoords, photoDistThresh*u.arcsec)
    #idx_psfimage are indexes into psfsourceCatCoords for the matched sources, while idx_psfmass are indexes into massCatCoords for the matched sources

    print('Found %d good cross-matches'%len(idx_psfmass))
    return good_cat_stars,cleanPSFSources,idx_psfmass,idx_psfimage

#good_cat_stars,cleanPSFSources,idx_psfmass,idx_psfimage = tables(Q,psfcatalogName)
#%%
#derive zero pt / put in swarped header
def zeropt(good_cat_stars,cleanPSFSources,idx_psfmass,idx_psfimage,imageName,filter,survey):
    if survey == '2MASS':
        psfoffsets = ma.array(good_cat_stars['%smag' % filter][idx_psfmass] - cleanPSFSources['MAG_POINTSOURCE'][idx_psfimage])
        psfoffsets = psfoffsets.data
    elif survey == 'VHS':
        psfoffsets = ma.array(good_cat_stars['%sap3' % filter][idx_psfmass] - cleanPSFSources['MAG_POINTSOURCE'][idx_psfimage])
        psfoffsets = psfoffsets.data
    else:
        print('Surveys other than 2MASS or VHS currently not supported')
    #Compute sigma clipped statistics
    zero_psfmean, zero_psfmed, zero_psfstd = sigma_clipped_stats(psfoffsets)
    print('PSF Mean ZP: %.2f\nPSF Median ZP: %.2f\nPSF STD ZP: %.2f'%(zero_psfmean, zero_psfmed, zero_psfstd))

    #header.set('ZP_ERR', zero_psfstd, 'PSF Zero Point STD', after='SATURATE')
    #header.set('ZP_MEAN', zero_psfmean, 'PSF Zero Point Mean', after='SATURATE')
    #header.set('ZP_MED', zero_psfmed, 'PSF Zero Point Median', after='SATURATE')

    #print('Header values added = ', header['ZP_MEAN'], header['ZP_MED'], header['ZP_ERR'])

    #fits.writeto(imageName,data,header=header,overwrite=True)
    #print('Swarped img header rewritten w/ zero pts, all done!')

    psfmag = zero_psfmed + cleanPSFSources['MAG_POINTSOURCE']
    psfmagerr = np.sqrt(cleanPSFSources['MAGERR_POINTSOURCE'] ** 2 + zero_psfstd ** 2)

    psfmagcol = Column(psfmag, name = '%sMAG_PSF' % filter,unit='mag')
    psfmagerrcol = Column(psfmagerr, name = 'e_%sMAG_PSF' % filter,unit='mag')
    cleanPSFSources.add_column(psfmagcol)
    cleanPSFSources.add_column(psfmagerrcol)
    cleanPSFSources.remove_column('VIGNET')

    cleanPSFSources.write('%s.%s.ecsv' % (imageName,survey),overwrite=True)
    print('%s.%s.ecsv written, CSV w/ corrected mags' % (imageName,survey))

#%% optional plots
def plots(directory,imageName,survey,filter,good_cat_stars,idx_psfmass,idx_psfimage):
    #appropriate mag column
    if survey == '2MASS':
        magcol = '%smag' % filter
    elif survey == 'VHS':
        magcol = '%sap3' % filter
    else:
        print('Surveys other than 2MASS or VHS currently not supported')

    # mag comparison plot
    PRIMEdata = ascii.read(directory + imageName + '.' + survey + '.ecsv')
    plt.figure(1,figsize=(8, 8))
    plt.plot(PRIMEdata['%sMAG_PSF' % filter][idx_psfimage],good_cat_stars['%s' % magcol][idx_psfmass],
             'r.', markersize=14, markeredgecolor='black',)
    plt.xlim(10, 22)
    plt.ylim(10, 22)
    plt.title('PRIME Mags vs %s Mags' % survey)
    plt.xlabel('PRIME Mags', fontsize=15)
    plt.ylabel('%s Mags' % survey, fontsize=15)
    plt.grid()
    plt.savefig('%s_mag_comp_plot.png' % survey)
    print('Saved mag comparison plot to dir!')

    #residual fit
    x = PRIMEdata['%sMAG_PSF' % filter][idx_psfimage]
    y = good_cat_stars['%s' % magcol][idx_psfmass]
    x_const = sm.add_constant(x)
    model = sm.OLS(y, x_const).fit()
    m = model.params[1]
    merr = model.bse[1]
    b = model.params[0]
    berr = model.bse[0]
    rsquare = model.rsquared
    rss = model.ssr
    res_err = np.std(model.resid)

    #residual plot
    plt.figure(2, figsize=(8, 6))
    plt.scatter(PRIMEdata['%sMAG_PSF' % filter][idx_psfimage], model.resid, color='red')
    plt.ylim(-1.5, 1.5)
    plt.title('PRIME vs %s Residuals' % survey)
    plt.ylabel('Residuals')
    plt.xlabel('Mags')
    plt.axhline(y=res_err, color='tab:orange', linestyle='--', linewidth=1)
    plt.axhline(y=res_err * 2, color='green', linestyle='--', linewidth=1)
    plt.axhline(y=-res_err, color='tab:orange', linestyle='--', linewidth=1)
    plt.axhline(y=-res_err * 2, color='green', linestyle='--', linewidth=1)
    plt.axhline(y=0, color='black', linestyle='--', linewidth=1)
    plt.legend(['Residuals', r'1 $\sigma$ = %.3f' % res_err, r'2 $\sigma$'], loc='lower left')
    info = ('slope = %.3f +/- %.3f' % (m,merr))+('\nintercept = %.3f +/- %.3f' % (b,berr))+('\nR$^{2}$ = %.3f' % rsquare)+('\nRSS = %d' % rss)
    plt.text(15,-1.25,info,bbox=dict(facecolor='none',edgecolor='black',pad=5.0))
    plt.savefig('%s_residual_plot.png' % survey,dpi=300)
    print('Saved residual plot to dir!')


#%%
defaults = dict(crop=300,survey='2MASS')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='runs sextractor and psfex on swarped img to get psf fit photometry, then '
                                                 'outputs ecsv w/ corrected mags')
    parser.add_argument('-exp_query', action='store_true', help='optional flag, exports ecsv of astroquery results '
                                                                'along with photometry')
    parser.add_argument('-plots', action='store_true', help='optional flag, exports mag comparison plot betw. PRIME and survey, '
                                                            'along with residual plot w/ statistics')
    parser.add_argument('-dir', type=str, help='[str], directory of stacked image')
    parser.add_argument('-name', type=str, help='[str], image name')
    parser.add_argument('-filter', type=str, help='[str], filter, ex. "J"')
    parser.add_argument('-survey', type=str, help='[str], survey to query, currently use either "2MASS" '
                                                  'or "VHS" (for VISTA in s. hemi. only), default = 2MASS')
    parser.add_argument('-crop', type=str, help='[int], # of pixels from edge of image to crop, default = 300',default=defaults["crop"])
    args = parser.parse_args()

    if args.exp_query and args.plots:
        data, header, w, raImage, decImage, width, height = img(args.dir, args.name, args.crop)
        Q = query(raImage, decImage, args.filter, width, height, args.survey)
        queryexport(Q,args.survey)
        catalogName = sex1(args.name)
        psfex(catalogName)
        psfcatalogName = sex2(args.name)
        good_cat_stars, cleanPSFSources, idx_psfmass, idx_psfimage = tables(Q, psfcatalogName, args.crop)
        zeropt(good_cat_stars, cleanPSFSources, idx_psfmass, idx_psfimage, args.name, args.filter, args.survey)
        plots(args.dir, args.name, args.survey, args.filter, good_cat_stars, idx_psfmass, idx_psfimage)
    elif args.exp_query:
        data, header, w, raImage, decImage, width, height = img(args.dir, args.name, args.crop)
        Q = query(raImage, decImage, args.filter, width, height, args.survey)
        queryexport(Q,args.survey)
        catalogName = sex1(args.name)
        psfex(catalogName)
        psfcatalogName = sex2(args.name)
        good_cat_stars, cleanPSFSources, idx_psfmass, idx_psfimage = tables(Q, psfcatalogName, args.crop)
        zeropt(good_cat_stars, cleanPSFSources, idx_psfmass, idx_psfimage, args.name, args.filter, args.survey)
    elif args.plots:
        data, header, w, raImage, decImage, width, height = img(args.dir, args.name, args.crop)
        Q = query(raImage, decImage, args.filter, width, height, args.survey)
        catalogName = sex1(args.name)
        psfex(catalogName)
        psfcatalogName = sex2(args.name)
        good_cat_stars, cleanPSFSources, idx_psfmass, idx_psfimage = tables(Q, psfcatalogName, args.crop)
        zeropt(good_cat_stars, cleanPSFSources, idx_psfmass, idx_psfimage, args.name, args.filter,args.survey)
        plots(args.dir,args.name,args.survey,args.filter,good_cat_stars,idx_psfmass,idx_psfimage)
    else:
        data, header, w, raImage, decImage, width, height = img(args.dir, args.name, args.crop)
        Q = query(raImage, decImage, args.filter, width, height, args.survey)
        catalogName = sex1(args.name)
        psfex(catalogName)
        psfcatalogName = sex2(args.name)
        good_cat_stars, cleanPSFSources, idx_psfmass, idx_psfimage = tables(Q, psfcatalogName, args.crop)
        zeropt(good_cat_stars, cleanPSFSources, idx_psfmass, idx_psfimage, args.name, args.filter,args.survey)

#%%

"""dir = '/Users/orion/Desktop/PRIME/GRB/preproc_swarp/'
imageName = 'coadd_cds.fits'
configFile = '/Users/orion/Desktop/PRIME/sex.config'
catalogName = imageName+'.cat'
paramName = '/Users/orion/Desktop/PRIME/tempsource.param'
psfName = imageName + '.psf'
psfcatalogName = imageName.replace('.fits','.psf.cat')
psfparamName = '/Users/orion/Desktop/PRIME/photomPSF.param'"""
