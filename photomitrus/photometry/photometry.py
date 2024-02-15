"""
Calibrates photometry for stacked image
"""

import numpy as np
import numpy.ma as ma
import argparse
import sys
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.stats import sigma_clipped_stats, sigma_clip
import subprocess
from astropy.io import fits
import os
from astropy.table import Column

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
    raImage, decImage = w.all_pix2world(data.shape[0]/2, data.shape[1]/2, 1)

    #get ra and dec of edge, subtract to get rough rad of image
    [raEdge, decEdge] = w.all_pix2world(data.shape[0]-int(crop), data.shape[1] / 2, 1)
    boxsize = (raImage - raEdge)*60 #rad in arcmin

    return data,header,w,raImage, decImage, boxsize

#data,header,w,raImage, decImage = img(directory,imageName)
#%%
# Use astroquery to get catalog search
from astroquery.vizier import Vizier
#Vizier.VIZIER_SERVER = 'vizier.ast.cam.ac.uk'

def query(raImage, decImage,boxsize):
    #boxsize = 33.4  # Set the box size to search for catalog stars
    maxmag = 18     # Magnitude cut-offs of sources to be cross-matched against
    catNum = 'II/246'  # changing to 2mass
    print('\nQuerying Vizier %s around RA %.4f, Dec %.4f with a radius of %.4f arcmin' % (
    catNum, raImage, decImage, boxsize))
    try:
        #You can set the filters for the individual columns (magnitude range, number of detections) inside the Vizier query
        v = Vizier(columns=['*'], column_filters={"gmag":"<%.2f"%maxmag, "Nd":">6", "e_gmag":"<1.086/3"}, row_limit=-1)
        Q = v.query_region(SkyCoord(ra = raImage, dec = decImage, unit = (u.deg, u.deg)), radius = str(boxsize)+'m', catalog=catNum, cache=False)
        #query vizier around (ra, dec) with a radius of boxsize
        #print(Q[0])
        print(Q[0].colnames)
    except:
        print('I cannnot reach the Vizier database. Is the internet working?')
    return Q

#Q = query(catNum, raImage, decImage)
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
        print('Executing command: %s' % command)
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
        print('Executing command: %s' % command)
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
        print("Executing command: %s" % command)
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
def zeropt(good_cat_stars,cleanPSFSources,idx_psfmass,idx_psfimage,imageName,filter):
    psfoffsets = ma.array(good_cat_stars['%smag' % filter][idx_psfmass] - cleanPSFSources['MAG_POINTSOURCE'][idx_psfimage]) # *TO WORK ON* currently need to change the filter here manually
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

    cleanPSFSources.write('%s.mag.ecsv' % imageName,overwrite=True)
    print('%s.mag.ecsv written, CSV w/ corrected mags' % imageName)

#%%
defaults = dict(crop=300)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='runs sextractor and psfex on swarped img to get psf fit photometry, then '
                                                 'writes zero pts to img header and outputs csv w/ corrected mags')
    parser.add_argument('-dir', type=str, help='[str], directory of stacked image')
    parser.add_argument('-name', type=str, help='[str], image name')
    parser.add_argument('-filter', type=str, help='[str], filter, ex. "J"')
    parser.add_argument('-crop', type=str, help='[int], # of pixels from edge of image to crop, default 300',default=defaults["crop"])
    args = parser.parse_args()

    data, header, w, raImage, decImage, boxsize = img(args.dir, args.name, args.crop)
    Q = query(raImage, decImage, boxsize)
    catalogName = sex1(args.name)
    psfex(catalogName)
    psfcatalogName = sex2(args.name)
    good_cat_stars, cleanPSFSources, idx_psfmass, idx_psfimage = tables(Q, psfcatalogName, args.crop)
    zeropt(good_cat_stars, cleanPSFSources, idx_psfmass, idx_psfimage, args.name, args.filter)

#%%

"""dir = '/Users/orion/Desktop/PRIME/GRB/preproc_swarp/'
imageName = 'coadd_cds.fits'
configFile = '/Users/orion/Desktop/PRIME/sex.config'
catalogName = imageName+'.cat'
paramName = '/Users/orion/Desktop/PRIME/tempsource.param'
psfName = imageName + '.psf'
psfcatalogName = imageName.replace('.fits','.psf.cat')
psfparamName = '/Users/orion/Desktop/PRIME/photomPSF.param'"""
