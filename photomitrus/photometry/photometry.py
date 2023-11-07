"""
Calibrates photometry for stacked image
"""

import numpy as np
import numpy.ma as ma
import argparse
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.stats import sigma_clipped_stats, sigma_clip
import subprocess
from astropy.io import fits
import os

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
dir = '/Users/orion/Desktop/PRIME/GRB/preproc_swarp/'
imageName = 'coadd_cds.fits'
# Set the box size to search for catalog stars
boxsize = 33.4  # arcminutes

# Magnitude cut-offs of sources to be cross-matched against
maxmag = 18
def img(dir,imageName):
    os.chdir(dir)
    f = fits.open(imageName)
    data = f[0].data  #This is the image array
    header = f[0].header

    #strong the image WCS into an object
    w = WCS(header)

    #Get the RA and Dec of the center of the image
    [raImage, decImage] = w.all_pix2world(data.shape[0]/2, data.shape[1]/2, 1)

    return data,header,w,[raImage, decImage]

data,header,w,[raImage, decImage] = img(dir,imageName)
#%%
# Use astroquery to get catalog search
from astroquery.vizier import Vizier
#Vizier.VIZIER_SERVER = 'vizier.ast.cam.ac.uk'

catNum = 'II/246'#changing to 2mass
def query(catNum, raImage, decImage):
    boxsize = 33.4  # arcminutes
    maxmag = 18
    print('\nQuerying Vizier %s around RA %.4f, Dec %.4f with a radius of %.4f arcmin' % (
    catNum, raImage, decImage, boxsize))
    try:
        #You can set the filters for the individual columns (magnitude range, number of detections) inside the Vizier query
        v = Vizier(columns=['*'], column_filters={"gmag":"<%.2f"%maxmag, "Nd":">6", "e_gmag":"<1.086/3"}, row_limit=-1)
        Q = v.query_region(SkyCoord(ra = raImage, dec = decImage, unit = (u.deg, u.deg)), radius = str(boxsize)+'m', catalog=catNum, cache=False)
        #query vizier around (ra, dec) with a radius of boxsize
        print(Q[0])
        print(Q[0].colnames)
    except:
        print('I cannnot reach the Vizier database. Is the internet working?')
    return Q

Q = query(catNum, raImage, decImage)
#%%
#run sextractor on swarped img to find sources
configFile = '/Users/orion/Desktop/PRIME/sex.config'
catalogName = imageName+'.cat'
paramName = '/Users/orion/Desktop/PRIME/tempsource.param'
def sex1(imageName, configFile, paramName):
    catalogName = imageName + '.cat'
    try:
        command = 'sex %s -c %s -CATALOG_NAME %s -PARAMETERS_NAME %s' % (imageName, configFile, catalogName, paramName)
        print('Executing command: %s' % command)
        rval = subprocess.run(command.split(), check=True)
    except subprocess.CalledProcessError as err:
        print('Could not run sextractor with exit error %s'%err)

#%%
#run psfex on sextractor LDAC from previous step
psfConfigFile = '/Users/orion/Desktop/PRIME/default.psfex'
os.chdir('/Users/orion/Desktop/PRIME/GRB/preproc_swarp/')
def psfex(catalogName,psfConfigFile):
    try:
        command = 'psfex %s -c %s' % (catalogName,psfConfigFile)
        print('Executing command: %s' % command)
        rval = subprocess.run(command.split(), check=True)
    except subprocess.CalledProcessError as err:
        print('Could not run psfex with exit error %s'%err)

#%%
#feed generated psf model back into sextractor w/ diff param (or you could use that param from the start but its slower)
psfName = imageName + '.psf'
psfcatalogName = psfName+'.cat'
psfparamName = '/Users/orion/Desktop/PRIME/photomPSF.param' #This is a new set of parameters to be obtained from SExtractor, including PSF-fit magnitudes
os.chdir('/Users/orion/Desktop/PRIME/GRB/preproc_swarp/')
def sex2(imageName, configFile, psfparamName):
    psfName = imageName + '.psf'
    psfcatalogName = psfName + '.cat'
    try:
        #We are supplying SExtactor with the PSF model with the PSF_NAME option
        command = 'sex %s -c %s -CATALOG_NAME %s -PSF_NAME %s -PARAMETERS_NAME %s' % (imageName, configFile, psfcatalogName, psfName, psfparamName)
        print("Executing command: %s" % command)
        rval = subprocess.run(command.split(), check=True)
    except subprocess.CalledProcessError as err:
        print('Could not run sextractor with exit error %s'%err)
    return psfcatalogName

#%%
#Rread in tables, clean edges (potentially unnecessary), crossmatch
def tables(Q,psfcatalogName):
    mass_imCoords = w.all_world2pix(Q[0]['RAJ2000'], Q[0]['DEJ2000'], 1)
    good_cat_stars = Q[0][np.where((mass_imCoords[0] > 215) & (mass_imCoords[0] < 4400) & (mass_imCoords[1] > 215) & (mass_imCoords[1] < 4400))]

    psfsourceTable = get_table_from_ldac(psfcatalogName)
    cleanPSFSources = psfsourceTable[(psfsourceTable['FLAGS']==0) & (psfsourceTable['FLAGS_MODEL']==0) & (psfsourceTable['XMODEL_IMAGE']<4400) & (psfsourceTable['XMODEL_IMAGE']>215) &(psfsourceTable['YMODEL_IMAGE']<4400) &(psfsourceTable['YMODEL_IMAGE']>215)]

    psfsourceCatCoords = SkyCoord(ra=cleanPSFSources['ALPHA_J2000'], dec=cleanPSFSources['DELTA_J2000'], frame='icrs', unit='degree')

    massCatCoords = SkyCoord(ra=good_cat_stars['RAJ2000'], dec=good_cat_stars['DEJ2000'], frame='icrs', unit='degree')
    #Now cross match sources
    #Set the cross-match distance threshold to 0.6 arcsec, or just about one pixel
    photoDistThresh = 0.6
    idx_psfimage, idx_psfmass, d2d, d3d = massCatCoords.search_around_sky(psfsourceCatCoords, photoDistThresh*u.arcsec)
    #idx_psfimage are indexes into psfsourceCatCoords for the matched sources, while idx_psfmass are indexes into massCatCoords for the matched sources

    print('Found %d good cross-matches'%len(idx_psfmass))
    return good_cat_stars,cleanPSFSources,idx_psfmass,idx_psfimage

good_cat_stars,cleanPSFSources,idx_psfmass,idx_psfimage = tables(Q,psfcatalogName)
#%%
#derive zero pt / put in swarped header
def zeropt(good_cat_stars,cleanPSFSources,idx_psfmass,idx_psfimage):
    psfoffsets = ma.array(good_cat_stars['Jmag'][idx_psfmass] - cleanPSFSources['MAG_POINTSOURCE'][idx_psfimage])
    #Compute sigma clipped statistics
    zero_psfmean, zero_psfmed, zero_psfstd = sigma_clipped_stats(psfoffsets)
    print('PSF Mean ZP: %.2f\nPSF Median ZP: %.2f\nPSF STD ZP: %.2f'%(zero_psfmean, zero_psfmed, zero_psfstd))

    header.set('ZP_ERR',zero_psfstd,'PSF Zero Point STD',after='SATURATE')
    header.set('ZP_MEAN',zero_psfmean,'PSF Zero Point Mean',after='SATURATE')
    header.set('ZP_MEDIAN',zero_psfmed,'PSF Zero Point Median',after='SATURATE')

    print(header['ZP_ERR'],header['ZP_MEAN'],header['ZP_MEDIAN'])

#%%
#%%
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='runs sextractor and psfex on swarped img to get psf fit photometry')
    parser.add_argument('args', nargs='*', type=str, metavar='a', help='Put in order: dir,imagename,configfile,psfconfigfile'
                                                                       'paramname,psfparamname')
    args = parser.parse_args()

    data, header, w, [raImage, decImage] = img(dir, imageName)
    Q = query(catNum, raImage, decImage)
    sex1(imageName, configFile, paramName)
    psfex(catalogName, psfConfigFile)
    psfcatalogName = sex2(imageName, configFile, psfparamName)
    good_cat_stars, cleanPSFSources, idx_psfmass, idx_psfimage = tables(Q, psfcatalogName)
    zeropt(good_cat_stars, cleanPSFSources, idx_psfmass, idx_psfimage)