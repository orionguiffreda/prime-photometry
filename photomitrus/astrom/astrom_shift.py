from astroquery.vizier import Vizier
from astropy.io import fits
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
import os
import subprocess
import argparse
from astropy.table import Table
import numpy as np
from astropy.table import Column

from settings import gen_config_file_name

#%% image info
def imaging(directory,imageName):
    os.chdir(directory)
    f = fits.open(imageName)
    data = f[0].data  #This is the image array
    header = f[0].header

    #strong the image WCS into an object
    w = WCS(header)

    #Get the RA and Dec of the center of the image
    [raImage, decImage] = w.all_pix2world(data.shape[0] / 2, data.shape[1] / 2, 1)

    return data,header,w,raImage,decImage

#data,header,w,raImage,decImage = imaging('/mnt/d/PRIME_photometry_test_files/astrom_shift/','01504054C1.sky.flat.fits')

#%% catalog query
def cat_query(raImage, decImage, filter):
    boxsize = 9
    catNum = 'II/246'  # changing to 2mass
    print('\nQuerying Vizier %s around RA %.4f, Dec %.4f, w/ box size %.2f' % (
        catNum, raImage, decImage, boxsize))
    try:
        # You can set the filters for the individual columns (magnitude range, number of detections) inside the Vizier query
        v = Vizier(columns=['*'], column_filters={"%smag" % filter: ">12", "Nd": ">6"}, row_limit=-1)
        Q = v.query_region(SkyCoord(ra=raImage, dec=decImage, unit=(u.deg, u.deg)), width=str(boxsize) + 'm', catalog=catNum, cache=False)
        # query vizier around (ra, dec) with a radius of boxsize
        # print(Q[0])
        print('Queried source total = ', len(Q[0]))
    except:
        print('I cannnot reach the Vizier database. Is the internet working?')
    return Q

#%% sextraction
def sex1(imageName):
    print('Running sextractor for psf...')
    configFile = gen_config_file_name('sex.config')
    paramName = gen_config_file_name('tempsource.param')
    catalogName = imageName + '.shift.cat'
    try:
        command = 'sex %s -c %s -CATALOG_NAME %s -PARAMETERS_NAME %s' % (imageName, configFile, catalogName, paramName)
        #print('Executing command: %s' % command)
        rval = subprocess.run(command.split(), check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError as err:
        print('Could not run sextractor with exit error %s'%err)
    return catalogName

def psfex(catalogName):
    print('Getting psf...')
    psfConfigFile = gen_config_file_name('default.psfex')
    try:
        command = 'psfex %s -c %s' % (catalogName,psfConfigFile)
        #print('Executing command: %s' % command)
        rval = subprocess.run(command.split(), check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError as err:
        print('Could not run psfex with exit error %s'%err)

def sex2(imageName):
    print('Sextracting sources...')
    psfname = imageName + '.shift.psf'
    configFile = gen_config_file_name('sex_shift.config')
    paramName = gen_config_file_name('astromshift.param')
    catname = imageName + '.psf.shift.cat'
    try:
        command = 'sex %s -c %s -CATALOG_NAME %s -PSF_NAME %s -PARAMETERS_NAME %s' % (imageName, configFile, catname, psfname ,paramName)
        #print("Executing command: %s" % command)
        rval = subprocess.run(command.split(), check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError as err:
        print('Could not run sextractor with exit error %s'%err)
    return catname

#%% creating & prepping tables for dist calc
def make_tables(directory, data, w, catname, Q, filter):
    print('Creating sorted catalog & prime tables...')
    sexcat = Table.read(directory+catname, hdu=2)

    crop = 1600
    max_x = data.shape[0]
    max_y = data.shape[1]
    mass_imCoords = w.all_world2pix(Q[0]['RAJ2000'], Q[0]['DEJ2000'], 1)
    inner_catsources = Q[0][np.where((mass_imCoords[0] > crop) & (mass_imCoords[0] < (max_x-crop)) & (mass_imCoords[1] > crop)
                                   & (mass_imCoords[1] < (max_y-crop)))]

    inner_primesources = sexcat[(sexcat['X_IMAGE'] < (max_x - crop)) & (sexcat['X_IMAGE'] > crop)
                                & (sexcat['Y_IMAGE'] < (max_y) - crop) & (
                                            sexcat['Y_IMAGE'] > crop) & (sexcat['FLUX_RADIUS'] > 1.5)]

    inner_primesources.sort('MAG_AUTO')

    inner_catsources.sort('%smag' % filter)
    cat_xy = w.all_world2pix(inner_catsources['RAJ2000'],inner_catsources['DEJ2000'],1)

    xs = Column(cat_xy[0], name='X_IMAGE',unit='pix')
    ys = Column(cat_xy[1], name='Y_IMAGE', unit='pix')
    inner_catsources.add_column(xs)
    inner_catsources.add_column(ys)
    return inner_primesources,inner_catsources


def prep_tables(inner_primesources,inner_catsources,num):
    first_primes_x = np.array(inner_primesources['X_IMAGE'][:num])
    first_primes_y = np.array(inner_primesources['Y_IMAGE'][:num])
    first_cats_x = np.array(inner_catsources['X_IMAGE'][:num])
    first_cats_y = np.array(inner_catsources['Y_IMAGE'][:num])

    first_primecoords = np.column_stack([first_primes_x,first_primes_y])
    first_catcoords = np.column_stack([first_cats_x,first_cats_y])
    return first_primecoords,first_catcoords

#%% distance calc
def calculate_distance(p1, p2):
    return np.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2)


def find_agreeing_distances(table1, table2, acc_range, length):
    """Find distances between the first 5 rows of two tables and check for agreeing distances within 3 pixels."""
    print('Computing distances...')
    # Extract the first 3 rows from each table
    points1 = table1
    points2 = table2

    distances = {}
    row_pairs = []

    # Compute distances between each pair of rows from table1 with each pair of rows from table2
    for i in range(len(points1)):
        for j in range(len(points2)):
            distance = calculate_distance(points1[i], points2[j])
            row_pairs.append(((i, j), distance))
            if distance not in distances:
                distances[distance] = []
            distances[distance].append((i, j))

    # Find distances that agree within 5 pixels
    agreeing_pairs = []
    sorted_distances = sorted(distances.keys())

    for i in range(len(sorted_distances)):
        for j in range(i + 1, len(sorted_distances)):
            if abs(sorted_distances[i] - sorted_distances[j]) <= acc_range:
                if sorted_distances[i] <= length:
                    agreeing_pairs.extend(distances[sorted_distances[i]])
                    agreeing_pairs.extend(distances[sorted_distances[j]])

    # Remove duplicates
    agreeing_pairs = list(set(agreeing_pairs))

    # Extract distances
    agreeing_distances = [pair for pair in row_pairs if pair[0] in agreeing_pairs]

    agreeing_pairs.sort()
    print(f"\nEuclidean distances that agree within {acc_range} pixels & < {length} pixels in length:")
    for pair, dist in zip(agreeing_pairs, agreeing_distances):
        print(f"Row from prime: {pair[0]}, Row from 2mass: {pair[1]}, Distance (pix): {dist[1]}")

    return agreeing_pairs, agreeing_distances


#%% final shift calc and header update
def xyshift(pairs, inner_primesources, inner_catsources, directory, header, data, imageName, stdev):
    indices_prime = [pair[0] for pair in pairs]
    indices_cat = [pair[1] for pair in pairs]

    filtered_prime = inner_primesources[indices_prime]
    filtered_cat = inner_catsources[indices_cat]

    #xy_shifts = []
    x_shiftarr = []
    y_shiftarr = []
    for i,j in zip(filtered_prime,filtered_cat):
        x_shift = i['X_IMAGE'] - j['X_IMAGE']
        y_shift = i['Y_IMAGE'] - j['Y_IMAGE']
        # xy_shift = tuple((x_shift,y_shift))
        x_shiftarr.append(x_shift)
        y_shiftarr.append(y_shift)
        # xy_shifts.append(xy_shift)

    # pruning outliers
    print('\npruning outliers >= %.2f stdev from median for x & y offsets...' % stdev)
    x_stdev = np.std(x_shiftarr)
    x_med = np.median(x_shiftarr)
    y_stdev = np.std(y_shiftarr)
    y_med = np.median(y_shiftarr)

    #x_shiftarr = [x for x in x_shiftarr if (x_med + x_stdev) >= x >= (x_med - x_stdev)]
    for i in range(len(x_shiftarr)):
        if x_shiftarr[i] <= (stdev*x_med - x_stdev) or x_shiftarr[i] >= (stdev*x_med + x_stdev):
            x_shiftarr[i] = 0
    #y_shiftarr = [y for y in y_shiftarr if (y_med + x_stdev) >= y >= (y_med - y_stdev)]
    for i in range(len(y_shiftarr)):
        if y_shiftarr[i] <= (stdev*y_med - y_stdev) or y_shiftarr[i] >= (stdev*y_med + y_stdev):
            y_shiftarr[i] = 0

    x_shiftarr_prune = []
    y_shiftarr_prune = []

    for i in range(len(x_shiftarr)):
        if x_shiftarr[i] != 0 and y_shiftarr[i] != 0:
            x_shiftarr_prune.append(x_shiftarr[i])
            y_shiftarr_prune.append(y_shiftarr[i])

    xy_shifts = np.column_stack([x_shiftarr_prune,y_shiftarr_prune])

    print('Sources considered: ')
    for pt in xy_shifts:
        print('X offset = %.3f, Y offset = %.3f' % (pt[0], pt[1]))

    xfinal_shift = np.median(x_shiftarr_prune)
    yfinal_shift = np.median(y_shiftarr_prune)
    print('\nfinal median shifts: x = %.3f, y = %.3f' % (xfinal_shift, yfinal_shift))

    crpix1 = header['CRPIX1']
    crpix2 = header['CRPIX2']

    header['CRPIX1'] = crpix1 + xfinal_shift
    header['CRPIX2'] = crpix2 + yfinal_shift

    imageshiftname = os.path.splitext(imageName)[0]
    imageshiftname = imageshiftname + '.shift.fits'

    print('Writing new FITS file w/ updated CRPIX: %s' % imageshiftname)
    newpath = os.path.join(directory, imageshiftname)
    oldpath = os.path.join(directory, imageName)

    fits.writeto(newpath, data, header, overwrite=True)

    return xfinal_shift, yfinal_shift
#%%


def change_all_files(xfinal_shift, yfinal_shift, directory):
    all_fits = [f for f in sorted(os.listdir(directory)) if f.endswith('.flat.fits')]
    all_fits = all_fits[1:]     # all files but first one (first one is completed already)

    print('Rewriting all FITS images w/ new CRPIX vals...')
    for f in all_fits:
        imgpath = os.path.join(directory,f)
        img = fits.open(imgpath)
        data = img[0].data
        header = img[0].header

        crpix1 = header['CRPIX1']
        crpix2 = header['CRPIX2']
        header['CRPIX1'] = crpix1 + xfinal_shift
        header['CRPIX2'] = crpix2 + yfinal_shift

        imageshiftname = os.path.splitext(f)[0]
        imageshiftname = imageshiftname + '.shift.fits'
        newpath = os.path.join(directory, imageshiftname)

        fits.writeto(newpath, data, header, overwrite=True)

    print('Moving old files to %sold/ directory and renaming shifted images...' % directory)
    os.mkdir(directory+'old/')
    all_fits_again = [f for f in sorted(os.listdir(directory)) if f.endswith('.flat.fits')]
    all_fits_shift = [f for f in sorted(os.listdir(directory)) if f.endswith('.shift.fits')]
    for f in all_fits_again:
        currentpath = os.path.join(directory,f)
        oldpath = os.path.join(directory+'old/',f)
        os.rename(currentpath,oldpath)

    for f in all_fits_shift:
        shiftpath = os.path.join(directory, f)
        imagerename = os.path.splitext(f)[0]
        fits_end = os.path.splitext(f)[1]
        imagerename = imagerename[:-6]
        imagerename = imagerename + fits_end
        newpath = os.path.join(directory,imagerename)
        os.rename(shiftpath,newpath)

#%% intermediate file removal
def removal(directory):
    fnames = ['.shift.cat','.psf']
    for f in os.listdir(directory):
        for name in fnames:
            if f.endswith(name):
                path = os.path.join(directory+f)
                try:
                    os.remove(path)
                    #print(f"Removed file: {path}")
                except Exception as e:
                    print(f"Error removing file: {path} - {e}")
#%%

defaults = dict(range=3,length=100,num=10,stdev=1)

#%%
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Corrects for translation in initial astrometry '
                                                 '(so astrom.net doesnt need to be used)')
    parser.add_argument('-pipeline', action='store_true', help='For use in the pipeline, runs on all fits imgs'
                                                               'in the given directory')
    parser.add_argument('-remove', action='store_true', help='removes intermediate catalogs')
    parser.add_argument('-dir', type=str, help='[str] path where input file is stored (should run on proc. image, '
                                               'so likely should be /C#_sub/)')
    parser.add_argument('-imagename', type=str, help='[str] input file name (should run on proc. image, '
                                                     'i.e. *.sky.flat.fits)')
    parser.add_argument('-filter', type=str, help='[str] filter used, ex. "J"')
    parser.add_argument('-range', type=float, help='[float] optional, range for eucl. dists to be considered '
                                                   'in agreement (pix), default = 3', default=defaults['range'])
    parser.add_argument('-length', type=float, help='[float] optional, value over which dists will not be '
                                                   'considered (pix), default = 100', default=defaults['length'])
    parser.add_argument('-num', type=int, help='[int] optional, # of sources, sorted by mag, to consider in '
                                                   'dist. calculation (dont make too large!), default = 10', default=defaults['num'])
    parser.add_argument('-stdev', type=float, help='[float] # of stdevs away from median to prune eucl. dists '
                                                   'default = 1', default=defaults['stdev'])
    args = parser.parse_args()

    data, header, w, raImage, decImage = imaging(args.dir, args.imagename)
    Q = cat_query(raImage, decImage, args.filter)
    catalogName = sex1(args.imagename)
    psfex(catalogName)
    catname = sex2(args.imagename)
    inner_primesources, inner_catsources = make_tables(args.dir, data, w, catname, Q, args.filter)
    first_primecoords, first_catcoords = prep_tables(inner_primesources, inner_catsources, args.num)
    agreeing_pairs, dists = find_agreeing_distances(first_primecoords, first_catcoords, args.range, args.length)
    xfinal_shift, yfinal_shift = xyshift(agreeing_pairs, inner_primesources, inner_catsources, args.dir, header, data,
                                         args.imagename, args.stdev)
    if args.remove:
        removal(args.dir)
    if args.pipeline:
        change_all_files(xfinal_shift, yfinal_shift, args.dir)

