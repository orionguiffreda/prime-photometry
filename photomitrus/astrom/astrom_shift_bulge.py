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
from itertools import combinations
from collections import defaultdict
import threading
import math
import sys

from photomitrus.settings import (gen_config_file_name, CHIP_ZPS)


def get_zp(chip, band, vary=False):
    zps = [v for k, v in CHIP_ZPS.items() if k == band]
    chosen_zp = zps[0][chip-1]

    if vary:
        num = 5  # num of zp iterations
        step = 1  # step size
        start = chosen_zp - (num // 2) * step
        return np.arange(start, start + 5 * step, step)

    return chosen_zp


def timed_input(prompt, timeout=60, default='Y'):
    answer = [default]

    def ask():
        ans = input(prompt)
        if ans:
            answer[0] = ans

    thread = threading.Thread(target=ask)
    thread.daemon = True
    thread.start()
    thread.join(timeout)
    return answer[0]


#%% image info


def imaging(directory, imageName, vary=False):
    os.chdir(directory)
    f = fits.open(os.path.join(directory, imageName))
    data = f[0].data  # This is the image array
    header = f[0].header

    # strong the image WCS into an object
    w = WCS(header)

    # Get the RA and Dec of the center of the image
    [raImage, decImage] = w.all_pix2world(data.shape[0] / 2, data.shape[1] / 2, 1)

    # Get zero point for image
    chip = header['CHIP']
    if header['FILTER2'] == 'Open':
        band = 'Z'
    else:
        band = header['FILTER2']

    zp = get_zp(chip, band, vary)

    # cat name
    pre = os.path.splitext(imageName)[0]
    catname = pre+'.cat'
    return data, header, w, raImage, decImage, zp, catname


#%% catalog query


def cat_query(raImage, decImage, band, boxsize, maglow=12.5, maghigh=14.5):
    columns = ['RAJ2000', 'DEJ2000', 'RAICRS', 'DEICRS', "%sPSF" % band.lower(), "%smag" % band, "%smag3" % band,
               'srcid']
    coords = SkyCoord(ra=raImage*u.degree, dec=decImage*u.degree, frame='fk5')
    galcoords = coords.galactic
    print('Using bulge field query, converting coords to galactic: l = %s, b = %s' % (galcoords.l.deg, galcoords.b.deg))
    catNum = 'II/348/vvv2'  # changing to VVV
    print('\nQuerying Vizier %s around l %.4f, b %.4f, w/ box size %.2f, mag lim of %s - %s' % (catNum,
        galcoords.l.deg, galcoords.b.deg, boxsize, maglow, maghigh))
    try:
        v = Vizier(columns=columns, column_filters={"%smag3" % band: "%s .. %s" % (maglow, maghigh),
                                                    "%sperrbits" % band: "<=16", "Nd": ">6"}, row_limit=-1)
        Q = v.query_region(SkyCoord(galcoords, unit=(u.deg, u.deg)), width=str(boxsize) + 'm',
                           catalog=catNum, cache=False, frame='galactic')
        print('Queried source total = ', len(Q[0]))
    except:
        print('Issue with bulge-specific query.')
    return Q

# Q = cat_query(266, -29, 'J', 15, bulge=True)
# print(Q[0])
# data, header, w, raImage, decImage = (
#     imaging('/mnt/photometry/supermaster_test/GB63_dither_20250406/J/C1_sub/','02284022C1.sky.flat.fits'))
# Q = cat_query(raImage, decImage, 'J', 15, bulge=True)
# f = fits.open('/mnt/photometry/supermaster_test/GB63_dither_20250406/J/C1_sub/02284022C1.sky.flat.fits')
# hdr = f[0].header
#
# # strong the image WCS into an object
# w = WCS(hdr)
# Coords = w.all_world2pix(Q[0]['RAJ2000'], Q[0]['DEJ2000'], 1)
# Path = '/mnt/photometry/supermaster_test/GB63_dither_20250406/J/C1_sub/catcoords.reg'
# #coords_prime = coords_prime.T
# newtext = open(Path, 'w+')
# for i,j in zip(Coords[0],Coords[1]):
#     newtext.write('\npoint(%f,%f) # point=circle 5' % (i,j))
#%% sextraction / psfex


def sex1(imageName):
    print('Running sextractor for psf...')
    configFile = gen_config_file_name('bulge_new.config')
    paramName = gen_config_file_name('tempsource.param')
    catname = imageName + '.cat'
    try:
        command = 'sex %s -c %s -CATALOG_NAME %s -PARAMETERS_NAME %s -CHECKIMAGE_TYPE NONE' % (imageName, configFile, catname, paramName)
        #print('Executing command: %s' % command)
        rval = subprocess.run(command.split(), check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError as err:
        print('Could not run sextractor with exit error %s'%err)
    return catname

# def psfex(catalogName):
#     print('Getting psf...')
#     psfConfigFile = gen_config_file_name('default.psfex')
#     try:
#         command = 'psfex %s -c %s' % (catalogName,psfConfigFile)
#         #print('Executing command: %s' % command)
#         rval = subprocess.run(command.split(), check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
#     except subprocess.CalledProcessError as err:
#         print('Could not run psfex with exit error %s'%err)
#
#
# def sex2(imageName):
#     print('Sextracting sources...')
#     psfname = imageName + '.shift.psf'
#     configFile = gen_config_file_name('sex_shift.config')
#     paramName = gen_config_file_name('astromshift.param')
#     catname = imageName + '.psf.shift.cat'
#     try:
#         command = 'sex %s -c %s -CATALOG_NAME %s -PSF_NAME %s -PARAMETERS_NAME %s' % (
#             imageName, configFile, catname, psfname, paramName)
#         # print("Executing command: %s" % command)
#         rval = subprocess.run(command.split(), check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
#     except subprocess.CalledProcessError as err:
#         print('Could not run sextractor with exit error %s'%err)
#     return catname

#%% creating & prepping tables for dist calc


def make_tables(directory, data, w, catname, Q, band, crop, zp, maglow=12.5, maghigh=14.5):
    print('Creating sorted catalog & prime tables...')
    sexcat = Table.read(os.path.join(directory, catname), hdu=2)
    colnames = Q[0].colnames

    max_x = data.shape[0]
    max_y = data.shape[1]
    mass_imCoords = w.all_world2pix(Q[0][colnames[0]], Q[0][colnames[1]], 1)
    inner_catsources = Q[0][np.where((mass_imCoords[0] > crop) & (mass_imCoords[0] < (max_x-crop)) & (mass_imCoords[1] > crop)
                                   & (mass_imCoords[1] < (max_y-crop)))]
    # inner_catsources = inner_catsources[(inner_catsources[colnames[2]] * 0.339 > 1)]

    inner_primesources = sexcat[(sexcat['FLAGS'] == 0) &
                                (sexcat['X_IMAGE'] < (max_x - crop)) & (sexcat['X_IMAGE'] > crop)
                                & (sexcat['Y_IMAGE'] < (max_y) - crop) & (
                                            sexcat['Y_IMAGE'] > crop) & (sexcat['FLUX_RADIUS'] * 0.498 > 1)]

    # zp correction
    print('applying zp correction to PRIME mags: %s' % zp)
    inner_primesources['MAG_AUTO'] = inner_primesources['MAG_AUTO'] + zp

    inner_primesources.write('prime_all.ecsv', overwrite=True)
    inner_primesources = inner_primesources[(inner_primesources['MAG_AUTO'] >= maglow) &
                                            (inner_primesources['MAG_AUTO'] <= maghigh)]

    # inner_primesources.sort('MAG_AUTO')
    # inner_catsources.sort(colnames[2])
    cat_xy = w.all_world2pix(inner_catsources[colnames[0]], inner_catsources[colnames[1]], 1)

    xs = Column(cat_xy[0], name='X_IMAGE', unit='pix')
    ys = Column(cat_xy[1], name='Y_IMAGE', unit='pix')
    inner_catsources.add_column(xs)
    inner_catsources.add_column(ys)

    inner_catsources.write('cat_all.ecsv', overwrite=True)

    Path = os.path.join(directory, 'catcoords_crop.reg')
    newtext = open(Path, 'w+')
    for i,j in zip(xs,ys):
        newtext.write('\npoint(%f,%f) # point=circle 5' % (i,j))

    Path = os.path.join(directory, 'primecoords_crop.reg')
    newtext = open(Path, 'w+')
    for i,j in zip(inner_primesources['X_IMAGE'],inner_primesources['Y_IMAGE']):
        newtext.write('\npoint(%f,%f) # point=circle 5' % (i,j))

    return inner_primesources, inner_catsources, colnames


def split_coordinates(n_segs, data, inner_primesources, inner_catsources, num, band, colnames, directory):
    # Calculate the number of segments along one dimension
    img_size = data.shape[0]
    segs = int(np.sqrt(n_segs))
    if segs * segs != n_segs:
        segs += 1

    # Calculate segment boundaries for x and y dimensions
    x_bounds = np.linspace(0, img_size, segs + 1, endpoint=True)
    y_bounds = np.linspace(0, img_size, segs + 1, endpoint=True)

    # Create a dictionary to hold the coordinates for each segment
    segments = {}

    for i in range(segs):
        for j in range(segs):
            # Define the boundaries for the current segment
            x_min = x_bounds[i]
            x_max = x_bounds[i + 1]
            y_min = y_bounds[j]
            y_max = y_bounds[j + 1]

            # Filter the coordinates within the current segment
            segment_primedata = inner_primesources[
                (inner_primesources['X_IMAGE'] >= x_min) & (inner_primesources['X_IMAGE'] < x_max) &
                (inner_primesources['Y_IMAGE'] >= y_min) & (inner_primesources['Y_IMAGE'] < y_max)
                ]

            segment_catdata = inner_catsources[
                (inner_catsources['X_IMAGE'] >= x_min) & (inner_catsources['X_IMAGE'] < x_max) &
                (inner_catsources['Y_IMAGE'] >= y_min) & (inner_catsources['Y_IMAGE'] < y_max)
                ]

            segment_primedata.sort('MAG_AUTO')
            segment_catdata.sort(colnames[2])

            # print(segment_primedata['MAG_AUTO', 'X_IMAGE', 'Y_IMAGE'])
            # print(segment_catdata)

            first_primes_x = np.array(segment_primedata['X_IMAGE'][:num])
            first_primes_y = np.array(segment_primedata['Y_IMAGE'][:num])
            first_cats_x = np.array(segment_catdata['X_IMAGE'][:num])
            first_cats_y = np.array(segment_catdata['Y_IMAGE'][:num])

            first_primecoords = np.column_stack([first_primes_x, first_primes_y])
            first_catcoords = np.column_stack([first_cats_x, first_cats_y])

            # Store the segment data in the dictionary
            segment_key = (i, j)
            segments[segment_key] = first_primecoords,first_catcoords,segment_primedata,segment_catdata

    # Path = os.path.join(directory, 'catcoords_crop.reg')
    # #coords_prime = coords_prime.T
    # newtext = open(Path, 'w+')
    # for i,j in zip(segment_catdata['X_IMAGE'],segment_catdata['Y_IMAGE']):
    #     newtext.write('\npoint(%f,%f) # point=circle 5' % (i,j))
    #
    # Path = os.path.join(directory, 'primecoords_crop.reg')
    # #coords_prime = coords_prime.T
    # newtext = open(Path, 'w+')
    # for i,j in zip(segment_primedata['X_IMAGE'],segment_primedata['Y_IMAGE']):
    #     newtext.write('\npoint(%f,%f) # point=circle 5' % (i,j))

    segnum = 0
    for key, values in segments.items():
        primecoords = values[0]
        pxs = []
        pys = []
        for coords in primecoords:
            pxs.append(coords[0])
            pys.append(coords[1])
        catcoords = values[1]
        cxs = []
        cys = []
        for coords in catcoords:
            cxs.append(coords[0])
            cys.append(coords[1])

        Path = os.path.join(directory, 'catcoords_sec_%s.reg' % segnum)
        #coords_prime = coords_prime.T
        newtext = open(Path, 'w+')
        for i,j in zip(cxs,cys):
            newtext.write('\npoint(%f,%f) # point=circle 5' % (i,j))

        Path = os.path.join(directory, 'primecoords_sec_%s.reg' % segnum)
        #coords_prime = coords_prime.T
        newtext = open(Path, 'w+')
        for i,j in zip(pxs,pys):
            newtext.write('\npoint(%f,%f) # point=circle 5' % (i,j))

        segnum += 1

    # segment_primedata.write('section_prime.ecsv', overwrite=True)

    return segments


def prep_tables(inner_primesources,inner_catsources,num):
    first_primes_x = np.array(inner_primesources['X_IMAGE'][:num])
    first_primes_y = np.array(inner_primesources['Y_IMAGE'][:num])
    first_cats_x = np.array(inner_catsources['X_IMAGE'][:num])
    first_cats_y = np.array(inner_catsources['Y_IMAGE'][:num])

    first_primecoords = np.column_stack([first_primes_x,first_primes_y])
    first_catcoords = np.column_stack([first_cats_x,first_cats_y])

    Path = '/mnt/photometry/supermaster_test/GB63_dither_20250406/J/C1_sub/catcoords_prune.reg'
    #coords_prime = coords_prime.T
    newtext = open(Path, 'w+')
    for i,j in zip(first_cats_x,first_cats_y):
        newtext.write('\npoint(%f,%f) # point=circle 5' % (i,j))

    Path = '/mnt/photometry/supermaster_test/GB63_dither_20250406/J/C1_sub/primecoords_prune.reg'
    #coords_prime = coords_prime.T
    newtext = open(Path, 'w+')
    for i,j in zip(first_primes_x,first_primes_y):
        newtext.write('\npoint(%f,%f) # point=circle 5' % (i,j))
    return first_primecoords,first_catcoords

#%% distance calc


def calculate_distance(p1, p2):
    return np.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2)


def find_agreeing_distances(table1, table2, acc_range, length):
    """Find distances between the first 5 rows of two tables and check for agreeing distances within 3 pixels."""
    print('Computing distances...')
    # Extract the first rows from each table
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

    # Find distances that agree w/in certain range, under certain length, etc.
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
        print(f"Row from prime: {pair[0]}, Row from catalog: {pair[1]}, Distance (pix): {dist[1]}")

    return agreeing_pairs, agreeing_distances


#%% final shift calc and header update


def xyshift(pairs, inner_primesources, inner_catsources, directory, header, data, imageName, stdev, segments=None):
    indices_prime = [pair[0] for pair in pairs]
    indices_cat = [pair[1] for pair in pairs]

    filtered_prime = inner_primesources[indices_prime]
    filtered_cat = inner_catsources[indices_cat]

    #xy_shifts = []
    x_shiftarr = []
    y_shiftarr = []
    for i, j in zip(filtered_prime,filtered_cat):
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
        if x_shiftarr[i] <= (x_med - stdev*x_stdev) or x_shiftarr[i] >= (x_med + stdev*x_stdev):
            x_shiftarr[i] = 0
    #y_shiftarr = [y for y in y_shiftarr if (y_med + x_stdev) >= y >= (y_med - y_stdev)]
    for i in range(len(y_shiftarr)):
        if y_shiftarr[i] <= (y_med - stdev*y_stdev) or y_shiftarr[i] >= (y_med + stdev*y_stdev):
            y_shiftarr[i] = 0

    x_shiftarr_prune = []
    y_shiftarr_prune = []

    for i in range(len(x_shiftarr)):
        if x_shiftarr[i] != 0 and y_shiftarr[i] != 0:
            x_shiftarr_prune.append(x_shiftarr[i])
            y_shiftarr_prune.append(y_shiftarr[i])

    if not x_shiftarr_prune:
        print('No pairs left after pruning!')
        xfinal_shift = 0
        yfinal_shift = 0
    else:

        xy_shifts = np.column_stack([x_shiftarr_prune,y_shiftarr_prune])

        print('Sources considered: ')
        for pt in xy_shifts:
            print('X offset = %.3f, Y offset = %.3f' % (pt[0], pt[1]))

        xfinal_shift = np.median(x_shiftarr_prune)
        yfinal_shift = np.median(y_shiftarr_prune)
        print('\nfinal median shifts: x = %.3f, y = %.3f' % (xfinal_shift, yfinal_shift))

        if not segments:
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
        else:
            pass

    return xfinal_shift, yfinal_shift
#%%


def segmentshiftcalcs(segments, acc_range, length, directory, header, data, imageName, stdev, segstd, segtrue):
    final_xshifts_arr_orig = []
    final_yshifts_arr_orig = []
    for seg in segments:
        table1 = segments[seg][0]
        table2 = segments[seg][1]
        print(f'Calculating agreeing pairs and dists for Section {seg}')
        pairs, distances = find_agreeing_distances(table1, table2, acc_range, length)

        print(f'Calculating median shifts for Section {seg}')
        xfinal_shift, yfinal_shift = xyshift(
            pairs, segments[seg][2], segments[seg][3], directory, header, data, imageName, stdev, segments=segtrue)
        final_xshifts_arr_orig.append(xfinal_shift)
        final_yshifts_arr_orig.append(yfinal_shift)

    final_xshifts_arr = final_xshifts_arr_orig.copy()
    final_yshifts_arr = final_yshifts_arr_orig.copy()

    x_finalshift_prune = final_xshifts_arr
    y_finalshift_prune = final_yshifts_arr

    return x_finalshift_prune, y_finalshift_prune


def segmentshift(segments, acc_range, length, directory, header, data, imageName, stdev, segstd, segtrue):
    final_xshifts_arr_orig = []
    final_yshifts_arr_orig = []
    for seg in segments:
        table1 = segments[seg][0]
        table2 = segments[seg][1]
        print(f'Calculating agreeing pairs and dists for Section {seg}')
        pairs, distances = find_agreeing_distances(table1, table2, acc_range, length)

        print(f'Calculating median shifts for Section {seg}')
        xfinal_shift, yfinal_shift = xyshift(
            pairs, segments[seg][2], segments[seg][3], directory, header, data, imageName, stdev, segments=segtrue)
        final_xshifts_arr_orig.append(xfinal_shift)
        final_yshifts_arr_orig.append(yfinal_shift)

    final_xshifts_arr = final_xshifts_arr_orig.copy()
    final_yshifts_arr = final_yshifts_arr_orig.copy()

    x_finalshift_prune = final_xshifts_arr
    y_finalshift_prune = final_yshifts_arr

    if not x_finalshift_prune:
        nonzero_idx = [idx for idx, val in enumerate(final_xshifts_arr_orig) if val != 0]
        if len(nonzero_idx) == 1:
            check = timed_input('Only 1 segment had a solution, do you want to go along with it? (Input Y or N): '
                                , timeout=60)
            if check.upper() == 'Y':
                idx = nonzero_idx[0]
                xfinal_shift = final_xshifts_arr_orig[idx]
                yfinal_shift = final_yshifts_arr_orig[idx]
            else:
                xfinal_shift = 0
                yfinal_shift = 0
        else:
            print('All segments disagree w/in %.2f stdevs! Is there an issue with the image?' % segstd)
            xfinal_shift = 0
            yfinal_shift = 0
    else:
        xy_shifts = np.column_stack([x_finalshift_prune,y_finalshift_prune])

        print('Final shifts considered: ')
        for pt in xy_shifts:
            print('X shift = %.3f, Y shift = %.3f' % (pt[0], pt[1]))

        xfinal_shift = np.median(x_finalshift_prune)
        yfinal_shift = np.median(y_finalshift_prune)
        print('\nfinal median shifts across all segments: x = %.3f, y = %.3f' % (xfinal_shift, yfinal_shift))

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


def change_all_files(xfinal_shift, yfinal_shift, directory, all_fits_arr=None):
    if xfinal_shift == 0:
        print('No agreement, thus cannot move forward with rewriting all files!')
    else:
        if all_fits_arr:
            all_fits = all_fits_arr
            # all_fits = all_fits[1:]
        else:
            all_fits = [f for f in sorted(os.listdir(directory)) if f.endswith('.flat.fits')]
            # all_fits = all_fits[1:]     # all files but first one (first one is completed already)

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

        print('Moving old files to %s/old/ directory and renaming shifted images...' % directory)
        old_storage_dir = os.path.join(directory, 'old')
        if not os.path.exists(old_storage_dir):
            os.mkdir(old_storage_dir)
        else:
            pass
        if all_fits_arr:
            all_fits_again = all_fits_arr
            all_fits_shift = []
            for f in all_fits_again:
                imageshiftname = os.path.splitext(f)[0]
                imageshiftname = imageshiftname + '.shift.fits'
                all_fits_shift.append(imageshiftname)
        else:
            all_fits_again = [f for f in sorted(os.listdir(directory)) if f.endswith('.flat.fits')]
            all_fits_shift = [f for f in sorted(os.listdir(directory)) if f.endswith('.shift.fits')]

        for f in all_fits_again:
            currentpath = os.path.join(directory, f)
            oldpath = os.path.join(old_storage_dir, f)
            os.rename(currentpath, oldpath)

        for f in all_fits_shift:
            shiftpath = os.path.join(directory, f)
            imagerename = os.path.splitext(f)[0]
            fits_end = os.path.splitext(f)[1]
            imagerename = imagerename[:-6]
            imagerename = imagerename + fits_end
            newpath = os.path.join(directory, imagerename)
            os.rename(shiftpath, newpath)

#%% intermediate file removal


def removal(directory):
    fnames = ['.shift.cat','.psf', '.reg']
    try:
        for f in os.listdir(directory):
            for name in fnames:
                if f.endswith(name):
                    path = os.path.join(directory, f)
                    try:
                        os.remove(path)
                        #print(f"Removed file: {path}")
                    except Exception as e:
                        print(f"Error removing file: {path} - {e}")
    except None as e:
        print('No files found to remove')
#%%


defaults = dict(range=3,length=100,num=15,stdev=1,segstd=2)

#%%

def boxchange(size):
    # size = size of box in arcmin
    boxsize_pix = size * 60 * (1/0.498)
    crop_pix = (4088 - boxsize_pix) / 2
    return size, crop_pix

def shiftiteration(
        directory, imagename, filter_used, crop, boxsize, n_segs, acc_range, length, num, stdev,
        segstd, maglow, maghigh, chosen_zp
):
    data, header, w, raImage, decImage, zp, catname = imaging(directory, imagename)
    Q = cat_query(raImage, decImage, filter_used, boxsize, maglow=maglow, maghigh=maghigh)
    inner_primesources, inner_catsources, colnames = make_tables(directory, data, w, catname, Q, filter_used, crop, chosen_zp,
                                                                 maglow=maglow, maghigh=maghigh)
    segments = split_coordinates(n_segs, data, inner_primesources, inner_catsources, num, filter_used, colnames, directory)
    xfinal_shifts, yfinal_shifts = segmentshiftcalcs(segments, acc_range, length, directory, header, data, imagename,
                                              stdev, segstd, segtrue=True)

    return xfinal_shifts, yfinal_shifts


def shift(
        directory, imagename, band, acc_range=3, length=100, num=15, stdev=1, segstd=2, test=False, vary=False
):
    if band == 'Y':
        filter_used = 'J'
    else:
        filter_used = band

    boxsize, crop = boxchange(3)
    n_segs = 4

    maglow = 12.5
    # maghigh = 13.5

    data, header, w, raImage, decImage, zp, catname = imaging(directory, imagename, vary)
    if vary:
        zp_vals = zp
    else:
        zp_vals = [zp]
    for zp in zp_vals:
        maghigh = 13.5
        print('\nVarying ZP = ',zp)
        Q = cat_query(raImage, decImage, filter_used, boxsize, maglow=maglow, maghigh=maghigh)
        # catname = sex1(imagename)
        # psfex(catalogName)
        # catname = sex2(imagename)
        inner_primesources, inner_catsources, colnames = make_tables(directory, data, w, catname, Q, filter_used, crop, zp,
                                                                     maglow=maglow, maghigh=maghigh)
        segments = split_coordinates(n_segs, data, inner_primesources, inner_catsources, num, filter_used, colnames, directory)
        prev_xfinal_shifts, prev_yfinal_shifts = segmentshiftcalcs(segments, acc_range, length, directory, header, data, imagename,
                                                  stdev, segstd, segtrue=True)

        # boxsize_arcmin = 3
        # boxsize, crop = boxchange(boxsize_arcmin)
        # maghigh = 13
        # print('\n Conducting 2nd iteration, mag range: %s - %s' % (maglow, maghigh))
        # xfinal_shifts, yfinal_shifts = shiftiteration(directory, imagename, filter_used, crop, boxsize, n_segs,
        #                                               acc_range, length, num, stdev, segstd, maglow, maghigh)

        first_iter_shifts = np.column_stack((prev_xfinal_shifts, prev_yfinal_shifts))
        # sec_iter_shifts = np.column_stack((xfinal_shifts, yfinal_shifts))
        master_shift_dict = {}
        master_shift_dict[0] = first_iter_shifts
        # master_shift_dict[1] = sec_iter_shifts

        # absdiff_xshifts = [abs(a - b) for a, b in zip(prev_xfinal_shifts, xfinal_shifts)]
        # absdiff_yshifts = [abs(a - b) for a, b in zip(prev_yfinal_shifts, yfinal_shifts)]
        #
        # like_shifts = [[], []]
        # threshold = 3
        #
        # for idx, val in enumerate(absdiff_xshifts):
        #     if val < threshold:
        #         like_shifts[0].append((idx, val))
        # for idx, val in enumerate(absdiff_yshifts):
        #     if val < threshold:
        #         like_shifts[1].append((idx, val))

        # matches = []
        # Check for matching indices and print message and values
        # for idx_x, val_x in like_shifts[0]:
        #     for idx_y, val_y in like_shifts[1]:
        #         if idx_x == idx_y:
        #             print(f"\nLike index found! Index: {idx_x}, dx: {val_x:.3f}, dy: {val_y:.3f}")
        #             matches.append((idx_x, (val_x, val_y)))

        master_idx_num = 1
        stop_flag = False
        chosen_final_shifts_x = []
        chosen_final_shifts_y = []

        while len(chosen_final_shifts_x) == 0:  # <- Continue looping until valid matches found
            all_matches = []

            if maghigh <= 12.75:
                print('High mag lim iteration has reached 12.5 mag, thus cannot continue with iterations!')
                stop_flag = True
                break

            maghigh -= 0.25
            print('\nNo matches found, iterating with mag range: %s - %s\n' % (maglow, maghigh))

            # shiftiteration function updates master_shift_dict[master_idx_num]
            xiter_shifts, yiter_shifts = shiftiteration(
                directory, imagename, filter_used, crop, boxsize, n_segs,
                acc_range, length, num, stdev, segstd, maglow=maglow,
                maghigh=maghigh, chosen_zp=zp
            )
            master_shift_dict[master_idx_num] = np.column_stack((xiter_shifts, yiter_shifts))

            threshold = 3

            # Compare newest set with previous sets
            for j in range(master_idx_num):
                set_i = master_shift_dict[j]
                set_j = master_shift_dict[master_idx_num]

                for idx_i, (x_i, y_i) in enumerate(set_i):
                    for idx_j, (x_j, y_j) in enumerate(set_j):
                        abs_dx = abs(x_i - x_j)
                        abs_dy = abs(y_i - y_j)

                        if abs_dx < threshold and abs_dy < threshold:
                            match_info = {
                                'pair': (j, master_idx_num),
                                'indices': (idx_i, idx_j),
                                'values_i': (x_i, y_i),
                                'values_j': (x_j, y_j),
                                'abs_dx': abs_dx,
                                'abs_dy': abs_dy
                            }
                            all_matches.append(match_info)
            for match in all_matches:
                print(f"Match between sets {match['pair']} at indices {match['indices']}:")
                print(
                    f"  Set {match['pair'][0]} [{match['indices'][0]}]: x = {match['values_i'][0]:.3f}, y = {match['values_i'][1]:.3f}")
                print(
                    f"  Set {match['pair'][1]} [{match['indices'][1]}]: x = {match['values_j'][0]:.3f}, y = {match['values_j'][1]:.3f}")
                print(f"  abs_dx = {match['abs_dx']:.3f}, abs_dy = {match['abs_dy']:.3f}\n")

            # Remove duplicate matches

            print('Removing duplicate matches by pair & index')
            match_counts = defaultdict(int)
            coord_key_map = {}  # Map key back to original match
            for m in all_matches:
                coord_i = (round(m['values_i'][0], 5), round(m['values_i'][1], 5))
                coord_j = (round(m['values_j'][0], 5), round(m['values_j'][1], 5))

                # Create order-independent key
                key = tuple(sorted([coord_i, coord_j]))

                match_counts[key] += 1
                coord_key_map[key] = m  # Save one copy of match for retrieval

            # Second pass: collect only truly unique matches
            unique_matches = [coord_key_map[key] for key, count in match_counts.items() if count == 1]
            # print(unique_matches)

            # extract values
            x_vals = []
            y_vals = []
            for match in unique_matches:
                x_vals.append(match['values_i'][0])
                y_vals.append(match['values_i'][1])
                x_vals.append(match['values_j'][0])
                y_vals.append(match['values_j'][1])

            # Median prune
            print('Pruning remaining matches outside of 3 pix from x and y medians...')
            if x_vals and y_vals:  # Check to avoid np.median([]) error

                # conglomerate all coord pairs
                all_coords = [(round(m['values_i'][0], 5), round(m['values_i'][1], 5)) for m in unique_matches] + \
                             [(round(m['values_j'][0], 5), round(m['values_j'][1], 5)) for m in unique_matches]

                # remove duplicates using a set
                unique_coords = list(set(all_coords))

                # split into separate x and y lists
                all_xs = [coord[0] for coord in unique_coords]
                all_ys = [coord[1] for coord in unique_coords]

                x_median = np.median(all_xs)
                print('X med: %.4f' % x_median)
                y_median = np.median(all_ys)
                print('Y med: %.4f' % y_median)

                chosen_final_shifts_x = []
                chosen_final_shifts_y = []

                # filter values w/in 3 pixels of median
                for x, y in zip(all_xs, all_ys):
                    if abs(x - x_median) <= 3 and abs(y - y_median) <= 3:
                        chosen_final_shifts_x.append(x)
                        chosen_final_shifts_y.append(y)

            if not chosen_final_shifts_x:
                print('No shifts satisfy criteria, continuing with loop...')
            else:
                print('Selected final shifts: ')
                for shiftx, shifty in zip(chosen_final_shifts_x, chosen_final_shifts_y):
                    print(' x: %.3f, y: %.3f' % (shiftx, shifty))
                ultimate_shift_x = np.median(chosen_final_shifts_x)
                ultimate_shift_y = np.median(chosen_final_shifts_y)
                print('\nUltimate final shift = X: %.3f, Y: %.3f' % (ultimate_shift_x, ultimate_shift_y))

            master_idx_num += 1

            # if len(all_matches) == 0:
            #     print('No real matches found, abs_dx & abs_dy = 0 for all, continuing...')
            # elif len(all_matches) > 1:
            #     print('Multiple shifts found. Checking for duplicates...')
            #     # Remove duplicates using a key of (pair, indices)
            #     seen_keys = set()
            #     unique_matches = []
            #     for m in all_matches:
            #         key = (m['pair'], m['indices'])
            #         if key not in seen_keys:
            #             seen_keys.add(key)
            #             unique_matches.append(m)
            #
            #
            #     # Extract all x and y values from the matches
            #     x_vals = []
            #     y_vals = []
            #     for match in unique_matches:
            #         x_vals.append(match['values_i'][0])
            #         y_vals.append(match['values_i'][1])
            #         x_vals.append(match['values_j'][0])
            #         y_vals.append(match['values_j'][1])
            #
            #     x_median = np.median(x_vals)
            #     y_median = np.median(y_vals)
            #
            #     # Filter matches based on proximity to median
            #     chosen_final_shifts_x = []
            #     chosen_final_shifts_y = []
            #
            #     for match in unique_matches:
            #         for x, y in [match['values_i'], match['values_j']]:
            #             if abs(x - x_median) <= 3 and abs(y - y_median) <= 3:
            #                 chosen_final_shifts_x.append(x)
            #                 chosen_final_shifts_y.append(y)
            #
            #     # Optional: Print final selections
            #     print("\nFiltered matches based on median proximity:")
            #     for x, y in zip(chosen_final_shifts_x, chosen_final_shifts_y):
            #         print(f"  x = {x:.3f}, y = {y:.3f}")
            #
            #     ultimate_shift_x = np.median(chosen_final_shifts_x)
            #     ultimate_shift_y = np.median(chosen_final_shifts_y)
            #     print('\nUltimate final shift = X: %.3f, Y: %.3f' % (ultimate_shift_x, ultimate_shift_y))
            # elif len(all_matches) == 1:
            #     print('Candidate matches found for below values!')
            #     # print(all_matches)
            #     chosen_pairs = [match for match in all_matches if match['abs_dx'] != 0]
            #     chosen_final_shifts_x = []
            #     chosen_final_shifts_y = []
            #     print(chosen_pairs)
            #     print(' x: %.3f, y: %.3f' % (chosen_pairs[0]['values_i'][0], chosen_pairs[0]['values_i'][1]))
            #     chosen_final_shifts_x.append(chosen_pairs[0]['values_i'][0])
            #     chosen_final_shifts_y.append(chosen_pairs[0]['values_i'][1])
            #
            #     print(' x: %.3f, y: %.3f' % (chosen_pairs[0]['values_j'][0], chosen_pairs[0]['values_j'][1]))
            #     chosen_final_shifts_x.append(chosen_pairs[0]['values_j'][0])
            #     chosen_final_shifts_y.append(chosen_pairs[0]['values_j'][1])
            #
            #     ultimate_shift_x = np.median(chosen_final_shifts_x)
            #     ultimate_shift_y = np.median(chosen_final_shifts_y)
            #     print('\nUltimate final shift = X: %.3f, Y: %.3f' % (ultimate_shift_x, ultimate_shift_y))
            #     break

            if stop_flag:
                ultimate_shift_x = 0
                ultimate_shift_y = 0

        # else:
        #     print('Like shifts include:')
        #     chosen_final_shifts_x = []
        #     chosen_final_shifts_y = []
        #     for chosen_set in master_shift_dict:
        #         chosen_shift = master_shift_dict[chosen_set][matches[0][0]]
        #         print(' x: %.3f, y: %.3f' % (chosen_shift[0], chosen_shift[1]))
        #         chosen_final_shifts_x.append(chosen_shift[0])
        #         chosen_final_shifts_y.append(chosen_shift[1])
        #
        #     ultimate_shift_x = np.median(chosen_final_shifts_x)
        #     ultimate_shift_y = np.median(chosen_final_shifts_y)
        #     print('\nUltimate final shift = X: %.3f, Y: %.3f' % (ultimate_shift_x, ultimate_shift_y))

    if not test:
        change_all_files(ultimate_shift_x, ultimate_shift_y, directory)



def main():
    parser = argparse.ArgumentParser(description='Corrects for translation in initial astrometry '
                                                 '(so astrom.net doesnt need to be used)')
    parser.add_argument('-test', action='store_true', help='optional flag to test for a successful solve, '
                                                           'doesnt write out any fits files.')
    parser.add_argument('-vary', action='store_true', help='optional flag to very zp')
    parser.add_argument('-dir', type=str, help='[str] path where input file is stored (should run on proc. image, '
                                               'so likely should be /C#_sub/)')
    parser.add_argument('-imagename', type=str, help='[str] input file name (should run on proc. image, '
                                                     'i.e. *.sky.flat.fits)')
    parser.add_argument('-band', type=str, help='[str] filter used, ex. "J"')
    parser.add_argument('-range', type=float, help='[float] optional, range for eucl. dists to be considered '
                                                   'in agreement (pix), default = 3', default=defaults['range'])
    parser.add_argument('-length', type=float, help='[float] optional, value over which dists will not be '
                                                   'considered (pix), default = 100', default=defaults['length'])
    parser.add_argument('-num', type=int, help='[int] optional, # of sources, sorted by mag, to consider in '
                                                   'dist. calculation (dont make too large!), default = 15', default=defaults['num'])
    parser.add_argument('-stdev', type=float, help='[float] # of stdevs away from median to prune eucl. dists '
                                                   'default = 1', default=defaults['stdev'])
    parser.add_argument('-segstd', type=float, help='[float] # of stdevs away from median to prune final shifts '
                                                    '(for use with the -segment flag), default = 2', default=defaults['segstd'])
    args, unknown = parser.parse_known_args()

    shift(args.dir, args.imagename, args.band, args.range, args.length, args.num, args.stdev, args.segstd, args.test,
          args.vary)


if __name__ == "__main__":
    main()