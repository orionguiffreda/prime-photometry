"""
Corrects initial astrometry for zero order translation error
"""

import os
import subprocess
import argparse
import threading
from datetime import datetime as dt

from astroquery.vizier import Vizier
from astropy.io import fits
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.wcs import utils
from astropy.table import Table
import numpy as np
from astropy.table import Column
from itertools import combinations
from collections import defaultdict
import math
import sys

from photometrus.settings import (gen_config_file_name, bulge_checker, CHIP_ZPS, PHOTOMETRY_QUERY_CATALOGS, AB_OFFSET_DICT,
                                  GB_QUERY_CATALOGS, set_vizier_mirror)
from photometrus.utils.defaults import ASTROM_DEFAULTS_NEW as defaults

def get_zp(band):
    zps = [v for k, v in CHIP_ZPS.items() if k == band]
    chosen_zp = zps[0]

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


def imaging(directory, imageName, x_offset=0, x_guess=None, y_guess=None):
    os.chdir(directory)
    f = fits.open(os.path.join(directory, imageName))
    data = f[0].data  # This is the image array
    header = f[0].header

    # strong the image WCS into an object
    w = WCS(header)

    # Get the RA and Dec of the center of the image
    if x_offset != 0:
        x_center = (data.shape[0] / 2) + x_offset
        [raImage, decImage] = w.all_pix2world(x_center, data.shape[1] / 2, 1)
    else:
        [raImage, decImage] = w.all_pix2world(data.shape[0] / 2, data.shape[1] / 2, 1)

    # Get zero point for image
    # chip = header['CHIP']
    if header['FILTER2'] == 'Open':
        band = 'Z'
    else:
        band = header['FILTER2']

    zp = get_zp(band)

    # detect if GB field
    case = header['OBJTYPE']
    bulge = bulge_checker(case)

    # if 'Bulge' in header['OBJTYPE']:
    #     bulge = True
    # else:
    #     bulge = False

    # cat name
    pre = os.path.splitext(imageName)[0]
    ext = os.path.splitext(imageName)[1]
    if 'shift' in pre:
        catname = pre[:-6] + '.cat'
    else:
        catname = pre+'.cat'

    # optional initial shift guess
    if x_guess and y_guess:
        print('Applying initial shift guesses: x = %s, y = %s' % (x_guess, y_guess))
        guess_hdr = header.copy()
        x_init = guess_hdr['CRPIX1']
        y_init = guess_hdr['CRPIX2']
        guess_hdr['CRPIX1'] = x_init + x_guess
        guess_hdr['CRPIX2'] = y_init + y_guess

        pre = pre.replace('flat','init.flat')
        init_imagename = pre+ext
        fits.writeto(os.path.join(directory, init_imagename), data, header, overwrite=True)  # preserve original fits img

        fits.writeto(os.path.join(directory, imageName), data, guess_hdr, overwrite=True)
    elif x_guess or y_guess:
        print('Please provide BOTH an X and Y guess! Proceeding with no initial guess...')


    return data, header, w, raImage, decImage, zp, catname, bulge


#%% catalog query


def cat_query(coords, band, boxsize, catNum, magcol, maglow=12.5, maghigh=14.5, errbits='<=16', bulge=False):
    columns = ['RAJ2000', 'DEJ2000', 'RAICRS', 'DEICRS', 'RA_ICRS', 'DE_ICRS', '%sap3' % band, '%s1ap3' % band,
               'e_%sap3' % band, '%smag' % band, "%smag3" % band, 'e_%smag' % band, '%smag' % band.lower()]

    # chosen coords
    if bulge:
        chosencoords = coords
        coords = chosencoords.galactic    # galactic conversion for bulge fields
        chosen_frame = 'galactic'
        frame_long = coords.l.deg
        frame_long_str = 'l = %.4f' % frame_long
        frame_lat = coords.b.deg
        frame_lat_str = 'b = %.4f' % frame_lat
    else:
        chosen_frame = 'fk5'
        frame_long = coords.ra.deg
        frame_long_str = 'RA: %.4f' % frame_long
        frame_lat = coords.dec.deg
        frame_lat_str = 'DEC: %.4f' % frame_lat

    if errbits != '<=16':
        errbits_vvv = errbits[0]
        errbits_2M = errbits[1]
    else:
        errbits_vvv = errbits
        errbits_2M = '!= null'
    print('Querying Vizier %s around %s, %s, boxwidth %.2f arcmin, mag lim of %s - %s, errbits constraints: %s, %s'
          % (catNum, frame_long_str, frame_lat_str, boxsize, maglow, maghigh, errbits_vvv, errbits_2M))
    v = Vizier(columns=columns, column_filters={"%s" % magcol: "%s .. %s" % (maglow, maghigh),
                                                 "%sperrbits" % band: errbits_vvv,
                                                "%s1perrb" % band: errbits_vvv,
                                                "%sflags" % band: errbits_vvv,
                                                "Cflg": errbits_2M,
                                                "Nd": ">6"}, row_limit=-1)
    Q = v.query_region(SkyCoord(coords, unit=(u.deg, u.deg)), width=str(boxsize) + 'm',
                       catalog=catNum, cache=False, frame=chosen_frame)
    print('Queried source total = ', len(Q[0]))
    return Q


def complex_query(raImage, decImage, band, boxsize, maglow=12, maghigh=14, bulge=False):
    # new automatic survey picking

    # current catalogs
    if bulge:
        print('Galactic bulge field detected!  Adjusting query parameters accordingly...')
        catalog_dict = GB_QUERY_CATALOGS
        catalogs = []
        if band == 'J' or band == 'H' or band == 'Y':
            for k, v in catalog_dict.items():
                if v[0] == 'J':
                    catalogs.append((k, v[1]))
        elif band == 'Z':
            for k, v in catalog_dict.items():
                if v[0] == 'Z':
                    catalogs.append((k, v[1]))
        else:
            print('Only J, H, Y, and Z band are supported!')
        coords = SkyCoord(ra=raImage * u.degree, dec=decImage * u.degree, frame='fk5')
        coords = coords.galactic    # galactic conversion for bulge fields
        chosen_frame = 'galactic'
        frame_long = coords.l.deg
        frame_long_str = 'l = %.4f' % frame_long
        frame_lat = coords.b.deg
        frame_lat_str = 'b = %.4f' % frame_lat
        print('Converting coords to galactic: %s, %s' % (frame_long_str, frame_lat_str))

        mag_high_cutoff = maghigh
        mag_low_cutoff = maglow
    else:
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
        coords = SkyCoord(ra=raImage * u.degree, dec=decImage * u.degree, frame='fk5')
        chosen_frame = 'fk5'
        frame_long = raImage
        frame_long_str = 'RA: %.4f' % frame_long
        frame_lat = decImage
        frame_lat_str = 'DEC: %.4f' % frame_lat

        print('Non-bulge field, implementing wide mag range...')
        mag_high_cutoff = 16
        mag_low_cutoff = maglow

    checkwidth = boxsize

    not_null = '!=null'

    # current columns
    v = Vizier(columns=['RAJ2000', 'DEJ2000', 'RAICRS', 'DEICRS', 'RA_ICRS', 'DE_ICRS', '%sap3' % band, '%s1ap3' % band,
                        '%smag' % band, "%smag3" % band, "%smag1" % band, '%smag' % band.lower(),
                        '%sPSF' % band.lower(), '%spmag' % band.lower()])
               # column_filters={'%sap3' % band: not_null, '%smag' % band: not_null, "%smag3" % band: not_null
               #     ,'%smag' % band.lower(): not_null, '%sPSF' % band.lower(): not_null, '%spmag' % band.lower(): not_null
               # })
    try:
        result = v.query_region(coords, width=str(checkwidth) + 'm', catalog=[f[1] for f in catalogs], frame=chosen_frame)
        test = result[0]
    except IndexError as e:
        # in case of strange failure in query, default to a 2mass query attempt
        try:
            result = v.query_region(coords, width=str(checkwidth) + 'm', catalog='II/246/',
                                    frame=chosen_frame)
        except IndexError as e:
            raise Exception(f'Sadly, no current surveys available in current area in {band} band: {e}')

    keys = result.format_table_list()

    print('%s, %s, Box Width: %s arcmin... '
          '\n%s band initial query resulting in: \n%s' % (frame_long_str, frame_lat_str, checkwidth, band, keys))

    keycheck = result.keys()

    errbitoptions_2mass = ['!= null','~0??', '~?0?']     # 2mass errbits column, 1st is for J and 2nd is for H
    if band == 'J':
        errbit_2mass = errbitoptions_2mass[1]
    elif band == 'H':
        errbit_2mass = errbitoptions_2mass[2]
    else:
        errbit_2mass = errbitoptions_2mass[0]

    if bulge:
        acc_source_num = 200    # total number of sources allowed in full query before trying smaller box
    else:
        acc_source_num = 400

    # contains changes in bounds to iterate through if too many sources:
    # format: [boxsize multiplier, mag lim scalar change, errbits column constraint]
    bounds_change_list = [[1.0, 0, 0, '<=16', '!= null'], [1.0, 0.5, 0, '<=16', '!= null'], [1.0, 0.5, 0, '<16', errbit_2mass],
                          [1.0, 0.5, 0.5, '<16', errbit_2mass], [1.0, 0.5, 0.75, '<16', errbit_2mass],
                          [1.0, 0.5, 1.0, '<16', errbit_2mass], [1.0, 0.5, 1.25, '<16', errbit_2mass],
                          [0.85, 0.5, 1.25, '<16', errbit_2mass]]

    success_flag = False
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

                    # AB surveys
                    if chosen_survey in ['DES_Y', 'DES_Z', 'Skymapper']:
                        print('Adjusting AB mags to VEGA...')
                        offset = AB_OFFSET_DICT.get(band, 0.0)
                        mag_high_cutoff -= offset

                    no_sources_flag = False

                    for bounds in bounds_change_list:
                        boxscale, low_lim_scalar, high_lim_scalar, errbits_constraint, errbits_2M = bounds
                        errbits = [errbits_constraint, errbits_2M]

                        effective_boxsize = boxsize * boxscale
                        eff_mag_low_cutoff = mag_low_cutoff + low_lim_scalar
                        eff_mag_high_cutoff = mag_high_cutoff - high_lim_scalar

                        print('\nQuerying Vizier %s around %s, %s, boxwidth %.2f arcmin, mag lim of %s - %s, '
                              'err constraints: %s, %s'
                              % (catNum, frame_long_str, frame_lat_str, effective_boxsize,
                                 eff_mag_low_cutoff, eff_mag_high_cutoff, errbits_constraint, errbits_2M))
                        try:
                            v = Vizier(columns=[cols[0], cols[1], cols[2]],
                                       column_filters={
                                           cols[2]: f"{eff_mag_low_cutoff:f}..{eff_mag_high_cutoff:f}",
                                           f"{band.lower()}Flag": "<4",
                                           f"{band}perrbits": errbits_constraint,
                                           f"{band}1perrb": errbits_constraint,
                                           f"{band}flags": errbits_constraint,
                                           "Cflg": errbits_2M,
                                           "Hclass": "== -1",
                                           "Class": "== 0",
                                           "Nd": ">6"
                                       }, row_limit=-1)

                            Q = v.query_region(SkyCoord(coords, unit=(u.deg, u.deg)),
                                               width=str(effective_boxsize) + 'm',
                                               catalog=catNum, cache=False,
                                               frame=chosen_frame)

                            if Q and len(Q[0]) > 0:
                                print('Queried source total = ', len(Q[0]))
                                if len(Q[0]) <= acc_source_num:
                                    success_flag = True  # Mark success
                                    break
                                else:
                                    print("Too many sources (>%i), trying different bounds..." % acc_source_num)
                                    no_sources_flag = True
                                    continue
                            else:
                                print(f"No sources found in {catNum}, trying fallback if available...")
                                break

                        except Exception as e:
                            print('Error in Vizier query.')
                            print(f"Error details: {e}")
                            continue

                    if success_flag:
                        break  # Break out of keycheck loop as well

        if success_flag:
            break  # Break out of catalogs loop

    return Q, coords, catNum, cols[2], eff_mag_low_cutoff, eff_mag_high_cutoff, effective_boxsize, errbits

#%% sextraction / psfex


def sex1(imageName, bulge=False):
    print('Running sextractor for psf...')
    if bulge:
        configFile = gen_config_file_name('sex2.config')
    else:
        configFile = gen_config_file_name('bulge_new.config')

    paramName = gen_config_file_name('astromshift_new.param')

    if imageName.endswith('.fits') or imageName.endswith('.new'):
        catname = os.path.splitext(imageName)[0] + '.cat'
    else:
        catname = imageName + '.cat'
    try:
        command = ('sex %s -c %s -CATALOG_NAME %s -PARAMETERS_NAME %s -CHECKIMAGE_TYPE NONE'
                   % (imageName, configFile, catname, paramName))
        #print('Executing command: %s' % command)
        rval = subprocess.run(command.split(), check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError as err:
        print('Could not run sextractor with exit error %s'%err)
    return catname


def psfex(catalogName):
    print('Running PSFex on sextrctr catalogue to generate psf for stars in the img...')
    psfConfigFile = gen_config_file_name('default.psfex')
    try:
        command = 'psfex %s -c %s' % (catalogName, psfConfigFile)
        # print('Executing command: %s' % command)
        subprocess.run(command.split(), check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError as err:
        print('Could not run psfex with exit error %s' % err)


def sex2(imageName):
    print('Running sextractor for psf w/ psfex model fits...')
    psfName = os.path.splitext(imageName)[0] + '.psf'
    newcatalogName = os.path.splitext(imageName)[0] + '.psf.cat'

    configFile = gen_config_file_name('sex2.config')
    paramName = gen_config_file_name('adv_shift.param')

    # We are supplying SExtactor with the PSF model with the PSF_NAME option
    command = 'sex %s -c %s -CATALOG_NAME %s -PSF_NAME %s -PHOT_FLUXFRAC 0.5 -PARAMETERS_NAME %s' % (
        imageName, configFile, newcatalogName, psfName, paramName)
    # print("Executing command: %s" % command)
    subprocess.run(command.split(), check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    return newcatalogName


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


def make_tables(directory, data, w, catname, Q, band, crop, zp, maglow=12, maghigh=14, shifted_cat=None, adv=False):
    print('Creating sorted catalog & prime tables...')
    if shifted_cat:
        sexcat = shifted_cat
    else:
        sexcat = Table.read(os.path.join(directory, catname), hdu=2)
    colnames = Q[0].colnames

    if isinstance(crop, float):
        left_crop = right_crop = y_crop = crop
    else:
        left_crop = crop[0]
        right_crop = crop[1]
        y_crop = crop[2]

    max_x = data.shape[0]
    max_y = data.shape[1]
    mass_imCoords = w.all_world2pix(Q[0][colnames[0]], Q[0][colnames[1]], 1)
    inner_catsources = Q[0][np.where((mass_imCoords[0] > left_crop) & (mass_imCoords[0] < (max_x-right_crop)) & (mass_imCoords[1] > y_crop)
                                   & (mass_imCoords[1] < (max_y-y_crop)))]
    # inner_catsources = inner_catsources[(inner_catsources[colnames[2]] * 0.339 > 1)]

    # print(sexcat.colnames)
    if adv:
        inner_primesources = sexcat[(sexcat['FLAGS'] < 2) &
                                    (sexcat['X_IMAGE'] < (max_x - right_crop)) & (sexcat['X_IMAGE'] > left_crop)
                                    & (sexcat['Y_IMAGE'] < (max_y) - y_crop) & (
                                                sexcat['Y_IMAGE'] > y_crop) & (sexcat['FLUX_RADIUS'] * 0.498 > 1)
                                    & (sexcat['CLASS_STAR'] > 0.5)]
    else:
        inner_primesources = sexcat[(sexcat['FLAGS'] <= 1) &
                                    (sexcat['X_IMAGE'] < (max_x - right_crop)) & (sexcat['X_IMAGE'] > left_crop)
                                    & (sexcat['Y_IMAGE'] < (max_y) - y_crop) & (
                                                sexcat['Y_IMAGE'] > y_crop) & (sexcat['FLUX_RADIUS'] * 0.498 > 1)
                                    & (sexcat['FLUX_MAX'] / sexcat['FLUX_AUTO'] < 0.25)
                                    ]

    # zp correction
    print('applying zp correction to PRIME mags: %s' % zp)
    inner_primesources['MAG_AUTO'] = inner_primesources['MAG_AUTO'] + zp

    # inner_primesources.write('prime_all.ecsv', overwrite=True)
    inner_primesources = inner_primesources[(inner_primesources['MAG_AUTO'] >= maglow) &
                                            (inner_primesources['MAG_AUTO'] <= maghigh)]

    # inner_primesources.sort('MAG_AUTO')
    # inner_catsources.sort(colnames[2])
    cat_xy = w.all_world2pix(inner_catsources[colnames[0]], inner_catsources[colnames[1]], 1)

    xs = Column(cat_xy[0], name='X_IMAGE', unit='pix')
    ys = Column(cat_xy[1], name='Y_IMAGE', unit='pix')
    inner_catsources.add_column(xs)
    inner_catsources.add_column(ys)

    inner_primesources.write('prime_all.ecsv', overwrite=True)

    Path = os.path.join(directory, 'catcoords_crop.reg')
    newtext = open(Path, 'w+')
    for i,j in zip(xs,ys):
        newtext.write('\npoint(%f,%f) # point=circle 5' % (i,j))
    # #
    Path = os.path.join(directory, 'primecoords_crop.reg')
    newtext = open(Path, 'w+')
    for i,j,k in zip(inner_primesources['X_IMAGE'],inner_primesources['Y_IMAGE'],inner_primesources['FLUX_RADIUS']):
        newtext.write(f'\ncircle({i}, {j}, {k}") # color=red')
    #
    # Path = os.path.join(directory, 'primecoords_all.reg')
    # newtext = open(Path, 'w+')
    # for i,j in zip(sexcat['X_IMAGE'],sexcat['Y_IMAGE']):
    #     newtext.write('\npoint(%f,%f) # point=circle 5' % (i,j))

    return inner_primesources, inner_catsources, colnames


def prep_tables(inner_primesources,inner_catsources,num):
    first_primes_x = np.array(inner_primesources['X_IMAGE'][:num])
    first_primes_y = np.array(inner_primesources['Y_IMAGE'][:num])
    first_cats_x = np.array(inner_catsources['X_IMAGE'][:num])
    first_cats_y = np.array(inner_catsources['Y_IMAGE'][:num])

    first_primecoords = np.column_stack([first_primes_x,first_primes_y])
    first_catcoords = np.column_stack([first_cats_x,first_cats_y])

    # Path = '/mnt/photometry/supermaster_test/GB63_dither_20250406/J/C1_sub/catcoords_prune.reg'
    # #coords_prime = coords_prime.T
    # newtext = open(Path, 'w+')
    # for i,j in zip(first_cats_x,first_cats_y):
    #     newtext.write('\npoint(%f,%f) # point=circle 5' % (i,j))
    #
    # Path = '/mnt/photometry/supermaster_test/GB63_dither_20250406/J/C1_sub/primecoords_prune.reg'
    # #coords_prime = coords_prime.T
    # newtext = open(Path, 'w+')
    # for i,j in zip(first_primes_x,first_primes_y):
    #     newtext.write('\npoint(%f,%f) # point=circle 5' % (i,j))

    return first_primecoords,first_catcoords

#%% distance calc


def calculate_distance(p1, p2):
    return np.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2)


def find_agreeing_distances(table1, table2, length):
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

    # Find distances under certain length
    agreeing_pairs = []

    for d in sorted(distances.keys()):
        if d <= length:
            agreeing_pairs.extend(distances[d])

    # Remove duplicates
    agreeing_pairs = list(set(agreeing_pairs))

    # Extract distances
    agreeing_distances = [pair for pair in row_pairs if pair[0] in agreeing_pairs]

    agreeing_pairs.sort()
    # print(f"\nEuclidean distances that are < {length} pixels in length:")
    # for pair, dist in zip(agreeing_pairs[:20], agreeing_distances[:20]):
    #     print(f"Row from prime: {pair[0]}, Row from catalog: {pair[1]}, Distance (pix): {dist[1]}")

    return agreeing_pairs, agreeing_distances


#%% final shift calc and header update


def xyshifts(pairs, inner_primesources, inner_catsources, iters):
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

    x_shiftarr = np.array(x_shiftarr)
    y_shiftarr = np.array(y_shiftarr)
    # array of all considered x and y shifts

    # histograms of binned shifts for x and y (1 pix bins)
    print('Binning and sorting all considered x and y shifts...')
    binwidth = 1
    x_bins = np.arange(min(x_shiftarr), max(x_shiftarr) + binwidth, binwidth)
    y_bins = np.arange(min(y_shiftarr), max(y_shiftarr) + binwidth, binwidth)

    x_counts, x_edges = np.histogram(x_shiftarr, bins=x_bins, density=True)
    y_counts, y_edges = np.histogram(y_shiftarr, bins=y_bins, density=True)

    # assign bin index
    x_indices = np.digitize(x_shiftarr, bins=x_edges) - 1
    y_indices = np.digitize(y_shiftarr, bins=y_edges) - 1

    # clip indices to valid range (in case some fall on rightmost edge)
    x_indices = np.clip(x_indices, 0, len(x_counts) - 1)
    y_indices = np.clip(y_indices, 0, len(y_counts) - 1)

    # map each shift value to its bin count
    x_bin_popularity = x_counts[x_indices]
    y_bin_popularity = y_counts[y_indices]

    # sort original arrays by bin popularity (descending)
    x_shifts_sorted = x_shiftarr[np.argsort(-x_bin_popularity)]
    x_shifts_sorted_pop = sorted(x_bin_popularity, reverse=True)

    y_shifts_sorted = y_shiftarr[np.argsort(-y_bin_popularity)]
    y_shifts_sorted_pop = sorted(y_bin_popularity, reverse=True)

    xy_shifts = np.column_stack([x_shifts_sorted, y_shifts_sorted])
    xy_shift_pops = np.column_stack([x_shifts_sorted_pop, y_shifts_sorted_pop])

    print(f'After binning and sorting, top {iters} shifts considered: ')
    for pt, pop in zip(xy_shifts[:iters], xy_shift_pops[:iters]):
        print('X offset = %.3f, Y offset = %.3f, X Bin Norm. Popularity = %.3f, Y Bin Norm. Popularity = %.3f' %
              (pt[0], pt[1], pop[0], pop[1]))

    # crpix1 = header['CRPIX1']
    # crpix2 = header['CRPIX2']
    #
    # header['CRPIX1'] = crpix1 + xfinal_shift
    # header['CRPIX2'] = crpix2 + yfinal_shift
    #
    # imageshiftname = os.path.splitext(imageName)[0]
    # imageshiftname = imageshiftname + '.shift.fits'
    #
    # print('Writing new FITS file w/ updated CRPIX: %s' % imageshiftname)
    # newpath = os.path.join(directory, imageshiftname)
    # oldpath = os.path.join(directory, imageName)
    #
    # fits.writeto(newpath, data, header, overwrite=True)

    return xy_shifts


#%% Apply shifts quickly to all sources (removes necessity for multiple sextractor runs)
def apply_shifts_to_cat(directory, catname, x_shift, y_shift):
    primecat = Table.read(os.path.join(directory, catname), hdu=2)
    shifted_inner_primesources = primecat.copy()
    # print(shifted_inner_primesources)

    # shifted_inner_primesources['X_IMAGE'] += x_shift
    # shifted_inner_primesources['Y_IMAGE'] += y_shift
    # print(shifted_inner_primesources)
    return shifted_inner_primesources


#%%  Rerunning sextractor and astroquery for new source positions
def shiftiteration(directory, imagename, filter_used, coords, maglow, maghigh, eff_boxsize, crop, catNum, magcol, errbits,
                   x_shift, y_shift, catname, adv
                   ):
    data, header, w, raImage, decImage, zp, catname_old, bulge = imaging(directory, imagename)

    # sex1(imagename, bulge=bulge)
    Q = cat_query(coords, filter_used, eff_boxsize, catNum, magcol, maglow, maghigh, errbits, bulge=bulge)
    # Q = complex_query(raImage, decImage, filter_used, boxsize, maglow=maglow, maghigh=maghigh, bulge=bulge)
    shifted_primecat = apply_shifts_to_cat(directory, catname, x_shift, y_shift)
    inner_primesources_iter, inner_catsources_iter, colnames = make_tables(directory, data, w, catname, Q, filter_used, crop, zp,
                                                                 maglow=maglow, maghigh=maghigh, shifted_cat=shifted_primecat,
                                                                           adv=adv)

    return inner_primesources_iter, inner_catsources_iter, colnames, w


#%% iterate through most likely binned shifts and check correctness w/ crossmatch

def iterate_and_test(
        xy_shifts, directory, header, data, imageName, filter_used, eff_boxsize, crop, coords, catNum, magcol,
        crsmtch_thresh_low, crsmtch_thresh_high, crsmtch_iters, maglow, maghigh, errbits, catname, bulge=False, adv=False
):
    print('Beginning iterative testing of sorted shifts...')
    best_completion = -1
    best_shift = None
    max_iterations = crsmtch_iters

    for idx, shift_pair in enumerate(xy_shifts[:max_iterations]):
        x_shift, y_shift = shift_pair

        # Copy header to avoid cumulative modification
        header_copy = header.copy()
        header_copy['CRPIX1'] = header['CRPIX1'] + x_shift
        header_copy['CRPIX2'] = header['CRPIX2'] + y_shift

        imageshiftname = os.path.splitext(imageName)[0] + '.shift.fits'
        newpath = os.path.join(directory, imageshiftname)

        print(f'\nTesting shift [{idx + 1}/{crsmtch_iters}]: Writing temporary new FITS file w/ updated CRPIX...')
        fits.writeto(newpath, data, header_copy, overwrite=True)

        # Rerunning sextractor and astroquery for new source positions
        inner_primesources_iter, inner_catsources_iter, colnames, wcs = (
            shiftiteration(directory, imageshiftname, filter_used, coords, maglow, maghigh, eff_boxsize, crop, catNum,
                           magcol, errbits, x_shift=x_shift, y_shift=y_shift, catname=catname, adv=adv))
        # Run crossmatch to determine successful solve
        SourceCatCoords = SkyCoord(ra=inner_catsources_iter[colnames[0]], dec=inner_catsources_iter[colnames[1]],
                                   frame='icrs', unit='degree')
        SourcePrimeCoords = utils.pixel_to_skycoord(inner_primesources_iter['X_IMAGE'], inner_primesources_iter['Y_IMAGE'],
                                                    wcs, origin=1)
        # SourcePrimeCoords = SkyCoord(ra=inner_primesources_iter['ALPHA_J2000'], dec=inner_primesources_iter['DELTA_J2000'],
        #                          frame='icrs', unit='degree')
        if bulge:
            photoDistThresh = 0.6 * u.arcsec    # looser crossmatch dist for all_sky
        else:
            photoDistThresh = 1.0 * u.arcsec
        idx_prime, idx_cat, d2d, d3d = SourceCatCoords.search_around_sky(SourcePrimeCoords, photoDistThresh)

        # Value to measure crossmatch completion, if high enough, then should be a successful solve

        print('PRIME Source Num = ',len(inner_primesources_iter))
        print('Catalog Source Num = ',len(inner_catsources_iter))
        print('Crossmatch Source Num = ',len(idx_prime))
        if len(inner_primesources_iter) > len(inner_catsources_iter):
            crsmtch_completion = len(idx_cat) / len(inner_catsources_iter)
        else:
            crsmtch_completion = len(idx_prime) / len(inner_primesources_iter)
        print(f'({round(x_shift,5)}, {round(y_shift,5)}) -> Crossmatch completion: {crsmtch_completion:.2f}')

        # Update best result if this run is better
        if crsmtch_completion > best_completion:
            best_completion = crsmtch_completion
            best_shift = shift_pair

        # Stop early if successful enough
        if crsmtch_completion > crsmtch_thresh_high:
            print('Sufficient match found!, stopping iteration')
            print(f'Final chosen shifts: ({round(best_shift[0],5)}, {round(best_shift[1],5)})')
            return shift_pair[0], shift_pair[1]

    # Final check: warn if best match is too low
    if best_completion < crsmtch_thresh_low:
        print('WARNING: No good match found. Best crossmatch completion was < 0.1.')
        return 0, 0

    print(f'Finished {crsmtch_iters} attempts. Best and final shifts: ({round(best_shift[0],5)}, {round(best_shift[1],5)}) '
          f'→ Completion: {best_completion:.2f}')
    return best_shift[0], best_shift[1]



def change_all_files(xfinal_shift, yfinal_shift, directory):
    if xfinal_shift == 0 or yfinal_shift == 0:
        print('No agreement, thus cannot move forward with rewriting all files!')
    else:
        all_fits = [f for f in sorted(os.listdir(directory)) if f.endswith('.flat.new')]
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

            header.set('X_SHIFT', xfinal_shift, 'X Value CRPIX1 shift', after='WCSAXES')
            header.set('Y_SHIFT', yfinal_shift, 'Y Value CRPIX2 shift', after='X_SHIFT')
            header.set('CRPIX1_OLD', crpix1, 'Initial CRPIX1 Value', after='Y_SHIFT')
            header.set('CRPIX2_OLD', crpix2, 'Initial CRPIX2 Value', after='CRPIX1_OLD')

            imageshiftname = os.path.splitext(f)[0]
            imageshiftname = imageshiftname + '.shift.new'
            newpath = os.path.join(directory, imageshiftname)

            fits.writeto(newpath, data, header, overwrite=True)

        print('Moving old files to %s/old/ directory and renaming shifted images...' % directory)
        old_storage_dir = os.path.join(directory, 'old')
        if not os.path.exists(old_storage_dir):
            os.mkdir(old_storage_dir)

        all_fits_again = [f for f in sorted(os.listdir(directory)) if f.endswith('.flat.new')]
        all_fits_shift = [f for f in sorted(os.listdir(directory)) if f.endswith('.shift.new')]

        for f in all_fits_again:
            currentpath = os.path.join(directory, f)
            oldpath = os.path.join(old_storage_dir, f)
            os.rename(currentpath, oldpath)

        catfiles = [f for f in sorted(os.listdir(directory)) if f.endswith('.flat.cat') or f.endswith('.psf.cat') or
                    f.endswith('.reg')]
        for cat in catfiles:
            currentpath = os.path.join(directory, cat)
            oldpath = os.path.join(old_storage_dir, cat)
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
    fnames = ['.shift.cat','.shift.fits','.psf','.psf.cat','.reg']
    print('Removing intermediate files')
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
        print('No intermediate files found to remove')
#%%


#%%

def boxchange(size, x_offset=0):
    # size = size of box in arcmin
    pixscale = 0.498
    imgsize = 4088

    boxsize_pix = size * 60 * (1/pixscale)
    crop_pix = (imgsize - boxsize_pix) / 2

    if x_offset != 0:
        center_pix = imgsize / 2

        new_center_x = center_pix + x_offset

        x_start = new_center_x - boxsize_pix / 2
        x_end = new_center_x + boxsize_pix / 2

        left_crop = x_start  # pixels from left edge
        right_crop = imgsize - x_end  # pixels from right edge
        crops = [left_crop, right_crop, crop_pix]
        return size, crops
    return size, crop_pix


def shift(
        directory, imagename, band, length=defaults['length'], num=defaults['num'], thresh_low=defaults['thresh_low'],
        thresh_high=defaults['thresh_high'], iters=defaults['iters'], test=False, adv_solve=defaults['adv_solve'],
        x=defaults['x_guess'],y=defaults['y_guess']
):

    x_offset = 0
    if band == 'Y':
        print('Switching Y band to J for ease of astrometry...')
        filter_used = 'J'
    else:
        filter_used = band

    start_time = dt.now()

    set_vizier_mirror()

    maglow = 12
    maghigh = 14

    data, header, w, raImage, decImage, zp, catname, bulge = imaging(directory, imagename, x_offset=x_offset, x_guess=x, y_guess=y)
    if bulge:
        boxsize, crop = boxchange(4)
        thresh_high = 0.25
    else:
        boxsize, crop = boxchange(10, x_offset=x_offset)

    if adv_solve:
        if os.path.isfile(os.path.splitext(imagename)[0] + '.psf.cat'):
            print('Previous psf cat detected, saving you some time and skipping all sextractor steps...')
            catname = os.path.splitext(imagename)[0] + '.psf.cat'
        else:
            catname = sex1(imagename, bulge=bulge)
            psfex(catname)
            catname = sex2(imagename)
    else:
        sex1(imagename, bulge=bulge)

    # Q = cat_query(raImage, decImage, filter_used, boxsize, maglow=maglow, maghigh=maghigh)
    Q, coords, catNum, magcol, mag_low_cutoff, mag_high_cutoff, eff_boxsize, errbits = complex_query(raImage, decImage,
                                                                                                     filter_used, boxsize,
                                                                         maglow=maglow, maghigh=maghigh, bulge=bulge)
    if eff_boxsize != boxsize:
        print('Adjusting crop for crossmatching with %.1f boxsize...' % eff_boxsize)
        eff_boxsize, crop = boxchange(eff_boxsize)

    inner_primesources, inner_catsources, colnames = make_tables(directory, data, w, catname, Q, filter_used, crop, zp,
                                                                 maglow=mag_low_cutoff, maghigh=mag_high_cutoff, adv=adv_solve)
    print('PRIME Source Num = ', len(inner_primesources))
    print('Catalog Source Num = ', len(inner_catsources))

    source_num_ratio = (np.min([len(inner_primesources), len(inner_catsources)]) /
                        np.max([len(inner_primesources), len(inner_catsources)]))
    if abs(source_num_ratio) < 0.1:
        print(f'\nLarge disparity in source number betw. PRIME & survey, min/max ratio: {round(source_num_ratio,3)}, '
              f'\nShift very likely to be inaccurate, assuming 0 shifts & skipping this image...')
        ultimate_shift_x = ultimate_shift_y = 0

    else:
        first_primecoords, first_catcoords = prep_tables(inner_primesources, inner_catsources, num)
        agreeing_pairs, dists = find_agreeing_distances(first_primecoords, first_catcoords, length)
        try:
            xy_shifts = xyshifts(agreeing_pairs, inner_primesources, inner_catsources, iters)

            ultimate_shift_x, ultimate_shift_y = iterate_and_test(xy_shifts, directory, header, data, imagename, filter_used,
                                                                  eff_boxsize, crop, coords, catNum, magcol, thresh_low, thresh_high,
                                                                  iters, maglow=mag_low_cutoff, maghigh=mag_high_cutoff, errbits=errbits,
                                                                  catname=catname, bulge=bulge, adv=adv_solve)
        except ValueError as e:
            print(f'Error in generating xy shifts, assuming 0, error: {e}')
            ultimate_shift_x = ultimate_shift_y = 0
        except Exception as e:
            print(f'Error in final shift generation, assuming 0, error: {e}')
            ultimate_shift_x = ultimate_shift_y = 0

    if not test:
        change_all_files(ultimate_shift_x, ultimate_shift_y, directory)
        if ultimate_shift_x != 0:
            removal(directory)

    end_time = dt.now()
    print('astrometric shift correction time:', (end_time - start_time).total_seconds())

    return ultimate_shift_x, ultimate_shift_y


def main():
    parser = argparse.ArgumentParser(description='Corrects for translation in initial astrometry '
                                                 '(so astrom.net doesnt need to be used)')
    parser.add_argument('-test', action='store_true', help='optional flag to test for a successful solve, '
                                                           'doesnt apply to all files in directory.', default=defaults['test'])
    parser.add_argument('-adv_solve', action='store_true', help='optional flag to run psfex and sextractor again for max'
                                                                ' astrometric & mag accuracy, designed for very dense '
                                                                'fields, but takes a long time.', default=defaults['adv_solve'])
    parser.add_argument('-dir', type=str, help='[str] path where input file is stored (should run on proc. image, '
                                               'so likely should be /C#_sub/)')
    parser.add_argument('-imagename', type=str, help='[str] input file name (should run on proc. image, '
                                                     'i.e. *.sky.flat.fits)')
    parser.add_argument('-band', type=str, help='[str] filter used, ex. "J"')
    parser.add_argument('-x', type=float, help='[float] optional initial guess for X, '
                                               '*PROVIDE BOTH AN X AND Y GUESS USING -x & -y*', default=defaults['x_guess'])
    parser.add_argument('-y', type=float, help='[float] optional initial guess for Y, '
                                               '*PROVIDE BOTH AN X AND Y GUESS USING -x & -y*', default=defaults['y_guess'])
    parser.add_argument('-length', type=float, help='[float] optional, value over which dists will not be '
                                                   'considered (pix), default = 100', default=defaults['length'])
    parser.add_argument('-num', type=int, help='[int] optional, # of sources, sorted by mag, to consider in '
                                                   'dist. calculation, default = 400', default=defaults['num'])
    parser.add_argument('-thresh_low', type=float, help='[float] Failure threshold for crossmatch betw. '
                                                    'PRIME & astroquery catalog, default = 0.1', default=defaults['thresh_low'])
    parser.add_argument('-thresh_high', type=float, help='[float] Successful threshold for crossmatch betw. '
                                                    'PRIME & astroquery catalog, default = 0.4', default=defaults['thresh_high'])
    parser.add_argument('-iters', type=int, help='[int] Number of iterations to test shift values before'
                                                 ' defaulting to highest crossmatch percentage, default = 10', default=defaults['iters'])
    args, unknown = parser.parse_known_args()

    shift(args.dir, args.imagename, args.band, args.length, args.num, args.thresh_low, args.thresh_high, args.iters,
          args.test, args.adv_solve, args.x, args.y)


if __name__ == "__main__":
    main()