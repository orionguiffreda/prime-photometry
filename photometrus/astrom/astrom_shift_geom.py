"""
Corrects initial astrometry for zero order translation error
"""

import os
import subprocess
import requests
import argparse
import threading
import signal
from datetime import datetime as dt

import astroalign as aa
from astroquery.vizier import Vizier
from astropy.io import fits
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.wcs import utils
from astropy.table import Table
import numpy as np
from astropy.table import Column
from scipy.spatial import cKDTree
from itertools import combinations
from collections import defaultdict
import math
import sys

from photometrus.query import query, twomass_query
from photometrus.settings import (gen_config_file_name, bulge_checker, CHIP_ZPS, ASTROM_QUERY_CATALOGS,
                                  AB_OFFSET_DICT, GB_QUERY_CATALOGS, AB_SURVEYS, set_vizier_mirror, local_query_box)
from photometrus.utils.defaults import ASTROM_DEFAULTS_NEW as defaults


class _AlignTimeout(Exception):
    pass


def get_zp(band):
    zps = [v for k, v in CHIP_ZPS.items() if k == band]
    chosen_zp = zps[0]

    return chosen_zp


# %% image info


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
        catname = pre + '.cat'

    # optional initial shift guess
    if x_guess and y_guess:
        print('Applying initial shift guesses: x = %s, y = %s' % (x_guess, y_guess))
        guess_hdr = header.copy()
        w = WCS(guess_hdr)

        crpix1 = guess_hdr['CRPIX1']
        crpix2 = guess_hdr['CRPIX2']

        guess_crval1, guess_crval2 = w.all_pix2world(
            [[crpix1 - x_guess, crpix2 - y_guess]], 1
        )[0]

        guess_hdr['CRVAL1'] = guess_crval1
        guess_hdr['CRVAL2'] = guess_crval2

        pre = pre.replace('flat', 'init.flat')
        init_imagename = pre + ext
        fits.writeto(os.path.join(directory, init_imagename), data, header,
                     overwrite=True)  # preserve original fits img

        fits.writeto(os.path.join(directory, imageName), data, guess_hdr, overwrite=True)
    elif x_guess or y_guess:
        print('Please provide BOTH an X and Y guess! Proceeding with no initial guess...')

    return data, header, w, raImage, decImage, zp, catname, bulge


# %% catalog query


def cat_query(coords, band, boxsize, catNum, magcol, maglow=12.5, maghigh=14.5, errbits='<=16', bulge=False):
    columns = ['RAJ2000', 'DEJ2000', 'RAICRS', 'DEICRS', 'RA_ICRS', 'DE_ICRS', '%sap3' % band, '%s1ap3' % band,
               'e_%sap3' % band, '%smag' % band, "%smag3" % band, 'e_%smag' % band, '%smag' % band.lower()]
    # chosen coords
    if bulge:
        chosencoords = coords
        coords = chosencoords.galactic  # galactic conversion for bulge fields
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

    try:
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
    except (requests.exceptions.ConnectionError, requests.exceptions.Timeout) as e:

        Q = local_query_box(ra_center=float(frame_long),
                            dec_center=float(frame_lat),
                            width=str(boxsize) + 'm',
                            dbname='prime_2mass_local',
                            tablename='twomass_local',
                            columns=["ra", "dec", f"{band.lower()}mag", f"e_{band.lower()}mag"],
                            column_filters={
                                f"{band.lower()}mag": f"BETWEEN {maglow:f} AND {maghigh:f}",
                                "cc_flg": errbits_2M,
                            }
                            )
    print('Queried source total = ', len(Q[0]))
    return Q


def complex_query(raImage, decImage, band, boxsize, maglow=12, maghigh=14, bulge=False):
    TWOMASS_CAT = 'II/246/'

    # ---- catalog selection -------------------------------------------------
    if bulge:
        catalog_dict = GB_QUERY_CATALOGS
        # In the bulge, Y is served from the J-band catalog list
        band_to_catband = {'J': 'J', 'H': 'J', 'Y': 'J', 'Z': 'Z'}
    else:
        catalog_dict = ASTROM_QUERY_CATALOGS
        band_to_catband = {'J': 'J', 'H': 'J', 'Y': 'Y', 'Z': 'Z'}

    if band not in band_to_catband:
        # Original only printed here and then crashed on catalogs[-1]
        raise ValueError('Only J, H, Y, and Z band are supported! Got %r' % (band,))

    catband = band_to_catband[band]
    catalogs = [(k, v[1]) for k, v in catalog_dict.items() if v[0] == catband]
    if not catalogs:
        raise ValueError('No catalogs configured for band %s (bulge=%s)' % (band, bulge))

    # ---- coordinates / magnitude limits ------------------------------------
    coords = SkyCoord(ra=raImage * u.degree, dec=decImage * u.degree, frame='fk5')
    if bulge:
        print('Galactic bulge field detected!  Adjusting query parameters accordingly...')
        coords = coords.galactic  # galactic conversion for bulge fields
        frame_long_str = 'l = %.4f' % coords.l.deg
        frame_lat_str = 'b = %.4f' % coords.b.deg
        print('Converting coords to galactic: %s, %s' % (frame_long_str, frame_lat_str))
        mag_high_cutoff = maghigh
        mag_low_cutoff = maglow
        acc_source_num = 500    # total sources allowed before trying a tighter query
    else:
        frame_long_str = 'RA: %.4f' % raImage
        frame_lat_str = 'DEC: %.4f' % decImage
        print('Non-bulge field, implementing wide mag range...')
        mag_high_cutoff = 16    # NB: maghigh is ignored for non-bulge fields
        mag_low_cutoff = maglow
        acc_source_num = 1500

    # 2MASS quality-flag options: [not-null, J-band, H-band]
    errbitoptions_2mass = ['!= null', '0??', '?0?']
    if band == 'J':
        errbit_2mass = errbitoptions_2mass[1]
    elif band == 'H':
        errbit_2mass = errbitoptions_2mass[2]
    else:
        errbit_2mass = errbitoptions_2mass[0]

    # format: [boxsize multiplier, mag-low shift, mag-high shift,
    #          errbits constraint, 2MASS errbits constraint]
    bounds_change_list = [[1.0, 0, 0, '<=16', '!= null'],
                          [1.0, 0.5, 0, '<=16', '!= null'],
                          [1.0, 0.5, 0, '<16', errbit_2mass],
                          [1.0, 0.5, 0.5, '<16', errbit_2mass],
                          [1.0, 0.5, 0.75, '<16', errbit_2mass],
                          [1.0, 0.5, 1.0, '<16', errbit_2mass],
                          [1.0, 0.5, 1.25, '<16', errbit_2mass],
                          [0.85, 0.5, 1.25, '<16', errbit_2mass]]

    # ---- query loop --------------------------------------------------------
    fallback = None  # last non-empty result that had too many sources

    for chosen_survey, catNum in catalogs:
        cat_mag_low = mag_low_cutoff
        cat_mag_high = mag_high_cutoff

        if chosen_survey in AB_SURVEYS:
            print('Adjusting AB mags to VEGA...')
            offset = AB_OFFSET_DICT.get(band, 0.0)
            cat_mag_high -= offset

        # Non-2MASS catalogs always use the tighter source limit
        max_sources = acc_source_num # if catNum == TWOMASS_CAT else 500

        for boxscale, low_shift, high_shift, errbits_constraint, errbits_2M in bounds_change_list:
            errbits = [errbits_constraint, errbits_2M]
            effective_boxsize = boxsize * boxscale
            eff_mag_low_cutoff = cat_mag_low + low_shift
            eff_mag_high_cutoff = cat_mag_high - high_shift

            print('\nQuerying %s around %s, %s, boxwidth %.2f arcmin, mag lim of %s - %s, '
                  'err constraints: %s, %s'
                  % (catNum, frame_long_str, frame_lat_str, effective_boxsize,
                     eff_mag_low_cutoff, eff_mag_high_cutoff, errbits_constraint, errbits_2M))
            try:
                if catNum == TWOMASS_CAT:
                    Q, _ = twomass_query(
                        coords=coords,
                        frame_long_str=frame_long_str,
                        frame_lat_str=frame_lat_str,
                        band=band,
                        width=effective_boxsize,
                        mag_low_cutoff=eff_mag_low_cutoff,
                        mag_high_cutoff=eff_mag_high_cutoff,
                        errbits=errbits_2M)
                    print('Queried source total = ', len(Q[0]))
                else:
                    chosen_survey = chosen_survey.split('_')[0]
                    Q, _, _ = query(
                        raImage=raImage,
                        decImage=decImage,
                        band=band,
                        survey=chosen_survey,
                        query_width=effective_boxsize,
                        mag_lower_lim=eff_mag_low_cutoff,
                        mag_upper_lim=eff_mag_high_cutoff,
                        errbits=errbits_constraint)

                if Q is None or len(Q) == 0 or len(Q[0]) == 0:
                    print(f"No sources found in {catNum}, trying fallback if available...")
                    break

                nsrc = len(Q[0])
                cols = Q[0].colnames
                result = (Q, coords, catNum, cols[2], eff_mag_low_cutoff,
                          eff_mag_high_cutoff, effective_boxsize, errbits)

                if nsrc <= max_sources:
                    return result

                print("Too many sources (>%i), trying different bounds..." % max_sources)
                fallback = result

            except Exception as e:
                print('Error in query.')
                print(f"Error details: {e}")
                continue

    if fallback is not None:
        print('WARNING: no catalog met the source-count limit; returning the last '
              'over-limit result.')
        return fallback

    raise RuntimeError('No sources found in any catalog for band %s at (%.4f, %.4f)'
                       % (band, raImage, decImage))


# %% sextraction / psfex


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
        # print('Executing command: %s' % command)
        rval = subprocess.run(command.split(), check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as err:
        print(f"SExtractor failed with exit code {err.returncode}")
        print(f"STDERR:\n{err.stderr}")
        print(f"STDOUT:\n{err.stdout}")
        raise Exception('Sextractor failed to run, is the stacked image quality adequate?')
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


def _xy_array(tbl, sort_col=None, is_mag=False):
    """(N,2) array of X_IMAGE/Y_IMAGE, brightest first.
    astroalign takes the first `max_control_points` of a point array as-is,
    so the ordering matters."""
    xy = np.column_stack([np.asarray(tbl['X_IMAGE'], dtype=float),
                          np.asarray(tbl['Y_IMAGE'], dtype=float)])
    if sort_col is not None and sort_col in tbl.colnames:
        vals = np.asarray(tbl[sort_col], dtype=float)
        order = np.argsort(vals if is_mag else -vals)
        xy = xy[order]
    return xy


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

    def cut_report(sexcat, zp, x_lo, x_hi, y_lo, y_hi, maglow, maghigh, peak=0.25):
        mag = sexcat['MAG_AUTO'] + zp
        cuts = [
            ('FLAGS <= 1', sexcat['FLAGS'] <= 1),
            ('inside crop', (sexcat['X_IMAGE'] > x_lo) & (sexcat['X_IMAGE'] < x_hi) &
             (sexcat['Y_IMAGE'] > y_lo) & (sexcat['Y_IMAGE'] < y_hi)),
            ('size > 1"', sexcat['FLUX_RADIUS'] * 0.498 > 1),
            (f'peak/flux < {peak}', sexcat['FLUX_MAX'] / sexcat['FLUX_AUTO'] < peak),
            (f'{maglow} <= mag <= {maghigh}', (mag >= maglow) & (mag <= maghigh)),
        ]
        keep = np.ones(len(sexcat), dtype=bool)
        print(f'{"cut":22s} {"alone":>7s} {"cumulative":>11s}')
        for name, m in cuts:
            keep &= m
            print(f'{name:22s} {m.sum():7d} {keep.sum():11d}')

    def size_report(sexcat, zp, x_lo, x_hi, y_lo, y_hi, maglow, maghigh, pixscale=0.498):
        mag = sexcat['MAG_AUTO'] + zp
        m = ((sexcat['FLAGS'] <= 1) &
             (sexcat['X_IMAGE'] > x_lo) & (sexcat['X_IMAGE'] < x_hi) &
             (sexcat['Y_IMAGE'] > y_lo) & (sexcat['Y_IMAGE'] < y_hi) &
             (mag >= maglow) & (mag <= maghigh))
        r = np.asarray(sexcat['FLUX_RADIUS'][m]) * pixscale
        print(f'{m.sum()} sources in crop/mag window (no size or peak cut)')
        print('FLUX_RADIUS (arcsec) percentiles 5/16/50/84/95:', np.round(np.percentile(r, [5, 16, 50, 84, 95]), 2))
        hist, edges = np.histogram(r, bins=np.arange(0, 3.01, 0.2))
        for h, e in zip(hist, edges[:-1]):
            print(f'{e:.1f}-{e + 0.2:.1f}" {h:4d} ' + '#' * (h // 2))

    # print(sexcat.colnames)
    min_size_arcsec = 0.6
    if adv:
        inner_primesources = sexcat[(sexcat['FLAGS'] == 0) &
                                    (sexcat['X_IMAGE'] < (max_x - right_crop)) & (sexcat['X_IMAGE'] > left_crop)
                                    & (sexcat['Y_IMAGE'] < (max_y) - y_crop) & (
                                                sexcat['Y_IMAGE'] > y_crop) & (sexcat['FLUX_RADIUS'] * 0.498 > min_size_arcsec)
                                    & (sexcat['FLUX_MAX'] / sexcat['FLUX_AUTO'] < 0.5)]
    else:
        inner_primesources = sexcat[(sexcat['FLAGS'] <= 1) &
                                    (sexcat['X_IMAGE'] < (max_x - right_crop)) & (sexcat['X_IMAGE'] > left_crop)
                                    & (sexcat['Y_IMAGE'] < (max_y) - y_crop) & (
                                                sexcat['Y_IMAGE'] > y_crop) & (sexcat['FLUX_RADIUS'] * 0.498 > min_size_arcsec)
                                    & (sexcat['FLUX_MAX'] / sexcat['FLUX_AUTO'] < 0.25)
                                    ]

    # zp correction
    print('applying zp correction to PRIME mags: %s' % zp)
    inner_primesources['MAG_AUTO'] = inner_primesources['MAG_AUTO'] + zp

    # inner_primesources.write('prime_all.ecsv', overwrite=True)
    inner_primesources = inner_primesources[(inner_primesources['MAG_AUTO'] >= maglow) &
                                            (inner_primesources['MAG_AUTO'] <= maghigh)]

    # cut_report(sexcat, zp, left_crop, max_x - right_crop, y_crop, max_y - y_crop, maglow, maghigh)
    # size_report(sexcat, zp, left_crop, max_x - right_crop, y_crop, max_y - y_crop, maglow, maghigh)

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


def _find_transform(cat_xy, prime_xy, n_ctrl, pixel_tol, timeout=150):
    def _handler(signum, frame):
        raise _AlignTimeout(f'exceeded {timeout} s')

    old_tol = aa.PIXEL_TOL
    aa.PIXEL_TOL = pixel_tol
    use_alarm = True
    try:
        old_handler = signal.signal(signal.SIGALRM, _handler)
    except ValueError:                 # not in the main thread: no timeout available
        use_alarm = False
    try:
        if use_alarm:
            signal.alarm(timeout)
        return aa.find_transform(cat_xy, prime_xy, max_control_points=n_ctrl)
    finally:
        if use_alarm:
            signal.alarm(0)
            signal.signal(signal.SIGALRM, old_handler)
        aa.PIXEL_TOL = old_tol


def effective_transform(M, ref_xy, use_scale=False):
    """M: 3x3 astroalign transform (cat -> prime, 1-based pixels). With use_scale=False, keep only
    the rotation plus the shift at CRPIX."""
    if use_scale:
        return M
    th = np.arctan2(M[1, 0], M[0, 0])
    A = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
    ref = np.asarray(ref_xy, float)
    Meff = np.eye(3)
    Meff[:2, :2] = A
    Meff[:2, 2] = M[:2, :2] @ ref + M[:2, 2] - A @ ref     # T(CRPIX) unchanged
    return Meff


def apply_transform_to_header(header, M):
    """Copy of header whose WCS satisfies new_world2pix = T(old_world2pix), T = M."""
    h = header.copy()
    w = WCS(header)
    ref = np.array([h['CRPIX1'], h['CRPIX2']], dtype=float)
    q0 = (np.linalg.inv(M) @ np.array([ref[0], ref[1], 1.0]))[:2]      # T^-1(CRPIX)
    h['CRVAL1'], h['CRVAL2'] = w.all_pix2world([q0], 1)[0]

    Ainv = np.linalg.inv(M[:2, :2])
    if 'CD1_1' in h:
        cd = np.array([[h['CD1_1'], h.get('CD1_2', 0.0)], [h.get('CD2_1', 0.0), h['CD2_2']]])
        (h['CD1_1'], h['CD1_2']), (h['CD2_1'], h['CD2_2']) = cd @ Ainv
    else:   # PC matrix + CDELT
        pc = np.array([[h.get('PC1_1', 1.0), h.get('PC1_2', 0.0)],
                       [h.get('PC2_1', 0.0), h.get('PC2_2', 1.0)]])
        (h['PC1_1'], h['PC1_2']), (h['PC2_1'], h['PC2_2']) = pc @ Ainv
    return h


def _thin(xy, cell_px):
    """xy sorted brightest-first; keep the brightest source in each cell_px x cell_px cell."""
    cx = np.floor((xy[:, 0] - xy[:, 0].min()) / cell_px).astype(int)
    cy = np.floor((xy[:, 1] - xy[:, 1].min()) / cell_px).astype(int)
    _, first = np.unique(cx * 100000 + cy, return_index=True)
    return xy[np.sort(first)]


def astroalign_shift(inner_primesources, inner_catsources, pixel_tol=2, max_control_points=400,
                     match_tol=1.0, max_shift=None, min_matches=10, min_sigma=5,
                     max_rot_deg=2.0, max_scale_dev=0.02, ref_xy=None,
                     prime_sort_col='MAG_AUTO', cat_sort_col=None):
    prime_xy = _xy_array(inner_primesources, prime_sort_col, is_mag=True)
    cat_xy = _xy_array(inner_catsources, cat_sort_col, is_mag=True)

    attempts = [(150, 20), (100, 30), (70, 45), (None, 180)]   # (cell size in px, timeout in s)
    result = None
    for cell, t in attempts:
        if cell is None:
            p_xy, c_xy = prime_xy, cat_xy
            n = min(max(len(prime_xy), len(cat_xy)), max_control_points)
        else:
            p_xy, c_xy = _thin(prime_xy, cell), _thin(cat_xy, cell)
            n = max(len(p_xy), len(c_xy))
        try:
            result = _find_transform(c_xy, p_xy, n, pixel_tol, timeout=t)
            print(f'astroalign succeeded (cell={cell}, {len(p_xy)}/{len(c_xy)} points)')
            break
        except (aa.MaxIterError, ValueError, TypeError, _AlignTimeout) as e:
            print(f'astroalign attempt (cell={cell}) failed: {type(e).__name__}: {e}')
    if result is None:
        return 0, 0, 0, None
    tform, (cat_matched, prime_matched) = result

    rot_deg = np.degrees(tform.rotation)
    if abs(rot_deg) > max_rot_deg or abs(tform.scale - 1) > max_scale_dev:
        print(f'astroalign solution rejected (rotation={rot_deg:.3f} deg, scale={tform.scale:.5f})')
        return 0, 0, 0, None

    # Re-match the FULL lists under the fitted transform
    M = tform.params
    cat_t = cat_xy @ M[:2, :2].T + M[:2, 2]
    dist, _ = cKDTree(prime_xy).query(cat_t, distance_upper_bound=match_tol)
    n_ok = int(np.isfinite(dist).sum())

    # Expected number of chance matches: catalog points landing inside the prime footprint,
    # times the prime source density, times the match area
    lo, hi = prime_xy.min(axis=0), prime_xy.max(axis=0)
    inside = np.all((cat_t > lo) & (cat_t < hi), axis=1)
    expected = inside.sum() * len(prime_xy) * np.pi * match_tol**2 / np.prod(hi - lo)
    sigma = (n_ok - expected) / np.sqrt(max(expected, 1.0))

    if ref_xy is not None:          # shift evaluated at the reference pixel (CRPIX)
        ref = np.asarray(ref_xy, dtype=float)
        x_shift, y_shift = (ref @ M[:2, :2].T + M[:2, 2]) - ref
    else:                           # fallback: median over astroalign's matched pairs
        d = np.asarray(prime_matched) - np.asarray(cat_matched)
        x_shift, y_shift = np.median(d[:, 0]), np.median(d[:, 1])

    print(f'astroalign: shift = ({x_shift:.3f}, {y_shift:.3f}) pix, rot = {rot_deg:.3f} deg, '
          f'scale = {tform.scale:.5f}\n'
          f'  {n_ok} full-list matches within {match_tol} pix '
          f'(chance expectation {expected:.1f}, {sigma:.0f} sigma)')

    if n_ok < min_matches or sigma < min_sigma:
        print('astroalign solution rejected: match count not significant')
        return 0, 0, n_ok, None
    if max_shift is not None and np.hypot(x_shift, y_shift) > max_shift:
        print(f'astroalign shift {np.hypot(x_shift, y_shift):.1f} pix exceeds max_shift, rejecting')
        return 0, 0, n_ok, None

    return x_shift, y_shift, n_ok, M, rot_deg


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


def check_shift(xy_shift, directory, header, data, imageName, filter_used, eff_boxsize, crop, coords,
                catNum, magcol, maglow, maghigh, errbits, catname, bulge=False, adv=False):
    """Apply the shift to a temp copy of the image, re-run the query/crossmatch, return completion fraction.
    (This is the body of the old iterate_and_test loop, for a single shift.)"""
    x_shift, y_shift = xy_shift
    header_copy = header.copy()
    w = WCS(header_copy)
    crpix1, crpix2 = header_copy['CRPIX1'], header_copy['CRPIX2']

    # subtract for CRVAL, as in the original
    header_copy['CRVAL1'], header_copy['CRVAL2'] = w.all_pix2world(
        [[crpix1 - x_shift, crpix2 - y_shift]], 1)[0]

    imageshiftname = os.path.splitext(imageName)[0] + '.shift.fits'
    fits.writeto(os.path.join(directory, imageshiftname), data, header_copy, overwrite=True)

    prime_iter, cat_iter, colnames, wcs = shiftiteration(
        directory, imageshiftname, filter_used, coords, maglow, maghigh, eff_boxsize, crop, catNum,
        magcol, errbits, x_shift=x_shift, y_shift=y_shift, catname=catname, adv=adv)

    cat_coords = SkyCoord(ra=cat_iter[colnames[0]], dec=cat_iter[colnames[1]], frame='icrs', unit='degree')
    prime_coords = utils.pixel_to_skycoord(prime_iter['X_IMAGE'], prime_iter['Y_IMAGE'], wcs, origin=1)

    thresh = (0.6 if bulge else 1.0) * u.arcsec
    idx_prime, idx_cat, _, _ = cat_coords.search_around_sky(prime_coords, thresh)

    print('PRIME Source Num = ', len(prime_iter))
    print('Catalog Source Num = ', len(cat_iter))
    print('Crossmatch Source Num = ', len(idx_prime))
    if len(prime_iter) > len(cat_iter):
        completion = len(idx_cat) / len(cat_iter)
    else:
        completion = len(idx_prime) / len(prime_iter)
    print(f'({x_shift:.5f}, {y_shift:.5f}) -> Crossmatch completion: {completion:.2f}')
    return completion


def check_shift_fast(x_shift, y_shift, header, inner_primesources, inner_catsources, colnames, bulge=False):
    h = header.copy()
    h['CRVAL1'], h['CRVAL2'] = WCS(h).all_pix2world(
        [[h['CRPIX1'] - x_shift, h['CRPIX2'] - y_shift]], 1)[0]
    w = WCS(h)

    prime = utils.pixel_to_skycoord(inner_primesources['X_IMAGE'], inner_primesources['Y_IMAGE'], w, origin=1)
    cat = SkyCoord(ra=inner_catsources[colnames[0]], dec=inner_catsources[colnames[1]],
                   frame='icrs', unit='degree')
    thresh = (0.6 if bulge else 1.0) * u.arcsec
    idx_prime, idx_cat, _, _ = cat.search_around_sky(prime, thresh)

    n_p, n_c = len(inner_primesources), len(inner_catsources)
    completion = (len(idx_cat) / n_c) if n_p > n_c else (len(idx_prime) / n_p)
    print(f'Crossmatch completion: {completion:.2f}')

    # for r in (1, 2, 3, 5):
    #     _, ic, _, _ = cat.search_around_sky(prime, r * u.arcsec)
    #     print(f'  completion within {r}": {len(np.unique(ic)) / n_c:.2f}')

    return completion


def change_all_files(xfinal_shift, yfinal_shift, directory, transform=None, rotation=None):
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

            w = WCS(header)

            crval1 = header['CRVAL1']
            crval2 = header['CRVAL2']
            crpix1 = header['CRPIX1']
            crpix2 = header['CRPIX2']

            if transform is not None:
                header = apply_transform_to_header(header, transform)
            else:
                # converting pixel shifts to world coord shifts
                new_crval1, new_crval2 = w.all_pix2world(
                    [[crpix1 - xfinal_shift, crpix2 - yfinal_shift]], 1
                )[0]

                header['CRVAL1'] = new_crval1
                header['CRVAL2'] = new_crval2

            header.set('X_SHIFT', xfinal_shift, 'X Value shift', after='WCSAXES')
            header.set('Y_SHIFT', yfinal_shift, 'Y Value shift', after='X_SHIFT')

            if rotation is not None:
                header.set('ROT_VAL', rotation, 'Corrected Rot. Value (deg)', before='X_SHIFT')

            header.set('OLD_CV1', crval1, 'Initial CRVAL1 Value', after='Y_SHIFT')
            header.set('OLD_CV2', crval2, 'Initial CRVAL2 Value', after='OLD_CV1')

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
    except OSError:
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
        thresh_high=defaults['thresh_high'], iters=defaults['iters'], single_solve=defaults['single_solve'], adv_solve=defaults['adv_solve'],
        x=defaults['x_guess'],y=defaults['y_guess']
):

    x_offset = 0
    if band == 'Y':
        print('Switching Y band to J for ease of astrometry...')
        filter_used = 'J'
    else:
        filter_used = band

    start_time = dt.now()

    # set_vizier_mirror()

    maglow = 12
    maghigh = 14

    data, header, w, raImage, decImage, zp, catname, bulge = imaging(directory, imagename, x_offset=x_offset, x_guess=x, y_guess=y)
    if bulge:
        boxsize, crop = boxchange(4)
        thresh_high = 0.25
    elif adv_solve:
        boxsize, crop = boxchange(6)
        thresh_high = 0.25
    else:
        boxsize, crop = boxchange(10, x_offset=x_offset)

    # if adv_solve:
    #     if os.path.isfile(os.path.splitext(imagename)[0] + '.psf.cat'):
    #         print('Previous psf cat detected, saving you some time and skipping all sextractor steps...')
    #         catname = os.path.splitext(imagename)[0] + '.psf.cat'
    #     else:
    #         catname = sex1(imagename, bulge=bulge)
    #         # psfex(catname)
    #         # catname = sex2(imagename)
    # else:
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
        transform = None
        rotation = None

    else:
        try:
            ref = (header['CRPIX1'], header['CRPIX2'])
            x_s, y_s, n_match, M, rot_deg = astroalign_shift(inner_primesources, inner_catsources, max_control_points=num,
                                                    max_shift=length, cat_sort_col=colnames[2], ref_xy=ref)
            transform = None
            rotation = None
            ultimate_shift_x = ultimate_shift_y = 0
            if M is not None:
                Meff = effective_transform(M, ref, use_scale=True)
                new_header = apply_transform_to_header(header, Meff)
                completion = check_shift_fast(0, 0, new_header, inner_primesources, inner_catsources,
                                              colnames, bulge=bulge)
                if completion >= thresh_low:
                    ultimate_shift_x, ultimate_shift_y, transform, rotation = x_s, y_s, Meff, rot_deg
                else:
                    print(f'WARNING: crossmatch completion {completion:.2f} < {thresh_low}, assuming 0 shift')
        except Exception as e:
            print(f'Error in final shift generation, assuming 0, error: {e}')
            ultimate_shift_x = ultimate_shift_y = 0

    if not single_solve:
        change_all_files(ultimate_shift_x, ultimate_shift_y, directory, transform=transform, rotation=rotation)
        if ultimate_shift_x != 0:
            removal(directory)

    proc_single_img_path = os.path.join(directory, imagename)
    end_time = dt.now()
    print('astrometric shift correction time:', (end_time - start_time).total_seconds())

    return ultimate_shift_x, ultimate_shift_y, proc_single_img_path


def main():
    parser = argparse.ArgumentParser(description='Corrects for translation in initial astrometry '
                                                 '(so astrom.net doesnt need to be used)')
    parser.add_argument('-single_solve', action='store_true', help='optional flag to test for a successful solve, '
                                                                   'doesnt apply to all files in directory.',
                        default=defaults['single_solve'])
    parser.add_argument('-adv_solve', action='store_true',
                        help='optional flag to run psfex and sextractor again for max'
                             ' astrometric & mag accuracy, designed for very dense '
                             'fields, but takes a long time.', default=defaults['adv_solve'])
    parser.add_argument('-dir', type=str, help='[str] path where input file is stored (should run on proc. image, '
                                               'so likely should be /C#_sub/)')
    parser.add_argument('-imagename', type=str, help='[str] input file name (should run on proc. image, '
                                                     'i.e. *.sky.flat.fits)')
    parser.add_argument('-band', type=str, help='[str] filter used, ex. "J"')
    parser.add_argument('-x', type=float, help='[float] optional initial guess for X, '
                                               '*PROVIDE BOTH AN X AND Y GUESS USING -x & -y*',
                        default=defaults['x_guess'])
    parser.add_argument('-y', type=float, help='[float] optional initial guess for Y, '
                                               '*PROVIDE BOTH AN X AND Y GUESS USING -x & -y*',
                        default=defaults['y_guess'])
    parser.add_argument('-length', type=float, help='[float] optional, value over which dists will not be '
                                                    'considered (pix), default = 100', default=defaults['length'])
    parser.add_argument('-num', type=int, help='[int] optional, # of sources, sorted by mag, to consider in '
                                               'dist. calculation, default = 400', default=defaults['num'])
    parser.add_argument('-thresh_low', type=float, help='[float] Failure threshold for crossmatch betw. '
                                                        'PRIME & astroquery catalog, default = 0.1',
                        default=defaults['thresh_low'])
    parser.add_argument('-thresh_high', type=float, help='[float] Successful threshold for crossmatch betw. '
                                                         'PRIME & astroquery catalog, default = 0.4',
                        default=defaults['thresh_high'])
    parser.add_argument('-iters', type=int, help='[int] Number of iterations to test shift values before'
                                                 ' defaulting to highest crossmatch percentage, default = 10',
                        default=defaults['iters'])
    args, unknown = parser.parse_known_args()

    shift(args.dir, args.imagename, args.band, args.length, args.num, args.thresh_low, args.thresh_high, args.iters,
          args.single_solve, args.adv_solve, args.x, args.y)


if __name__ == "__main__":
    main()