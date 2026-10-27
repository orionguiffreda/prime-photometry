"""Filtering functions designed to screen out bad images / observations before stacking"""

import os
import shutil
import argparse
import subprocess
import sys
from astropy.io import fits
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.wcs import WCS
import numpy as np
import fnmatch
import matplotlib.pyplot as plt
import pandas as pd
from astropy.stats import SigmaClip

from photometrus.photometry.photom_utils import get_table_from_ldac
from photometrus.query import query, gaia_query
from photometrus.stack.stack import get_astrom_files

from photometrus.utils.defaults import PROCESSING_DEFAULTS as defaults


def generate_img_qual_csv(img_list, data_vals_arr, colname, csv_path=None):
    """
    Generates csv for all images and specific data value
    """

    # Add to or create CSV of values for each file
    new_df = pd.DataFrame({
        'img_filepath': img_list,
        colname: data_vals_arr,
    })

    csv_name = f'IMG_QUAL_STATS_{os.path.basename(img_list[0])}-{os.path.basename(img_list[-1])}.csv'
    if not csv_path:
        dirname = os.path.dirname(img_list[0]).rstrip('/')
        parent_dirname = os.path.dirname(dirname)
        csv_path = os.path.join(parent_dirname, csv_name)

    if os.path.exists(csv_path):
        existing_df = pd.read_csv(csv_path)

        if 'img_filepath' not in existing_df.columns:
            raise ValueError(f"Existing CSV at {csv_path} has no 'img_filepath' column to merge on.")

        # Drop old column if present so it gets cleanly overwritten
        if colname in existing_df.columns:
            existing_df = existing_df.drop(columns=[colname])

        merged_df = existing_df.merge(new_df, on='img_filepath', how='outer')
    else:
        merged_df = new_df

    merged_df.to_csv(csv_path, index=False)

    return csv_path


def bckg_filter(full_img_list, stdev_limit=3, csv_path=None):
    """
    Screens out images with abnormal sky scaling factors (SKY_FAC) by looking at SKY_FAC and comparing to median,
    imgs w/ total source #s below specified sig of the median are omitted from stack
    """

    print('\nRunning filter for images w/ bad background..')

    sky_fac_arr = []
    for img in full_img_list:
        hdr = fits.getheader(img)

        try:
            sky_fac = hdr['SKY_FAC']
        except KeyError:
            print('SKY_FAC header val not found in image! Defaulting to 1')
            sky_fac = 1

        print(f'     Sky Scaling Factor: {sky_fac:.1f}')
        sky_fac = round(sky_fac, 1)
        sky_fac_arr.append(sky_fac)

    med_sky_fac = np.median(sky_fac_arr, axis=0)
    stdev_sky_fac = np.std(sky_fac_arr)
    low_thresh = med_sky_fac - (stdev_limit * stdev_sky_fac)
    high_thresh = med_sky_fac + (stdev_limit * stdev_sky_fac)
    thresh = [low_thresh, high_thresh]
    print(f'     Median: {med_sky_fac:.1f}, stdev: {stdev_sky_fac:.1f}, thresh ({stdev_limit} sig): '
          f'{low_thresh:.1f} - {high_thresh:.1f}')

    # Add to or create CSV of values for each file
    csv_path = generate_img_qual_csv(img_list=full_img_list, data_vals_arr=sky_fac_arr, colname='SKY_FAC',
                                     csv_path=csv_path)

    return sky_fac_arr, thresh, csv_path


def shallow_filter(full_img_list, img_cat_list, crop=defaults['crop'], stdev_limit=3, csv_path=None):
    """
    Screens out shallow images by looking at total # of sources and comparing to median, imgs w/ total source #s
    below specified sig of the median are omitted from stack
    """

    print('\nRunning filter for shallow images..')
    img_data = fits.getdata(full_img_list[0])
    max_x = img_data.shape[0]
    max_y = img_data.shape[1]

    source_num_arr = []
    for img,cat in zip(full_img_list, img_cat_list):
        tbl = get_table_from_ldac(cat)
        good_tbl = tbl[(tbl['XWIN_IMAGE'] < (max_x - crop)) & (tbl['XWIN_IMAGE'] > crop) &
                             (tbl['YWIN_IMAGE'] < (max_y) - crop) & (tbl['YWIN_IMAGE'] > crop)]

        indiv_source_num = len(good_tbl)
        print(f'     Source #: {indiv_source_num}')
        source_num_arr.append(indiv_source_num)

    med_source_num = np.median(source_num_arr, axis=0)
    stdev_source_num = np.std(source_num_arr)
    low_thresh = med_source_num - (stdev_limit * stdev_source_num)
    print(f'     Median: {med_source_num:.0f}, stdev: {stdev_source_num:.0f}, thresh ({stdev_limit} sig): {low_thresh:.0f}')

    # Add to or create CSV of values for each file
    csv_path = generate_img_qual_csv(img_list=full_img_list, data_vals_arr=source_num_arr, colname='Source #',
                                     csv_path=csv_path)

    return source_num_arr, low_thresh, csv_path


def crossmatch_img_to_survey(input_surv_tbl, input_img_tbl, input_wcs, input_data, crop):
    """
    Crossmatches input image catalog to input survey catalog to determine completion
    """

    max_x = input_data.shape[0]
    max_y = input_data.shape[1]

    tbl = input_img_tbl
    Q = input_surv_tbl

    try:
        # print(f'    Querying {survey} and crossmatching to determine catalog completion..')
        # Q, _, _ = query(raImage=input_ra, decImage=input_dec, band=band, w=input_wcs, data=input_data,
        #                 survey=survey)

        # getting source coords and cropping corners for both input img catalog and 2mass query

        good_img_stars = tbl[(tbl['XWIN_IMAGE'] < (max_x - crop)) & (tbl['XWIN_IMAGE'] > crop) &
                             (tbl['YWIN_IMAGE'] < (max_y) - crop) & (tbl['YWIN_IMAGE'] > crop)]

        img_RA = 'ALPHA_J2000'
        img_DEC = 'DELTA_J2000'

        query_colnames = Q[0].colnames
        Q_RA = query_colnames[0]
        Q_DEC = query_colnames[1]

        Q_imCoords = input_wcs.all_world2pix(Q[0][Q_RA], Q[0][Q_DEC], 1)

        good_Q_stars = Q[0][
            np.where((Q_imCoords[0] > crop) & (Q_imCoords[0] < (max_x - crop)) & (Q_imCoords[1] > crop) & (
                    Q_imCoords[1] < (max_y - crop)))]


        ImgCatCoords = SkyCoord(ra=good_img_stars[img_RA], dec=good_img_stars[img_DEC], frame='icrs', unit='degree')
        QueryCatCoords = SkyCoord(ra=good_Q_stars[Q_RA], dec=good_Q_stars[Q_DEC], frame='icrs', unit='degree')

        # print(' Img cropped source total = ', len(good_img_stars))
        # print(f' Chosen survey cropped source total = ', len(good_Q_stars))

        crsmtch_thresh = 1.0
        idx_img, idx_query, d2d, d3d = QueryCatCoords.search_around_sky(ImgCatCoords, crsmtch_thresh * u.arcsec)

        df = pd.DataFrame({
            'idx_img': idx_img,
            'idx_query': idx_query,
            'd2d': d2d.to(u.arcsec).value  # example in arcsec
        })

        # Sort by separation and drop duplicates, keeping the closest
        query_matches_closest = df.sort_values('d2d').drop_duplicates('idx_query', keep='first')

        # print(f'    Crossmatched image source num = {len(img_matches_closest)}')
        completion = len(query_matches_closest) / len(good_Q_stars)
        print('     Completion = %.2f' % completion)
    except Exception as e:
        print(f'    Error in Vizier query & Survey Crossmatch: {e}')
        print('     Assuming bad completion!')
        completion = 0

    return completion


def bad_astrom_filter(full_img_list, img_cat_list, band, survey='GAIA', crop=defaults['crop'], stdev_limit=3,
                      csv_path=None):
    """
    Screens out images w/ bad astrometry by querying 2mass in area & cross-matching, imgs w/ total source #s below
    specified sig of the median are omitted from stack
    """

    completion_arr = []

    print('\nRunning filter for images w/ bad astrometry or completion..')

    image_hdrs = [fits.getheader(img) for img in full_img_list]

    image_ras = [(float(hdr['CRVAL1'])) for hdr in image_hdrs]
    image_decs = [(float(hdr['CRVAL2'])) for hdr in image_hdrs]
    med_ra = np.nanmedian(image_ras)
    med_dec = np.nanmedian(image_decs)

    print(f' Querying {survey} and crossmatching to determine catalog completion..')
    # Q, _, _ = query(raImage=med_ra, decImage=med_dec, band=band, survey=survey)
    coords = SkyCoord(ra=med_ra * u.degree, dec=med_dec * u.degree, frame='fk5')
    Q, _ = gaia_query(coords=coords, width=48, chosen_frame='fk5')

    for img, cat in zip(full_img_list, img_cat_list):
        input_img = fits.open(img)
        input_hdr = input_img[0].header
        input_data = input_img[0].data
        input_wcs = WCS(input_hdr)

        tbl = get_table_from_ldac(cat)

        completion = crossmatch_img_to_survey(input_surv_tbl=Q, input_img_tbl=tbl, input_wcs=input_wcs,
                                              input_data=input_data, crop=crop)

        completion_arr.append(completion)

    completion_arr = np.array(completion_arr)
    print('     Removing broken astrometric solns < 0.05 completion from med & stdev calc!')
    valid_mask = completion_arr > 0.05

    # Uses floor above to resist outliers of obviously broken astrometric solutions
    med_completion = np.median(completion_arr[valid_mask])
    stdev_completion = np.std(completion_arr[valid_mask])
    low_thresh = med_completion - (stdev_limit * stdev_completion)
    print(f'     Median: {med_completion:.2f}, stdev: {stdev_completion:.2f}, thresh ({stdev_limit} sig): {low_thresh:.2f}')

    # Add to or create CSV of values for each file
    completion_arr = [round(val, 3) for val in completion_arr]
    csv_path = generate_img_qual_csv(img_list=full_img_list, data_vals_arr=completion_arr, colname='Completion %',
                                     csv_path=csv_path)

    # NLC failure condition
    acc_2m_compl_lim = 0.1
    if med_completion <= acc_2m_compl_lim:
        raise Exception(f'*WARNING* Median completion below {acc_2m_compl_lim}!: {med_completion:.1f}\n'
                        f'Arresting stacking due to unacceptable image quality!\n'
                        f'Image quality likely poor across obs.! Is this an error with NLC / is it cloudy? '
                        f'Check individual images!')

    return completion_arr, low_thresh, csv_path


def bad_focus_filter(full_img_list, img_cat_list, crop=defaults['crop'], stdev_limit=4, csv_path=None):
    """
    Screens out images that are out of focus by looking at flux radius and comparing to median, imgs w/ total source #s
    outside specified sig of the median are omitted from stack
    """

    print('\nRunning filter for badly focused images..')
    img_data = fits.getdata(full_img_list[0])
    max_x = img_data.shape[0]
    max_y = img_data.shape[1]

    flx_rad_arr = []
    for img,cat in zip(full_img_list, img_cat_list):
        tbl = get_table_from_ldac(cat)
        good_tbl = tbl[(tbl['XWIN_IMAGE'] < (max_x - crop)) & (tbl['XWIN_IMAGE'] > crop) &
                        (tbl['YWIN_IMAGE'] < (max_y) - crop) & (tbl['YWIN_IMAGE'] > crop)
                       & (tbl['FLUX_RADIUS'] >= 1 / 0.498)]

        indiv_flx_rad = np.nanmedian(good_tbl['FLUX_RADIUS'])
        print(f'    Flux rad: {indiv_flx_rad:.3f}')
        flx_rad_arr.append(indiv_flx_rad)

    med_flx_rad = np.median(flx_rad_arr, axis=0)
    stdev_flx_rad = np.std(flx_rad_arr)
    # low_thresh = med_flx_rad - (stdev_limit * stdev_flx_rad)
    high_thresh = med_flx_rad + (stdev_limit * stdev_flx_rad)
    print(f'    Median: {med_flx_rad:.3f}, stdev: {stdev_flx_rad:.3f}, thresh ({stdev_limit} sig): {high_thresh:.3f}')

    # Add to or create CSV of values for each file
    flx_rad_arr = [round(val, 3) for val in flx_rad_arr]
    csv_path = generate_img_qual_csv(img_list=full_img_list, data_vals_arr=flx_rad_arr, colname='Flux Radius',
                                     csv_path=csv_path)

    return flx_rad_arr, high_thresh, csv_path


def prune_offending_images(full_img_list, full_img_cat_list, sky_fac_arr, sky_fac_thresh, source_num_arr,
                           source_num_thresh, completion_arr, compl_thresh, flx_rad_arr, flx_rad_thresh):

    n = len(full_img_list)
    arrays = {'sky_fac': sky_fac_arr, 'source_num': source_num_arr,
              'completion': completion_arr, 'flx_rad': flx_rad_arr}
    for name, arr in arrays.items():
        if len(arr) != n or len(full_img_cat_list) != n:
            raise ValueError(f'Length mismatch: {name} has {len(arr)} entries, '
                             f'expected {n}')

    sky = np.asarray(sky_fac_arr)
    sky_lo, sky_hi = sky_fac_thresh

    # (label, mask of images that PASS), all computed on original indices
    filters = [
        ('background',  (sky >= sky_lo) & (sky <= sky_hi)),
        ('source num',  np.asarray(source_num_arr) >= source_num_thresh),
        ('completion',  np.asarray(completion_arr) >= compl_thresh),
        ('psf size',    np.asarray(flx_rad_arr) <= flx_rad_thresh),
    ]

    keep = np.ones(n, dtype=bool)
    for label, passed in filters:
        print(f'Pruning for {label}...')
        newly_removed = np.flatnonzero(keep & ~passed)
        if newly_removed.size:
            pruned = [full_img_list[i] for i in newly_removed]
            print(f' {newly_removed.size} filtered out! Pruned images: {pruned}')
        keep &= passed

    good = np.flatnonzero(keep)
    good_img_list = [full_img_list[i] for i in good]
    good_img_cat_list = [full_img_cat_list[i] for i in good]
    return good_img_list, good_img_cat_list
