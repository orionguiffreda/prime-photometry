"""
List of stacking functions to use in main pipeline
"""
import os
import sys
import subprocess
import argparse
import fnmatch
import shutil
from astropy.io import fits
import numpy as np
from astropy.convolution import interpolate_replace_nans, Gaussian2DKernel
from datetime import datetime as dt

from photometrus.settings import gen_pipeline_file_name, gen_mflat_file_name, gen_sflat_file_name
from photometrus.preprocess import astromangle_new
from photometrus.preprocess import astromangle_wcs
from photometrus.preprocess import gen_astrometry
from photometrus.preprocess import flatfield
from photometrus.sky import gen_sky
from photometrus.sky import sky_sub
from photometrus.astrom import astrom_shift
from photometrus.astrom import astrom_shift_new
from photometrus.astrom import astrometry
from photometrus.astrom.astrometry import sex, local_scamp
from photometrus.stack import stack
from photometrus.stack.proc_img_filtering import (shallow_filter, bad_astrom_filter, bad_focus_filter, bckg_filter,
                                                  prune_offending_images)

from photometrus.settings import (bulge_checker, auto_bulge_detect, gen_mask_file_name)
from photometrus.utils.defaults import PROCESSING_DEFAULTS as defaults


# image quality filtering for individual images


def img_qual_filtering(astromdir, band, stdev_limit=None, **kwargs):
    if stdev_limit is None:
        stdev_limit = [3, 3, 4]
        shallow_lim = stdev_limit[0]
        bad_astr_lim = stdev_limit[1]
        bad_foc_lim = stdev_limit[2]
    else:
        shallow_lim = bad_astr_lim = bad_foc_lim = stdev_limit
    print('\nRunning quality filters over all processed images!')

    img_list, img_cat_list = sex(imgdir=astromdir, param_file='astrom_world.param')

    sky_fac_arr, sky_fac_thresh, csv_path = bckg_filter(full_img_list=img_list,
                                                        stdev_limit=shallow_lim)

    source_num_arr, source_num_thresh, csv_path = shallow_filter(full_img_list=img_list,
                                                                 img_cat_list=img_cat_list,
                                                                 stdev_limit=shallow_lim,
                                                                 csv_path=csv_path)

    completion_arr, compl_thresh, csv_path = bad_astrom_filter(full_img_list=img_list,
                                                               img_cat_list=img_cat_list,
                                                               band=band,
                                                               stdev_limit=bad_astr_lim,
                                                               csv_path=csv_path)

    flx_rad_arr, flx_rad_thresh, csv_path = bad_focus_filter(full_img_list=img_list,
                                                             img_cat_list=img_cat_list,
                                                             stdev_limit=bad_foc_lim,
                                                             csv_path=csv_path)

    good_img_list, good_img_cat_list = prune_offending_images(img_list, img_cat_list, sky_fac_arr, sky_fac_thresh,
                                                              source_num_arr,
                                                              source_num_thresh, completion_arr, compl_thresh,
                                                              flx_rad_arr, flx_rad_thresh,
                                                              )

    print('Filtering complete!\n')

    return good_img_list


# stacking


def stacking(astromdir, stackdir, chip, **kwargs):
    os.chdir(gen_pipeline_file_name())
    print('stacking all images using SWARP...')

    print(
        '\nEquivalent argparse cmd: photometrus process stack -sub %s -stack %s -chip %i' % (astromdir, stackdir, chip))

    stack.stack(subpath=astromdir, stackpath=stackdir, chip=chip)

    return None


STACK_FCTNS = {
    'img_qual_filtering': {'func': img_qual_filtering, 'output_key': 'astromdir'},
    'stacking': {'func': stacking, 'output_key': None},
}


def stack_fctns(commands=None, **kwargs):
    """
    Runs all stacking function steps defined in `commands` (defaults to default_stack_proc).
    kwargs is a superset of everything any step might need, EACH FUNCTION SHOULD HAVE THE SAME OUTPUT.
    Each step function only pulls out what named vars it needs.
    Anything else falls into its own **kwargs and is ignored.
    """
    if commands is None:
        commands = defaults['default_stack_proc']
    for cmd in commands:
        entry = STACK_FCTNS[cmd]
        func, output_key = entry['func'], entry['output_key']
        result = func(**kwargs)
        if output_key is not None and result is not None:
            kwargs[output_key] = result

    return kwargs
