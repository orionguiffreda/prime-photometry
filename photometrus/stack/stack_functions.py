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
from photometrus.stack import stack

from photometrus.settings import (bulge_checker, auto_bulge_detect, gen_mask_file_name)
from photometrus.utils.defaults import PROCESSING_DEFAULTS as defaults

# stacking


def stacking(astromdir, stackdir, chip):
    os.chdir(gen_pipeline_file_name())
    print('stacking all images using SWARP...')

    print('\nEquivalent argparse cmd: photometrus process stack -sub %s -stack %s -chip %i' % (astromdir, stackdir, chip))

    stack.stack(subpath=astromdir, stackpath=stackdir, chip=chip)


STACK_FCTNS = {
    'stacking': stacking,
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
        STACK_FCTNS[cmd](**kwargs)