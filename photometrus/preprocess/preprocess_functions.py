"""
List of preprocess functions to use in main pipeline
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


def flatfielding(imagepaths, FF_imagepaths, band, chip, date=None, **kwargs):
    """
    Runs flatfielding for normal pipeline runs

    Parameters
    ----------
    imagepaths: list
        List of input images paths in str form
    FF_imagepaths: str
        Directory to storage of flatfielded images (ex. '/J/C1_FF')
    band: str
        Band of observation
    chip: int
        Number of detector
    date: str
        Date of observation in YYYYMMDD format (ex. '20260606')

    Returns
    -------
    flatpath: str
        Filepath to utilized flat field
    ff_img_paths
        Path or paths to flatfielded images
    """

    os.chdir(gen_pipeline_file_name())
    print('using master flat to flat field ramp image..')
    flatpath = gen_mflat_file_name(band, chip, date, sflat=True)
    # flatpath = gen_sflat_file_name(band, chip, date)

    print('\nEquivalent argparse cmd: photometrus process flatfield -in_path %s -out_path %s'
          ' -flat_path %s' % (imagepaths, FF_imagepaths, flatpath))

    ff_img_paths = flatfield.flat_field_cmd(in_path=imagepaths, out_path=FF_imagepaths, flat_path=flatpath)

    return ff_img_paths, {'flatpath': flatpath}


def _flatfielding_step(imagepaths, **kwargs):
    """Wrapper: runs flatfielding but returns only the new imagepaths."""
    flatpath, ff_img_paths = flatfielding(imagepaths, **kwargs)
    return ff_img_paths


def gen_bad_pix_mask_by_image(imagepaths, astrom_dir, chip, **kwargs):
    """
    Generates final bad pixel mask per image using bad pix mask (per chip).  Bad pixels will be turned into nan values
    on individual science images, while the binary final bad pix mask will be written to astrom_dir to be used with swarp

    Parameters
    ----------
    imagepaths: list or str
        List of input images paths in str form, or single image path
    astrom_dir: str
        Directory to full processed image directory (ex. './J/C1_astrom/'), final bad pix masks are placed here for
        swarp to use
    chip: int
        Number of detector

    Returns
    -------
    imagepaths
        Path or paths to nan masked images
    """

    badpixmap = fits.getdata(gen_mask_file_name('badpixmask_new_c%i.fits' % chip))
    # satulimMap = fits.getdata(gen_mask_file_name('satulim.C%i.fits' % chip))[4:4092, 4:4092]

    if type(imagepaths) is list:
        print(f'\nRunning bad pix & saturation map masking on FF images')
        for imagepath in imagepaths:
            image = fits.getdata(imagepath)
            header = fits.getheader(imagepath)

            image[badpixmap == 0] = np.nan
            # satmasked = image >= satulimMap
            # image[satmasked] = np.nan

            image_bad_map = (~np.isnan(image)).astype('uint16')
            imagename = os.path.basename(imagepath)
            base_imagename = imagename.split('.')[0]
            maskpath = os.path.join(astrom_dir, f"{base_imagename}.weightmap.fits")

            fits.HDUList(fits.PrimaryHDU(header=header, data=image)).writeto(imagepath, overwrite=True)
            fits.HDUList(fits.PrimaryHDU(header=header, data=image_bad_map)).writeto(maskpath, overwrite=True)

    elif os.path.isfile(imagepaths):
        image = fits.getdata(imagepaths)
        header = fits.getheader(imagepaths)

        image[badpixmap == 0] = np.nan
        # satmasked = image >= satulimMap
        # image[satmasked] = np.nan

        image_bad_map = (~np.isnan(image)).astype('uint8')
        imagename = os.path.basename(imagepaths)
        base_imagename = imagename.split('.')[0]
        maskpath = os.path.join(astrom_dir, f"{base_imagename}.weightmap.fits")

        fits.HDUList(fits.PrimaryHDU(header=header, data=image)).writeto(imagepaths, overwrite=True)
        fits.HDUList(fits.PrimaryHDU(header=header, data=image_bad_map)).writeto(maskpath, overwrite=True)

    else:
        raise TypeError('imagepaths type is not acceptable!  Must be list of input filepaths or str of single filepath')

    print(f'All final bad pix masks deposited at: {astrom_dir}')
    print('Nan masking of science images & final mask creations completed!\n')
    return imagepaths


def interpolate_bad_pix(imagepaths, chip, kernel_stdev=2, **kwargs):
    """
    Using bad pix mask (already generated per chip) & saturation map (per chip), bad pixels are converted to nans in
    individual images.  A gaussian interpolation is then run on each image to interpolate the nan values

    Parameters
    ----------
    imagepaths: list or str
        List of input images paths in str form, or single image path
    chip: int
        Number of detector
    kernel_stdev: int
        Size of gaussian kernels in stdevs

    Returns
    -------
    imagepaths
        Path or paths to nan masked images
    """

    badpixmap = fits.getdata(gen_mask_file_name('badpixmask_new_c%i.fits' % chip))
    satulimMap = fits.getdata(gen_mask_file_name('satulim.C%i.fits' % chip))[4:4092, 4:4092]

    if type(imagepaths) is list:
        print(f' Running masking and interpolation on FF images')
        for imagepath in imagepaths:
            image = fits.getdata(imagepath)
            header = fits.getheader(imagepath)

            image[badpixmap == 0] = np.nan
            satmasked = image >= satulimMap
            image[satmasked] = np.nan

            kernel = Gaussian2DKernel(x_stddev=kernel_stdev)
            nan_percent = 100 * np.count_nonzero(np.isnan(image)) / (image.shape[0] * image.shape[1])
            print(' start_nan_percentage: {}'.format(nan_percent))
            print(' Running kernel interpolation...')
            filled = interpolate_replace_nans(image, kernel)
            remaining_nans = np.isnan(filled).sum()
            end_nan_percent = 100 * remaining_nans / (image.shape[0] * image.shape[1])
            print(' end_nan_percentage: {}'.format(end_nan_percent))

            fits.HDUList(fits.PrimaryHDU(header=header, data=image)).writeto(imagepath, overwrite=True)

    elif os.path.isfile(imagepaths):
        image = fits.getdata(imagepaths)
        header = fits.getheader(imagepaths)

        image[badpixmap == 0] = np.nan
        satmasked = image >= satulimMap
        image[satmasked] = np.nan

        kernel = Gaussian2DKernel(x_stddev=kernel_stdev)
        nan_percent = 100 * np.count_nonzero(np.isnan(image)) / (image.shape[0] * image.shape[1])
        print(' start_nan_percentage: {}'.format(nan_percent))
        print(' Running kernel interpolation...')
        filled = interpolate_replace_nans(image, kernel)
        remaining_nans = np.isnan(filled).sum()
        end_nan_percent = 100 * remaining_nans / (image.shape[0] * image.shape[1])
        print(' end_nan_percentage: {}'.format(end_nan_percent))

        fits.HDUList(fits.PrimaryHDU(header=header, data=image)).writeto(imagepaths, overwrite=True)

    return imagepaths


PREPROCESS_FCTNS = {
    'flatfielding': flatfielding,
    'interpolate_bad_pix': interpolate_bad_pix,
    'gen_bad_pix_mask_by_image': gen_bad_pix_mask_by_image,
}


def preprocess_fctns(imagepaths, commands=None, **kwargs):
    """
    Runs all preprocessing function steps defined in `commands` (defaults to default_preproc).
    kwargs is a superset of everything any step might need, EACH FUNCTION SHOULD HAVE THE SAME OUTPUT.
    Each step function only pulls out what named vars it needs.
    Anything else falls into its own **kwargs and is ignored.
    """
    if commands is None:
        commands = defaults['default_preproc']
    for cmd in commands:
        result = PREPROCESS_FCTNS[cmd](imagepaths, **kwargs)
        if isinstance(result, tuple):
            imagepaths, extra = result
            kwargs.update(extra)
        else:
            imagepaths = result
    return imagepaths
