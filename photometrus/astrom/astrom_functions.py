"""
List of astrometry functions to use in main pipeline
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


# %% initial astrometry


def astrom_angle(astromdir, parentdir, chip, rot_val=48, **kwargs):
    os.chdir(gen_pipeline_file_name())
    print('running initial astrometry on input imgs...')
    inputpath = parentdir

    outputpath = astromdir
    print(
        '\nEquivalent argparse cmd: photometrus astrom angle -input %s -output %s' % (
            inputpath, outputpath))

    # astromangle_new.astrom_angle(input_field=ramppath, output_dir=astrompath, rot_val=rot_val)
    if os.path.isfile(inputpath):
        init_astrm_path = astromangle_wcs.astrom_angle(input_field=inputpath, output_dir=outputpath, rot_val=rot_val)
        return init_astrm_path

    astromangle_wcs.astrom_angle(input_field=inputpath, output_dir=outputpath, rot_val=rot_val)


def shift(subdir, band, adv=False, old=False, single_solve=False, **kwargs):
    os.chdir(gen_pipeline_file_name())
    print('Shifting astrometry...')

    if os.path.isfile(subdir):
        all_fits = [os.path.split(subdir)[1]]
        appl_file_range = all_fits
        subpath = os.path.split(subdir)[0]
        err_msg = 'Individual fits file did not successfully solve w/ astrom_shift!  Assuming no shift!'
    else:
        all_fits = [f for f in sorted(os.listdir(subdir)) if f.endswith('.flat.fits') or f.endswith('.flat.new')]
        appl_file_range = all_fits[:len(all_fits) // 2]
        err_msg = ('Examined half of applicable fits files in directory for astrom_shift w/ no success!  '
                   'Assuming no shift!')

    if old:
        imgname = all_fits[0]
        astrom_shift.shift(directory=subdir, imagename=imgname, band=band)
    else:
        shift_x = shift_y = None
        good_fits = None
        proc_single_img_path = None
        for imgname in appl_file_range:
            print('\nEquivalent argparse cmd: photometrus astrom shift -dir %s -imagename %s -band %s' %
                  (subdir, imgname, band))
            shift_x, shift_y, proc_single_img_path = astrom_shift_new.shift(directory=subdir, imagename=imgname, band=band, adv_solve=adv,
                                                      single_solve=single_solve)
            if (shift_x != 0) or (shift_y != 0):
                good_fits = fits
                break

            print('\nEquivalent argparse cmd: photometrus astrom shift -dir %s -imagename %s -band %s -adv' %
                  (subdir, imgname, band))
            shift_x, shift_y, proc_single_img_path = astrom_shift_new.shift(directory=subdir, imagename=imgname, band=band, adv_solve=True,
                                                      single_solve=single_solve)
            if (shift_x != 0) or (shift_y != 0):
                good_fits = fits
                break

        if good_fits is None:
            print(err_msg)

        return shift_x, shift_y, proc_single_img_path


# %% astrometric verification


def verify_astrom(astromdir, subdir, chip, band, rot_val, bulge=False, **kwargs):
    initial_rot_val = rot_val
    chosen_rot_val = initial_rot_val
    attempted_angles = set()

    while True:
        # Initial shift attempt
        astrom_angle(astromdir, subdir, chip, chosen_rot_val)

        shift_x, shift_y, proc_single_img_path = shift(astromdir, band)
        if shift_x == 0 and shift_y == 0:
            print("Shift algorithm (initial) could not solve")
            print("Skipping to ROTOFF verification...")
        else:
            break

        # check if shift succeeded
        if os.path.isfile(astromdir):
            shift_fail_check = os.path.exists(os.path.splitext(astromdir)[0] + '.shift.fits')
            all_fits = ['']
        else:
            shift_fail_check = [f for f in os.listdir(astromdir) if f.endswith('.shift.fits')]
            all_fits = [f for f in sorted(os.listdir(astromdir)) if f.endswith('.flat.fits') or f.endswith('.flat.new')]
        if not shift_fail_check:
            break

        print('\nShift astrometry failed! Verifying ROTOFF val in an image FITS header...')
        if not all_fits:
            print('No suitable FITS files found in dir!')
            break

        imgpath = os.path.join(astromdir, all_fits[0])

        # Check if rotoff is real
        rotoff_real = False
        try:
            img = fits.open(imgpath)
            imghdr = img[0].header
            rotoff_check = int(imghdr['ROTOFF'])
            rotoff_real = True
            print('ROTOFF value is real: %i...' %
                  rotoff_check)
        except (ValueError, KeyError):
            print('ROTOFF value is not real! Varying ROTOFF value by +90 deg...')

        if rotoff_real:
            raise Exception('\n*PROCESSING ARRESTED, CHECK FIELD*'
                            'No suitable shift values found across half the total images! Arresting processing to prevent'
                            'bad image generation / hanging, are the images bad quality or very dense?')

        # if rotoff_real:
        #     # Advanced shift
        #     try:
        #         shift(astromdir, band, adv=True)
        #     except (IndexError, ValueError) as e:
        #         print("Shift algorithm (advanced) encountered an error: %s" % e)
        #
        #     shift_fail_check = [f for f in os.listdir(astromdir) if f.endswith('.shift.fits')]
        #     if not shift_fail_check:
        #         break
        #
        #     # old algo fallback
        #     if not bulge:
        #         print('Improved shift algorithm failed... Attempting old shift algorithm. *MAY HAVE INACCURACY*')
        #         try:
        #             shift(astromdir, band, old=True)
        #         except (IndexError, ValueError) as e:
        #             print("Shift algorithm (old) encountered an error: %s" % e)
        #
        #     # No rotation variation if rotoff was real, even if these failed
        #     break

        else:
            # Rotation variation block for if rotoff is bad
            attempted_angles.add(chosen_rot_val)
            chosen_rot_val = (chosen_rot_val + 90) % 360
            print('New ROTOFF value: %i' % chosen_rot_val)

            if chosen_rot_val in attempted_angles:
                print("All rotations failed. Moving on, but astrometry is likely to fail, so examine images further!")
                break

        return shift_x, shift_y


def verify_astrom_indiv(inputimgpath, outputpath, chip, band, rot_val=defaults['rot_val'], input_shift_x=None,
                        input_shift_y=None, bulge=False, **kwargs):
    initial_rot_val = rot_val
    chosen_rot_val = initial_rot_val
    attempted_angles = set()

    # if shift already found previously
    if input_shift_x:
        init_astrm_path = astrom_angle(outputpath, inputimgpath, chip, chosen_rot_val)

        print(f'\nApplying known shift: {input_shift_x:.2f}, {input_shift_y:.2f} to {init_astrm_path}')

        with fits.open(init_astrm_path, mode='update') as hdul:
            header = hdul[0].header
            crpix1 = header['CRPIX1']
            crpix2 = header['CRPIX2']
            header['CRPIX1'] = crpix1 + input_shift_x
            header['CRPIX2'] = crpix2 + input_shift_y
            hdul.close()

        proc_single_img_path = init_astrm_path
        return input_shift_x, input_shift_y, proc_single_img_path

    while True:
        # Initial shift attempt
        init_astrm_path = astrom_angle(outputpath, inputimgpath, chip, chosen_rot_val)

        shift_x, shift_y, proc_single_img_path = shift(init_astrm_path, band, single_solve=True)
        if shift_x == 0 and shift_y == 0:
            print("Shift algorithm (initial) could not solve")
            print("Skipping to ROTOFF verification...")
        else:
            print()
            break

        # check if shift succeeded
        shift_fail_check = os.path.exists(os.path.splitext(init_astrm_path)[0] + '.shift.fits')
        all_fits = ['']
        if not shift_fail_check:
            break

        print('\nShift astrometry failed! Verifying ROTOFF val in an image FITS header...')

        # Check if rotoff is real
        rotoff_real = False
        try:
            img = fits.open(init_astrm_path)
            imghdr = img[0].header
            rotoff_check = int(imghdr['ROTOFF'])
            rotoff_real = True
            print('ROTOFF value is real: %i...' %
                  rotoff_check)
        except (ValueError, KeyError):
            print('ROTOFF value is not real! Varying ROTOFF value by +90 deg...')

        if rotoff_real:
            raise Exception('\n*PROCESSING ARRESTED, CHECK FIELD*'
                            'No suitable shift values found across half the total images! Arresting processing to prevent'
                            'bad image generation / hanging, are the images bad quality or very dense?')

        else:
            # Rotation variation block for if rotoff is bad
            attempted_angles.add(chosen_rot_val)
            chosen_rot_val = (chosen_rot_val + 90) % 360
            print('New ROTOFF value: %i' % chosen_rot_val)

            if chosen_rot_val in attempted_angles:
                print("All rotations failed. Moving on, but astrometry is likely to fail, so examine images further!")
                break

    return shift_x, shift_y, proc_single_img_path


def astromatic_astrometry(astromdir, band=None, sex=None, **kwargs):
    os.chdir(gen_pipeline_file_name())
    print('using SXTRCTR and SCAMP to generate better astrometry...')

    print(f'\nEquivalent argparse cmd: photometrus astrom astromatic -improved -path {astromdir} -band {band}')

    astrometry.astrometry(path=astromdir, band=band, improved=True)


ASTROM_FCTNS = {
    'astrom_angle': astrom_angle,
    'shift': shift,
    'verify_astrom': verify_astrom,
    'verify_astrom_indiv': verify_astrom_indiv,
    'astromatic_astrometry': astromatic_astrometry,
}


def astrom_fctns(commands=None, **kwargs):
    """
    Runs all astrom function steps defined in `commands` (defaults to default_astrom_proc).
    kwargs is a superset of everything any step might need, EACH FUNCTION SHOULD HAVE THE SAME OUTPUT.
    Each step function only pulls out what named vars it needs.
    Anything else falls into its own **kwargs and is ignored.
    """
    if commands is None:
        commands = defaults['default_astrom_proc']
    for cmd in commands:
        ASTROM_FCTNS[cmd](**kwargs)