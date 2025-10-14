"""
improves astrometry betw. 2 epochs for image sub purposes
"""

import os
import argparse
import shutil
import subprocess
from astropy.io import fits
from astropy.table import Table

from photometrus.stack.stack import swarp_sx, swarp_missfits
from photometrus.settings import gen_config_file_name
from photometrus.photometry.photometry import photometry


def multi_epoch_scamp(input_epoch_cat_path, base_epoch_cat_path):
    """
    Runs scamp on input epoch w/ base epoch's sextractor .cat as a local ref catalog

    Parameters
    ----------
    input_epoch_cat_path: str
        Full filepath to epoch you want to match the base epoch's catalog to
    base_epoch_cat_path: str
        Full filepath to sextractor catalog of base epoch
    """

    sc = gen_config_file_name('scamp.conf')
    command = (f'scamp {input_epoch_cat_path} -c {sc} -ASTREF_CATALOG FILE -ASTREFCAT_NAME {base_epoch_cat_path} '
               f'-ASTREFMAG_LIMITS -99.0,99.0 '
               f'-ASTREFCENT_KEYS ALPHA_J2000,DELTA_J2000 -ASTREFERR_KEYS ERRAWIN_WORLD,ERRBWIN_WORLD,ERRTHETAWIN_WORLD '
               f'-ASTREFMAG_KEY MAG_AUTO -ASTREFMAGERR_KEY MAGERR_AUT')

    subprocess.run(command.split(), check=True)


def multi_epoch_astrom(base_epoch_path, matching_epoch_path):
    """
    Runs sextractor on base epoch stacked image, uses the resulting .cat file as a supplier to
    scamp, which is run on the other epoch.  This is in order to match astrometry between both.
    The matching epoch image is copied to the base epoch directory, and all processes are run there.

    Parameters
    ----------
    base_epoch_path: str
        Full filepath to epoch you want to base the astrometry on
    matching_epoch_path: str
        Full filepath to epoch you want to match the base epoch astrometry to
    """

    base_epoch_dir, base_epoch_name = os.path.split(base_epoch_path)
    os.chdir(base_epoch_dir)
    base_epoch_hdr = fits.getheader(base_epoch_path)
    chip = base_epoch_hdr['CHIP']
    # if len(base_epoch_hdr['FILTER2']) > 1:
    #     band = 'Z'
    # else:
    #     band = base_epoch_hdr['FILTER2']

    match_epoch_dir, match_epoch_name = os.path.split(matching_epoch_path)
    match_epoch_new_path = os.path.join(base_epoch_dir, match_epoch_name)
    shutil.copyfile(matching_epoch_path, match_epoch_new_path)

    # print(f'Running photometry on {base_epoch_path}')
    # photometry(full_filename=base_epoch_path, band=band)
    # base_ecsv = [f for f in os.listdir(base_epoch_dir) if f.endswith('.ecsv') and base_epoch_name in f]
    # base_ecsv_path = os.path.join(base_epoch_dir, base_ecsv[0])
    # base_cat = Table.read(base_ecsv_path)
    # base_cat_path = base_ecsv_path.replace('.ecsv', '.fits')
    # base_cat.write(base_cat_path, format='fits', overwrite=True)
    #
    # print(f'\nRunning photometry on {match_epoch_new_path}')
    # photometry(full_filename=match_epoch_new_path, band=band)
    # match_ecsv = [f for f in os.listdir(base_epoch_dir) if f.endswith('.ecsv') and match_epoch_name in f]
    # match_ecsv_path = os.path.join(base_epoch_dir, match_ecsv[0])
    # match_cat = Table.read(match_ecsv_path)
    # match_cat_path = match_ecsv_path.replace('.ecsv', '.fits')
    # match_cat.write(match_cat_path, format='fits', overwrite=True)

    print(f'Sextracting base epoch: {base_epoch_name}...')
    base_cat_path = swarp_sx(imgpath=base_epoch_path, chip=chip)
    print(f'Sextracting matching epoch: {match_epoch_name}...')
    match_cat_path = swarp_sx(imgpath=match_epoch_new_path, chip=chip)
    multi_epoch_scamp(input_epoch_cat_path=match_cat_path, base_epoch_cat_path=base_cat_path)
    swarp_missfits(imgpath=match_epoch_new_path, chip=chip)


def main():
    parser = argparse.ArgumentParser(description='improves astrometry betw. 2 epochs for image sub purposes, all work'
                                                 'is done in the base epoch directory')
    parser.add_argument(
        '-base_epoch', type=str, help='[str] Full filepath to epoch you want to base the astrometry on')
    parser.add_argument(
        '-match_epoch', type=str, help='[str] Full filepath to epoch you want to match the base epoch astrometry to')
    args, unknown = parser.parse_known_args()

    multi_epoch_astrom(base_epoch_path=args.base_epoch, matching_epoch_path=args.match_epoch)


if __name__ == "__main__":
    main()
