import pandas as pd
import os
import subprocess
from photomitrus.settings import gen_pipeline_file_name
from photomitrus.preprocess import gen_flat
import shutil
import argparse

#%%


def flatdatadownload(directory, band, date, chip=None):
    print('Creating parent dir: %s' % directory)
    if not os.path.isdir(directory):
        os.mkdir(directory)
    else:
        pass
    if chip:
        if band == 'Z':
            try:
                command = ('python ./getdata.py -d %s -i 192.168.212.22 -f ramp --objname FLAT '
                           '-c %s --filter1 Z --filter2 Open %s') % (directory, chip, date)
                print('Executing command: %s' % command)
                rval = subprocess.run(command.split(), check=True)
            except subprocess.CalledProcessError as err:
                print('Could not run with exit error %s' % err)
        else:
            try:
                command = ('python ./getdata.py -d %s -i 192.168.212.22 -f ramp --objname FLAT '
                           '-c %s --filter1 Open --filter2 %s %s') % (directory, chip, band, date)
                print('Executing command: %s' % command)
                rval = subprocess.run(command.split(), check=True)
            except subprocess.CalledProcessError as err:
                print('Could not run with exit error %s' % err)
    else:
        if band == 'Z':
            try:
                command = ('python ./getdata.py -d %s -i 192.168.212.22 -f ramp --objname FLAT '
                           '-c 1,2,3,4 --filter1 Z --filter2 Open %s') % (directory, date)
                print('Executing command: %s' % command)
                rval = subprocess.run(command.split(), check=True)
            except subprocess.CalledProcessError as err:
                print('Could not run with exit error %s' % err)
        else:
            try:
                command = ('python ./getdata.py -d %s -i 192.168.212.22 -f ramp --objname FLAT '
                           '-c 1,2,3,4 --filter1 Open --filter2 %s %s') % (directory, band, date)
                print('Executing command: %s' % command)
                rval = subprocess.run(command.split(), check=True)
            except subprocess.CalledProcessError as err:
                print('Could not run with exit error %s' % err)


def auto_flat_creation(directory, band, date, chip=None):
    if not chip:
        chips = [1,2,3,4]
    else:
        chips = [chip]

    mflat_storage_dir = gen_pipeline_file_name() + '/mflats/'

    flat_filter_list = []
    mflat_paths = []
    for f in chips:

        print('\nBeginning mflat gen for chip %s' % f)

        start_images_names_1, start_images_names_2,end_images_names_1, end_images_names_2, flat_filter = gen_flat.flatlists(directory,f)
        save_name_start, save_name_end = gen_flat.flatprocessing(directory, start_images_names_1, start_images_names_2,
                                                        end_images_names_1, end_images_names_2)
        if save_name_start:
            save_name = save_name_start
        elif save_name_end:
            save_name = save_name_end
        else:
            print('No flats on night!')

        flat_filter_list.append(flat_filter)

        datename = 'mflat.%s.%s.C%s.fits' % (band, date, f)
        # datename_end = 'mflat.end.%s.%s.fits' % (band, date)

        current_mflat_path = os.path.join(directory + 'mflats/', save_name)

        mflatpath = os.path.join(directory + 'mflats/', datename)
        # mflatpath_end = os.path.join(directory + 'mflats/', datename_end)

        os.rename(current_mflat_path, mflatpath)
        # os.rename(os.path.join(directory + 'mflats/', save_name_end), mflatpath_end)
        mflat_paths.append(mflatpath)

        print('Storing mflat file: %s...' % datename)

        pipeline_mflat_path = os.path.join(mflat_storage_dir, datename)
        shutil.copyfile(mflatpath, pipeline_mflat_path)

    return flat_filter_list[0], mflat_paths

"""
def mflat_storage_and_pick(flat_filter, mflat_paths, band, chip):

    if band == flat_filter:
        print('Tonight"s flats taken in matching filter! \nProceeding with master flat...')

        mflat_for_obs = current_mflat_path

    else:
        print('Not taken in matching filter... Defaulting to most recent master flat!')
        mflat_list = [f for f in os.listdir(mflat_storage_dir)]
"""
#%%


def autoflatgen(directory, date, band, chip=None, no_download=False):
    if no_download:
        auto_flat_creation(directory, band, date, chip)
    else:
        flatdatadownload(directory, band, date, chip)
        auto_flat_creation(directory, band, date, chip)


def main():
    parser = argparse.ArgumentParser(description='Downloads data and generates master flats for an observation, '
                                                 'storing them in prime-photometry')
    parser.add_argument('-no_download', action='store_true', help='optional flag, use if you already have'
                                                                  ' the data')
    parser.add_argument('-dir', type=str, help='[str], directory')
    parser.add_argument('-date', type=str, help='[str] date of observation, in yyyymmdd format')
    parser.add_argument('-band', type=str, help='[str] filter, ex. "J"')
    parser.add_argument('-chip', type=int, help='[int], optional, to only generate mflat for 1 detector',
                        default=None)
    args, unknown = parser.parse_known_args()
    autoflatgen(args.dir, args.date, args.band, args.chip, args.no_download)


if __name__ == "__main__":
    main()
