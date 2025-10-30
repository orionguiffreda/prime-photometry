"""
Functions for automatic master / super flat generation
"""

import pandas as pd
import os
import sys
from pandas import to_datetime
from datetime import datetime, timedelta
import calendar
import argparse
from astropy.io import fits
import numpy as np

from photometrus.settings import gen_pipeline_file_name, gen_flat_dir, mflat_checker, gen_mflat_file_name
from photometrus.preprocess import gen_flat
from photometrus.getfiles import get_data_files

#%%


def auto_flat_creation(date, band=None, chip=None):
    try:
        datetime = to_datetime(date)
        date = datetime.strftime('%Y%m%d')
    except AttributeError:
        sys.exit('Date of observation not extracted, cannot move forward with flat gen & processing! '
                 'Did you specify a "-date" field correctly?')

    if not chip:
        chips = [1,2,3,4]
    else:
        chips = [chip]

    # mflat_storage_dir = gen_pipeline_file_name() + '/mflats/'
    mflat_storage_dir = gen_flat_dir() + os.path.sep

    for f in chips:

        print('\nBeginning mflat gen for chip %s' % f)

        save_name_start, save_name_end, flat_filter = gen_flat.flatgen(mflat_storage_dir, date, f, band)

        if save_name_start:
            save_name = save_name_start
        elif save_name_end:
            save_name = save_name_end
        else:
            print('Will take closest mflat during processing, continuing..')
            save_name = 0

        if save_name:
            datename = 'mflat.%s.%s.C%s.fits' % (flat_filter, date, f)
            datename_end = 'mflat.end.%s.%s.C%s.fits' % (flat_filter, date, f)

            current_chosen_path = os.path.join(mflat_storage_dir, save_name)
            renamed_chosen_path = os.path.join(mflat_storage_dir, datename)
            os.rename(current_chosen_path, renamed_chosen_path)
            print('Renamed to', renamed_chosen_path)

            if save_name_start:
                if save_name_end:
                    current_auxiliary_path = os.path.join(mflat_storage_dir, save_name_end)
                    renamed_auxiliary_path = os.path.join(mflat_storage_dir+'/auxiliary_mflats/', datename_end)
                    os.rename(current_auxiliary_path, renamed_auxiliary_path)
                    print('Auxiliary mflat:', renamed_auxiliary_path)
                else:
                    print('No auxiliary mflat, no flat data at end of night.')
        else:
            pass


"""
        current_mflat_path = os.path.join(directory + 'mflats/', save_name)

        mflatpath = os.path.join(directory + 'mflats/', datename)
        # mflatpath_end = os.path.join(directory + 'mflats/', datename_end)

        os.rename(current_mflat_path, mflatpath)
        # os.rename(os.path.join(directory + 'mflats/', save_name_end), mflatpath_end)
        mflat_paths.append(mflatpath)

        print('Storing mflat file: %s...' % datename)

        pipeline_mflat_path = os.path.join(mflat_storage_dir, datename)
        shutil.copyfile(mflatpath, pipeline_mflat_path)
"""
#%%


def autoflatgen(date, band=None, chip=None):
    auto_flat_creation(date, band, chip)


def new_mflat_checker(date, band=None, sflat=False):
    base_dir = gen_pipeline_file_name()
    flat_dir = gen_flat_dir()
    # flat_dir = '/home/alex/PycharmProjects/prime-photometry/photometrus/mflats/'
    if not band:
        mflat_list = [
            f for f in sorted(os.listdir(flat_dir)) if f.endswith('.fits') if '.%s.' % date in f]

        if not mflat_list:
            print('\nNo master flat currently generated for this date!')
            return False
        else:
            print('\nMaster flats exist for this date!')
        return True

    else:
        # checking for any current mflats
        mflat_list = [
            f for f in sorted(os.listdir(flat_dir)) if f.endswith('.fits') if '.%s.%s.' % (band, date) in f]

        # checking if flat data was taken in band
        try:
            obs_flatlists, m_lists = get_data_files(date=date, objname='FLAT', filter1='Open',
                                                filter2=band)

            if any(lst for lst in m_lists):
                print('Missing files!: ', m_lists)
        except FileNotFoundError:
            obs_flatlists = [[],[]]

        if mflat_list or all(len(sub) != 0 for sub in obs_flatlists):
            if sflat:
                print('Master flats exist or flat data was taken for this date and band!')
            return True
        else:
            if sflat:
                print('No master flat currently generated or flat data taken for this date and band!')
            return False


def gen_mflat_file_name_new(band, chip, date=None):
    flat_dir = gen_flat_dir()
    mflat_list = [
        f for f in sorted(os.listdir(flat_dir)) if f.endswith('.fits') if '.%s.' % band in f if 'C%s' % chip in f
        if len(f) <= 24]

    if date:

        # get mflat closest to obs date, if there is a tie, it picks the earlier one to be safe
        def closest_file(file_list, target_date):
            target = datetime.strptime(target_date, "%Y%m%d")

            def extract_date(file):
                return datetime.strptime(file.split(".")[2], "%Y%m%d")

            filename = min(file_list, key=lambda f: (abs((extract_date(f) - target).days), extract_date(f)))
            print('Getting mflat closest to given date: ', filename)
            return filename

        # filename = closest_file(mflat_list, date)
        # print('Getting mflat closest to given date: ', filename)

        def iterate_until_mflat(date):
            target_date = datetime.strptime(date, "%Y%m%d")

            if new_mflat_checker(date, band=band):
                chosen_date = date
            else:
                offset = 1
                chosen_date = None
                while chosen_date is None:
                    for delta in (-offset, offset):
                        check_date = target_date + timedelta(days=delta)
                        check_date_str = check_date.strftime("%Y%m%d")
                        print('\nChecking %s for existing mflat in %s band...' % (check_date_str, band) )
                        if new_mflat_checker(check_date_str, band=band):
                            chosen_date = check_date_str
                            break
                    offset += 1

            filename = [m for m in mflat_list if '.%s.%s.' % (band, chosen_date) in m]
            if filename:
                print('Matching mflat found! : ', filename)
            else:
                autoflatgen(chosen_date)

            new_mflat_list = [
                f for f in sorted(os.listdir(flat_dir)) if f.endswith('.fits') if '.%s.' % band in f if
                                                                                  'C%s' % chip in f
                if len(f) <= 24]
            filename = [m for m in new_mflat_list if '.%s.%s.C%s' % (band, chosen_date, chip) in m]
            return filename[0]

        filename = iterate_until_mflat(date)

    else:
        # if no date given, get latest mflat
        mflat_list = sorted(mflat_list, reverse=True)
        filename = mflat_list[0]
        print('No date, getting latest mflat: ', filename)
    return os.path.join(flat_dir, filename)


# generates a super master flat from all mfs in a date's month for that band
def superflatgen(date, band, chip=None):
    flat_dir = gen_flat_dir() + os.path.sep

    try:
        # Try single-date format
        dt = datetime.strptime(date, "%Y%m%d")
        first_day = dt.replace(day=1)
        last_day = dt.replace(day=calendar.monthrange(dt.year, dt.month)[1])
    except ValueError:
        # Try date-range format 'yyyymmdd-yyyymmdd'
        try:
            start_str, end_str = date.split('-')
            first_day = datetime.strptime(start_str, "%Y%m%d")
            last_day = datetime.strptime(end_str, "%Y%m%d")
        except Exception as e:
            raise ValueError(f"Date format not recognized: {date}") from e

    print('Searching for all existing mflats between %s - %s' % (first_day, last_day))

    matched_dates = []
    current = first_day
    while current <= last_day:
        check_date_str = current.strftime("%Y%m%d")
        if mflat_checker(check_date_str, band=band, sflat=True):
            matched_dates.append(check_date_str)
        current += timedelta(days=1)

    print('%s band mflats found: %s' % (band, matched_dates))
    all_matches = []
    for date in matched_dates:
        if chip:
            mflat = gen_mflat_file_name(band, chip, date=date)
            all_matches.append(mflat)
        else:
            for i in range(1,5):
                mflat = gen_mflat_file_name(band, i, date=date)
                all_matches.append(mflat)

    if chip:
        matches = [mf for mf in all_matches if '.C%s' % chip in mf]
        match_paths = [os.path.join(flat_dir, mf) for mf in matches]
        match_data = [fits.getdata(img) for img in match_paths]
        fdata_stack = np.stack(match_data)
        super_flat_data = np.nanmedian(fdata_stack, axis=0)

        start_date = matches[0].split('.')[3]
        end_date = matches[-1].split('.')[3]
        sfname = f'sflat.{band}.{start_date}-{end_date}.C{i}.fits'
        fits.writeto(os.path.join(flat_dir, sfname), super_flat_data, overwrite=True)

        print('sflat generated at: ',os.path.join(flat_dir, sfname))

    else:
        for i in range(1, 5):
            matches = [mf for mf in all_matches if '.C%s' % i in mf]
            match_paths = [os.path.join(flat_dir, mf) for mf in matches]
            match_data = [fits.getdata(img) for img in match_paths]
            fdata_stack = np.stack(match_data)
            super_flat_data = np.nanmedian(fdata_stack, axis=0)

            start_date = matches[0].split('.')[3]
            end_date = matches[-1].split('.')[3]
            sfname = f'sflat.{band}.{start_date}-{end_date}.C{i}.fits'
            fits.writeto(os.path.join(flat_dir, sfname), super_flat_data, overwrite=True)

            print('sflat generated at: ', os.path.join(flat_dir, sfname))




def main():
    parser = argparse.ArgumentParser(description='Downloads data and generates master flats for an observation, '
                                                 'storing them in prime-photometry')
    parser.add_argument('-sflat', action='store_true', help='use to generate super flat (median of all mflats'
                                                            ' from that month)')
    parser.add_argument('-date', type=str, help='[str] date of observation, in yyyymmdd format, or when '
                                                'using -sflat, can be 2 dates to generate sflats betw., w/ the format '
                                                '"yyyymmdd-yyyymmdd"')
    parser.add_argument('-band', type=str, help='[str], optional, specify filter ex. "J", otherwise it will'
                                                ' try to generate mflats regardless of filter', default=None)
    parser.add_argument('-chip', type=int, help='[int], optional, to only generate mflat for 1 detector',
                        default=None)
    args, unknown = parser.parse_known_args()

    if args.sflat:
        superflatgen(args.date, args.band, args.chip)
    else:
        autoflatgen(args.date, args.band, args.chip)


if __name__ == "__main__":
    main()
