import pandas as pd
import os
from photometrus.settings import gen_pipeline_file_name, gen_flat_dir
from photometrus.preprocess import gen_flat
from pandas import to_datetime

import argparse
import sys

#%%


def auto_flat_creation(date, band=None, chip=None):
    datetime = to_datetime(date)
    date = datetime.strftime('%Y%m%d')

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


def main():
    parser = argparse.ArgumentParser(description='Downloads data and generates master flats for an observation, '
                                                 'storing them in prime-photometry')
    parser.add_argument('-date', type=str, help='[str] date of observation, in yyyymmdd format')
    parser.add_argument('-band', type=str, help='[str], optional, specify filter ex. "J", otherwise it will'
                                                ' try to generate mflats regardless of filter', default=None)
    parser.add_argument('-chip', type=int, help='[int], optional, to only generate mflat for 1 detector',
                        default=None)
    args, unknown = parser.parse_known_args()

    autoflatgen(args.date, args.band, args.chip)


if __name__ == "__main__":
    main()
