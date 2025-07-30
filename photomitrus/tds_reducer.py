import pandas as pd
import argparse
import os

from photomitrus.multi_master import (multi_master, parentcreation)
from photomitrus.photometry.multi_photom import mastermultiphotom


def reduce_and_photom_all(csvpath, date, bands, chips=None, reduction_only=False, photometry_only=False):
    pointings = pd.read_csv(csvpath)
    if not chips:
        chosenchips = '1,2,3,4'
    else:
        chosenchips = chips

    if date:
        chosen_date = date
    else:
        chosen_date = csvpath[-16:-8]

    if bands:
        chosen_bands = [bands]
    else:
        chosen_bands = ['J','Z']

    for band in chosen_bands:
        if 'Filter1' in pointings:
            if band == 'Z':
                prunelist = pointings.loc[pointings['Filter1'] == 'Z']
                prunelist = prunelist['FieldNumber'].tolist()
            else:
                prunelist = pointings.loc[pointings['Filter1'] == 'Open']
                prunelist = prunelist['FieldNumber'].tolist()
        else:
            prunelist = pointings['FieldNumber'].tolist()
        for field in prunelist:
            if reduction_only:
                multi_master(target=field, date=chosen_date, band=band, chip=chosenchips, removal=True)
            elif photometry_only:
                parent_dir, field_dir = parentcreation(target=field, date=chosen_date, band=band)
                stackpath = os.path.join(parent_dir, 'stack')
                print('Stack directory: ', stackpath)
                mastermultiphotom(stackpath=stackpath, band=band, chip=chosenchips)
            else:
                multi_master(target=field, date=chosen_date, band=band, chip=chosenchips, removal=True)
                parent_dir, field_dir = parentcreation(target=field, date=chosen_date, band=band)
                stackpath = os.path.join(parent_dir, 'stack')
                print('Stack directory: ', stackpath)
                mastermultiphotom(stackpath=stackpath, band=band, chip=chosenchips)

#%%


def main():
    parser = argparse.ArgumentParser(description='data reducer / photometric runner for large time domain surveys')
    parser.add_argument('-csvpath', type=str, help='[str] input path to csv with all pointings for certain'
                                                   ' date: csv should have simple format of 1 column w/ column name '
                                                   '"FieldNumber" and field names below in column')
    parser.add_argument('-date', type=str, help='[str] date in yyyymmdd, leave blank if csv named w/ format '
                                                'ex. "Durbak_20241224_TDS.csv"')
    parser.add_argument('-chips', type=str, help='[str] specific chip or chips to process')
    parser.add_argument('-bands', type=str, help='[str] specific band or bands to process')
    parser.add_argument('-reduction_only', action='store_true', help='optional flag, use to reduce all fields'
                                                                     ' W/O running photometry')
    parser.add_argument('-photometry_only', action='store_true', help='optional flag, use to run photometry '
                                                                      'on all fields when ALREADY reduced')
    args, unknown = parser.parse_known_args()  # TODO: get default arguments from defaults dict

    reduce_and_photom_all(args.csvpath, args.date, args.bands, args.chips, args.reduction_only, args.photometry_only)


if __name__ == "__main__":
    main()
