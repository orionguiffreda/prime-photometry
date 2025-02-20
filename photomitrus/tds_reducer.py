import pandas as pd
import argparse
from photomitrus.multi_master import multi_master


def reduce_all(csvpath, date, chips, bands):
    pointings = pd.read_csv(csvpath)
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
            multi_master(target=field, date=chosen_date, band=band, chip=chips, removal=True)

#%%


def main():
    parser = argparse.ArgumentParser(description='data reducer for large time domain surveys')
    parser.add_argument('-csvpath', type=str, help='[str] input path to csv with all pointings for certain'
                                                   ' date')
    parser.add_argument('-date', type=str, help='[str] date in yyyymmdd, leave blank if csv named w/ format '
                                                'ex. "Durbak_20241224_TDS.csv"')
    parser.add_argument('-chips', type=str, help='[str] specific chip or chips to process')
    parser.add_argument('-bands', type=str, help='[str] specific band or bands to process')
    args, unknown = parser.parse_known_args()

    reduce_all(args.csvpath, args.date, args.chips, args.bands)


if __name__ == "__main__":
    main()
