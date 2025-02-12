import pandas as pd
import argparse
from photomitrus.multi_master import multi_master


def reduce_all(csvpath, chips, bands):
    pointings = pd.read_csv(csvpath)
    date = csvpath[-16:-8]
    if bands:
        chosen_bands = [bands]
    else:
        chosen_bands = ['J','Z']

    for band in chosen_bands:
        if band == 'Z':
            prunelist = pointings.loc[pointings['Filter1'] == 'Z']
            prunelist = prunelist['FieldNumber'].tolist()
        else:
            prunelist = pointings.loc[pointings['Filter1'] == 'Open']
            prunelist = prunelist['FieldNumber'].tolist()
        for field in prunelist:
            multi_master(target=field, date=date, band=band, chip=chips, removal=True)

#%%


def main():
    parser = argparse.ArgumentParser(description='data reducer for large time domain surveys')
    parser.add_argument('-csvpath', type=str, help='[str] input path to csv with all pointings for certain'
                                                   ' date')
    parser.add_argument('-chips', type=str, help='[str] specific chip or chips to process')
    parser.add_argument('-bands', type=str, help='[str] specific band or bands to process')
    args, unknown = parser.parse_known_args()

    reduce_all(args.csvpath, args.chips, args.bands)


if __name__ == "__main__":
    main()
