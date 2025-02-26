"""
Runs photometry on whole observation
"""

import os
import argparse

from photomitrus.photometry import photometry

#%%
defaults = dict(chip=[1,2,3,4])


def fetchstacks(stackpath, chip=None):
    if not chip:
        chips = defaults["chip"]
    else:
        if type(chip) is list:
            chips = chip
        elif type(chip) is int:
            chips = [chip]
        else:
            chips = chip.split(',')
            chips = [int(f) for f in chips]

    allstacks = [f for f in os.listdir(stackpath) if f.endswith('.fits') and f.startswith('coadd.Open-J')]
    matchingstacks = sorted([img for img in allstacks if any(f'C{chosenchips}' in img for chosenchips in chips)])
    if not matchingstacks:
        raise FileNotFoundError('No stacks matching format and given chip(s) are found!')
    else:
        print('Matching stacked images: ', matchingstacks)

    return matchingstacks


def multiphotom(stackpath, matchingstacks, band, survey=None):
    for img in matchingstacks:
        wholeimgpath = os.path.join(stackpath, img)
        print('\nRunning photometry on: %s' % wholeimgpath)
        if survey:
            photometry.photometry(full_filename=wholeimgpath, band=band, survey=survey)
        else:
            photometry.photometry(full_filename=wholeimgpath, band=band)


#%%


def mastermultiphotom(stackpath, band, chip=None, survey=None):
    matchingstacks = fetchstacks(stackpath, chip)
    multiphotom(stackpath, matchingstacks, band, survey)


def main():
    parser = argparse.ArgumentParser(description='Use to run photometry on whole observations (all chips)')
    parser.add_argument('-stackpath', type=str, help='[str], full directory path to stack image directory')
    parser.add_argument('-band', type=str, help='[str] filter, ex. "J"')
    parser.add_argument('-chip', type=str, help='[str] Optional, use to process specific chips, use "1,2,3,4"'
                                                ' format., default = all 4 chips',default=None)
    parser.add_argument('-survey', type=str, help='Specify specific survey to query for photometry (default'
                                                  ' picks for you), see photometrus photometry single -h for list of '
                                                  'available surveys')
    args, unknown = parser.parse_known_args()

    mastermultiphotom(args.stackpath, args.band, args.chip, args.survey)


if __name__ == "__main__":
    main()