"""
Runs photometry on whole observation
"""

import os
import argparse

from photomitrus.photometry import photometry
from photomitrus.photometry import ellipticity_logger

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

    allstacks = [f for f in os.listdir(stackpath) if f.endswith('.fits') and f.startswith('coadd.Open')]
    matchingstacks = sorted([img for img in allstacks if any(f'C{chosenchips}' in img for chosenchips in chips)])
    if not matchingstacks:
        raise FileNotFoundError('No stacks matching format and given chip(s) are found!')
    else:
        print('Matching stacked images: ', matchingstacks)

    return matchingstacks


def multiphotom(stackpath, matchingstacks, band, survey=None, grb_ra=None, grb_dec=None, grb_coordlist=None,
                grb_radius=None):
    for img in matchingstacks:
        wholeimgpath = os.path.join(stackpath, img)
        print('\nRunning photometry on: %s' % wholeimgpath)
        photometry.photometry(full_filename=wholeimgpath, band=band, survey=survey, grb_ra=grb_ra, grb_dec=grb_dec,
                           grb_coordlist=grb_coordlist, grb_radius=grb_radius)

    ellipticity_logger.logger(directory=stackpath)

#%%


def mastermultiphotom(stackpath, band, chip=None, survey=None, grb_ra=None, grb_dec=None, grb_coordlist=None,
                grb_radius=None):
    matchingstacks = fetchstacks(stackpath, chip)
    multiphotom(stackpath, matchingstacks, band, survey, grb_ra, grb_dec, grb_coordlist, grb_radius)


def main():
    parser = argparse.ArgumentParser(description='Use to run photometry on whole observations (all chips)')
    parser.add_argument('-stackpath', type=str, help='[str], full directory path to stack image directory')
    parser.add_argument('-band', type=str, help='[str] filter, ex. "J"')
    parser.add_argument('-chip', type=str, help='[str] Optional, use to process specific chips, use "1,2,3,4"'
                                                ' format., default = all 4 chips',default=None)
    parser.add_argument('-survey', type=str, help='Specify specific survey to query for photometry (default'
                                                  ' picks for you), see photometrus photometry single -h for list of '
                                                  'available surveys', default=None)
    parser.add_argument('-grb_ra', type=str, help='[str], RA for GRB source, either in hh:mm:ss or decimal'
                                                  '*NOTE* When using sexagesimal, use "-grb_ra=value_here" NOT "-grb_ra '
                                                  'value_here", as argparse doesnt like negative sexagesimals',
                        default=None)
    parser.add_argument('-grb_dec', type=str, help='[str], DEC for GRB source, either in dd:mm:ss or decimal'
                                                   '*NOTE* When using sexagesimal, use "-grb_dec=value_here" NOT "-grb_dec '
                                                  'value_here", as argparse doesnt like negative sexagesimals',
                        default=None)
    parser.add_argument('-grb_coordlist', type=str, nargs='+',
                        help='[float] Used to check multiple GRB locations.  Input RA and DECs of locations '
                             'with the format: -coordlist 123,45 -123,-45 etc..  *DONT USE -RA '
                             '& -DEC BUT INCLUDE -grb_radius*', default=None)
    parser.add_argument('-grb_radius', type=str,
                        help='[str], # of arcsec diameter to search for GRB, default = 4.0".  You can specify arcsec,'
                             ' arcmin, or deg w/ an underscore.  Ex. "-grb_radius 3_arcmin" will specify an area of 3 '
                             'arcminutes.  If just a number is applied, it defaults to arcsec.',
                        default='4.0')
    args, unknown = parser.parse_known_args()

    mastermultiphotom(args.stackpath, args.band, args.chip, args.survey, args.grb_ra, args.grb_dec, args.grb_coordlist,
                      args.grb_radius)


if __name__ == "__main__":
    main()