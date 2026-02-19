"""
Runs photometry on whole observation
"""

import os
import argparse

from photometrus.photometry import photometry
from photometrus.photometry import ellipticity_logger
from photometrus.utils.defaults import PROCESSING_DEFAULTS as defaults

#%%


def fetchstacks(stackpath, chip=None):
    if not chip:
        chips = defaults['chip']
    else:
        if type(chip) is list:
            chips = chip
        elif type(chip) is int:
            chips = [chip]
        else:
            chips = chip.split(',')
            chips = [int(f) for f in chips]

    allstacks = [f for f in os.listdir(stackpath) if f.endswith('.fits') and f.startswith('coadd.')]
    matchingstacks = sorted([img for img in allstacks if any(f'C{chosenchips}' in img for chosenchips in chips)])
    if not matchingstacks:
        raise FileNotFoundError('No stacks matching format and given chip(s) are found!')
    else:
        print('Matching stacked images: ', matchingstacks)

    return matchingstacks


def multiphotom(stackpath, matchingstacks, band, survey, grb_ra, grb_dec,
                grb_coordlist, grb_radius, grb_only, grb_name, no_int_cal, keep, det_cut, no_plots):
    for img in matchingstacks:
        wholeimgpath = os.path.join(stackpath, img)
        print('\nRunning photometry on: %s' % wholeimgpath)
        photometry.photometry(full_filename=wholeimgpath, band=band, survey=survey, grb_ra=grb_ra, grb_dec=grb_dec,
                           grb_coordlist=grb_coordlist, grb_radius=grb_radius, grb_only=grb_only, grb_name=grb_name, no_int_cal=no_int_cal,
                              keep=keep, det_cut=det_cut, no_plots=no_plots)

    if len(matchingstacks) == 4:
        ellipticity_logger.logger(directory=stackpath)
    else:
        pass

#%%


def mastermultiphotom(stackpath=defaults['stackpath'], band=defaults['band'], chip=defaults['chip'], survey=defaults['survey'],
                      grb_ra=defaults['grb_ra'], grb_dec=defaults['grb_dec'], grb_coordlist=defaults['grb_coordlist'],
                      grb_radius=defaults['grb_radius'], grb_only=defaults['grb_only'], grb_name=defaults['grb_name'], no_int_cal=defaults['no_int_cal'],
                      keep=defaults['keep'], det_cut=defaults['det_cut'], no_plots=defaults['no_plots']):

    cmd_str = f'\nEquivalent argparse cmd: photometrus photometry -stackpath {stackpath} -band {band} -chip {chip}'
    if grb_ra:
        grb_str = f' -grb_ra {grb_ra} -grb_dec {grb_dec} -grb_radius {grb_radius}'
    else:
        grb_str = ''
    print(cmd_str + grb_str)

    matchingstacks = fetchstacks(stackpath, chip)
    multiphotom(stackpath, matchingstacks, band, survey, grb_ra, grb_dec, grb_coordlist, grb_radius, grb_only, grb_name, no_int_cal,
                keep, det_cut, no_plots)


def main():
    parser = argparse.ArgumentParser(description='Use to run photometry on whole observations (all chips)')
    parser.add_argument('-stackpath', type=str, help='[str], full directory path to stack image directory',
                        default=defaults['stackpath'])
    parser.add_argument('-band', type=str, help='[str] filter, ex. "J"', default=defaults['band'])
    parser.add_argument('-chip', type=str, help='[str] Optional, use to process specific chips, use "1,2,3,4"'
                                                ' format., default = all 4 chips',default=defaults['chip'])
    parser.add_argument('-survey', type=str, help='Specify specific survey to query for photometry (default'
                                                  ' picks for you), see photometrus photometry single -h for list of '
                                                  'available surveys', default=defaults['survey'])
    parser.add_argument('-grb_ra', type=str, help='[str], RA for GRB source, either in hh:mm:ss or decimal'
                                                  '*NOTE* When using sexagesimal, use "-grb_ra=value_here" NOT "-grb_ra '
                                                  'value_here", as argparse doesnt like negative sexagesimals',
                        default=defaults['grb_ra'])
    parser.add_argument('-grb_dec', type=str, help='[str], DEC for GRB source, either in dd:mm:ss or decimal'
                                                   '*NOTE* When using sexagesimal, use "-grb_dec=value_here" NOT "-grb_dec '
                                                  'value_here", as argparse doesnt like negative sexagesimals',
                        default=defaults['grb_dec'])
    parser.add_argument('-grb_coordlist', type=str, nargs='+',
                        help='[float] Used to check multiple GRB locations.  Input RA and DECs of locations '
                             'with the format: -coordlist 123,45 -123,-45 etc..  *DONT USE -RA '
                             '& -DEC BUT INCLUDE -grb_radius*', default=defaults['grb_coordlist'])
    parser.add_argument('-grb_radius', type=str,
                        help='[str], # of arcsec diameter to search for GRB, default = 4.0".  You can specify arcsec,'
                             ' arcmin, or deg w/ an underscore.  Ex. "-grb_radius 3_arcmin" will specify an area of 3 '
                             'arcminutes.  If just a number is applied, it defaults to arcsec.',
                        default=defaults['grb_radius'])
    parser.add_argument('-grb_name', type=str,
                        help='[str] optional name for grb-related data products, default = "GRB"',
                        default=defaults["grb_name"])
    parser.add_argument('-grb_only', action='store_true',
                        help='optional flag, use if running -grb again on already created catalog',
                        default=defaults['grb_only'])
    parser.add_argument('-no_int_cal', action='store_true',
                        help='optional flag, use to STOP photometric fit intercept calibration, taking only the initial calc.',
                        default=defaults["no_int_cal"])
    parser.add_argument('-no_plots', action='store_true',
                        help='optional flag, stops creation of mag comparison plot betw. PRIME and survey, '
                             'along with residual plot w/ statistics, lim mag plot', default=defaults['no_plots'])
    parser.add_argument('-keep', action='store_true',
                        help='optional flag, use if you DONT want to remove intermediate products after getting photom,'
                             ' i.e. the ".cat" and ".psf" files', default=defaults['keep'])
    parser.add_argument('-det_cut', type=float, help='[float], num of median image sigma to cut off sources'
                                                     ' (ex. det_thresh of 2 => cutoff = med - 2*sigma',
                        default=defaults["det_cut"])

    args, unknown = parser.parse_known_args()
    # print(args)

    mastermultiphotom(args.stackpath, args.band, args.chip, args.survey, args.grb_ra, args.grb_dec, args.grb_coordlist,
                      args.grb_radius, args.grb_only, args.grb_name, args.no_int_cal, args.keep, args.det_cut, args.no_plots)


if __name__ == "__main__":
    main()