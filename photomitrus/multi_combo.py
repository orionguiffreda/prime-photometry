from photomitrus import multi_master
from photomitrus.photometry import multi_photom
from photomitrus.photometry import ellipticity_logger
from photomitrus.photometry.photometry import defaults
from photomitrus.settings import PIPELINE_DEFAULT_DIR

import os
import shutil
import argparse

#%%


def combo(target, date, band, chip=None, parentdir=False, rot_val=48, no_shift=False, astromnet=False,
          sky_override_path=False, removal=False, no_get_files=False, no_download=False, no_mflat=False, bulge=False,
          survey=None, grb_ra=defaults['RA'], grb_dec=defaults['DEC'], grb_coordlist=None, grb_radius=defaults['thresh'],
          auto_mode=False, input_ramp_lists=None
          ):

    if parentdir:
        chosen_parent = parentdir
    else:
        chosen_parent, field_dir = multi_master.parentcreation(target, date, band)
    stackpath = os.path.join(chosen_parent, 'stack')

    if not chip:
        chips = [1, 2, 3, 4]
    else:
        chips = chip.split(',')
        chips = [int(f) for f in chips]

    for f in chips:
        multi_master.multi_master(target, date, band, f, chosen_parent, rot_val, no_shift, astromnet, no_download,
                                  sky_override_path, removal, no_get_files, no_mflat, bulge, auto_mode=auto_mode,
                                  input_ramp_lists=input_ramp_lists)
        multi_photom.mastermultiphotom(stackpath, band, f, survey, grb_ra, grb_dec, grb_coordlist, grb_radius)

    if len(chips) == 4:
        catcheck = [f for f in sorted(os.listdir(stackpath)) if f.endswith('.ecsv') and f.startswith('coadd')]
        if len(catcheck) == 4:
            ellipticity_logger.logger(directory=stackpath)

#%%


def main():
    parser = argparse.ArgumentParser(description='Use to process and run photometry on whole observations (all chips).'
                                                 ' Can use all options (except -parallel & file download) from photomitrus '
                                                 'pipeline & photomitrus photometry')
    parser.add_argument('-target', type=str, help='[str] target field, objname in log, ex. "field1234"')
    parser.add_argument('-date', type=str, help='[str] date of observation, in yyyymmdd format')
    parser.add_argument('-band', type=str, help='[str] filter, ex. "J"')
    parser.add_argument('-chip', type=str,
                        help='[str] Optional, use to process specific chips, use "1,2,3,4"'
                             ' format.', default=None)
    parser.add_argument('-no_get_files', action='store_true', help='optional flag, use if you want to download '
                                                                     'through old method (scp), new method passes file '
                                                                     'paths (new saves space and time)')
    parser.add_argument('-no_download', action='store_true', help='optional flag, use if you *ALREADY* have the data'
                                                                  'downloaded, *NOT* to use new file path method')
    parser.add_argument('-no_shift', action='store_true', help='optional flag, DO NOT use astrometric shift'
                                                               ' script in place of astrom.net, will not use either (shift is default)')
    parser.add_argument('-astromnet', action='store_true',
                        help='optional flag, use astrom.net to reinforce astrometry')
    parser.add_argument('-removal', action='store_true',
                        help='optional flag, used to remove intermediate subdirectories and data, leaving only the '
                             'stacks & skies; intended for space saving in large nights of observation')
    parser.add_argument('-parent', type=str, help='[str] *NOW OPTIONAL* specify parent directory to '
                                                  'store all data products, otherwise it will automatically generate w/'
                                                  'the format "/target_date/band/"', default=None)
    parser.add_argument('-rot_val', type=float, help='[float] optional, put in your rot angle in deg,'
                                                     ' if you had a non-default rotation angle in your obs'
                                                     ' (default = 48 deg or 172800")', default=48)
    parser.add_argument('-sky_override', type=str, help='[str], Optional path to specify already generated '
                                                        'sky to use in sky sub, skipping sky gen. Input full file path.',
                        default=None)
    parser.add_argument('-no_mflat', action='store_true', help='optional flag, use if you *DO NOT* want to'
                                                               ' automatically generate mflats for this night if none'
                                                               ' exist')
    parser.add_argument('-bulge', action='store_true',
                        help='optional flag, utilize setup specifically designed for bulge fields.  Hopefully we can'
                             ' automate this in the future')
    parser.add_argument('-survey', type=str, help='Specify specific survey to query for photometry (default'
                                                  ' picks for you), see photometrus single_photometry -h for list of '
                                                  'available surveys')
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

    combo(args.target, args.date, args.band, args.chip, args.parent, args.rot_val, args.no_shift, args.astromnet,
          args.sky_override, args.removal, args.no_get_files, args.no_download, args.no_mflat, args.bulge, args.survey,
          args.grb_ra, args.grb_dec, args.grb_coordlist, args.grb_radius)


if __name__ == "__main__":
    main()

