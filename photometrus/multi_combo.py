import os
import shutil
import argparse

from photometrus import multi_master
from photometrus.photometry import multi_photom
from photometrus.photometry import ellipticity_logger
from photometrus.settings import PIPELINE_DEFAULT_DIR
from photometrus.utils.defaults import PROCESSING_DEFAULTS as defaults

#%%


def combo(target=defaults['target'], date=defaults['date'], band=defaults['band'], chip=defaults['chip'],
          parentdir=defaults['parent'], rot_val=defaults['rot_val'], no_shift=defaults['no_shift'], astromnet=defaults['astromnet'],
          sky_override_path=defaults['sky_override_path'], removal=defaults['removal'], no_get_files=defaults['no_get_files'],
          no_download=defaults['no_download'], no_mflat=defaults['no_mflat'], bulge=defaults['bulge'],survey=defaults['survey'],
          grb_ra=defaults['grb_ra'], grb_dec=defaults['grb_dec'], grb_coordlist=defaults['grb_coordlist'],
          grb_radius=defaults['grb_radius'],auto_mode=defaults['automode'], input_ramp_lists=defaults['ramplist']
          ):
    print('combo date', date)
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
                                                 ' Can use all options (except -parallel & file download) from photometrus '
                                                 'pipeline & photometrus photometry')
    parser.add_argument('-target', type=str, help='[str] target field, objname in log, ex. "field1234"',
                        default=defaults['target'])
    parser.add_argument('-date', type=str, help='[str] date of observation, in yyyymmdd format',
                        default=defaults['date'])
    parser.add_argument('-band', type=str, help='[str] filter, ex. "J"', default=defaults['band'])
    parser.add_argument('-chip', type=str,
                        help='[str] Optional, use to process specific chips, use "1,2,3,4"'
                             ' format.', default=defaults['chip'])
    parser.add_argument('-no_get_files', action='store_true', help='optional flag, use if you want to download '
                                                                     'through old method (scp), new method passes file '
                                                                     'paths (new saves space and time)',
                        default=defaults['no_get_files'])
    parser.add_argument('-no_download', action='store_true', help='optional flag, use if you *ALREADY* have the data'
                                                                  'downloaded, *NOT* to use new file path method',
                        default=defaults['no_download'])
    parser.add_argument('-no_shift', action='store_true', help='optional flag, DO NOT use astrometric shift'
                                                               ' script in place of astrom.net, will not use either (shift is default)',
                        default=defaults['no_shift'])
    parser.add_argument('-astromnet', action='store_true',
                        help='optional flag, use astrom.net to reinforce astrometry',
                        default=defaults['astromnet'])
    parser.add_argument('-removal', action='store_true',
                        help='optional flag, used to remove intermediate subdirectories and data, leaving only the '
                             'stacks & skies; intended for space saving in large nights of observation')
    parser.add_argument('-parent', type=str, help='[str] *NOW OPTIONAL* specify parent directory to '
                                                  'store all data products, otherwise it will automatically generate w/'
                                                  'the format "/target_date/band/"', default=defaults['parent'])
    parser.add_argument('-rot_val', type=float, help='[float] optional, put in your rot angle in deg,'
                                                     ' if you had a non-default rotation angle in your obs'
                                                     ' (default = 48 deg or 172800")', default=defaults['rot_val'])
    parser.add_argument('-sky_override', type=str, help='[str], Optional path to specify already generated '
                                                        'sky to use in sky sub, skipping sky gen. Input full file path.',
                        default=defaults['sky_override_path'])
    parser.add_argument('-no_mflat', action='store_true', help='optional flag, use if you *DO NOT* want to'
                                                               ' automatically generate mflats for this night if none'
                                                               ' exist', default=defaults['no_mflat'])
    parser.add_argument('-bulge', action='store_true',
                        help='optional flag, utilize setup specifically designed for bulge fields.  Hopefully we can'
                             ' automate this in the future', default=defaults['bulge'])
    parser.add_argument('-survey', type=str, help='Specify specific survey to query for photometry (default'
                                                  ' picks for you), see photometrus single_photometry -h for list of '
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
    args, unknown = parser.parse_known_args()  # TODO: get these default arguments from defaults dict
    print('main', args.date)
    print('args', args)
    print('unknown', unknown)
    combo(args.target, args.date, args.band, args.chip, args.parent, args.rot_val, args.no_shift, args.astromnet,
          args.sky_override, args.removal, args.no_get_files, args.no_download, args.no_mflat, args.bulge, args.survey,
          args.grb_ra, args.grb_dec, args.grb_coordlist, args.grb_radius, defaults['automode'], defaults['ramplist'])


if __name__ == "__main__":
    main()

