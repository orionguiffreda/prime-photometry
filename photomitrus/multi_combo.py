from photomitrus import multi_master
from photomitrus.photometry import multi_photom
from photomitrus.settings import PIPELINE_DEFAULT_DIR

import os
import shutil
import argparse

#%%


def combo(target, date, band, chip=None, parentdir=False, rot_val=48, no_shift=False, astromnet=False,
          sky_override_path=False, removal=False, survey=None):

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

    full_ramp_list, m_list = multi_master.datalistdownload(chosen_parent, target, band, date)
    if all(not lst for lst in full_ramp_list):
        print('Error finding files, No data! Or perhaps wrong date or target?')
        if parentdir:
            pass
        else:
            print('Removing default directory: %s' % field_dir)
            try:
                os.chdir(PIPELINE_DEFAULT_DIR)
                shutil.rmtree(field_dir, ignore_errors=True)
            except FileNotFoundError:
                print('Directory already no longer exists.')
    for f in chips:
        if astromnet:
            multi_master.refineprocess(chosen_parent, f, band, rot_val, sky_override_path, removal=removal,
                          fullramplist=full_ramp_list)
            multi_photom.mastermultiphotom(stackpath, band, f, survey)
        elif not no_shift:
            multi_master.shiftprocess(chosen_parent, f, band, rot_val, sky_override_path, removal=removal,
                         fullramplist=full_ramp_list)
            multi_photom.mastermultiphotom(stackpath, band, f, survey)
        else:
            multi_master.baseprocess(chosen_parent, f, band, rot_val, sky_override_path, removal=removal,
                        fullramplist=full_ramp_list)
            multi_photom.mastermultiphotom(stackpath, band, f, survey)

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
    parser.add_argument('-survey', type=str, help='Specify specific survey to query for photometry (default'
                                                  ' picks for you), see photometrus photometry single -h for list of '
                                                  'available surveys')
    args, unknown = parser.parse_known_args()

    combo(args.target, args.date, args.band, args.chip, args.parent, args.rot_val, args.no_shift, args.astromnet,
          args.sky_override, args.removal, args.survey)


if __name__ == "__main__":
    main()

