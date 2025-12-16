"""
Runs processing on whole observations
"""

import os
import sys
import subprocess
import argparse
from subprocess import Popen
import shutil

from photometrus.preprocess import auto_flat
from photometrus.getdata import download_data
from photometrus.getfiles import get_data_files
from photometrus.master import master
from photometrus.preprocess.auto_flat import gen_mflat_file_name_new
from photometrus.utils.defaults import PROCESSING_DEFAULTS as defaults

#%%


def auto_mflat_gen(date, band, chip):
    file = gen_mflat_file_name_new(band, chip, date, sflat=True)

    # print("auto_mflat_gen date", date)
    # check = mflat_checker(date)
    # print('Running auto-MFLAT generation...')
    # if check is False:
    #     auto_flat.autoflatgen(date)
    # else:
    #     pass


def parentcreation(target, date, band):
    field_dir_name = '%s_%s' % (target, date)
    field_dir = os.path.join(defaults['parent'], field_dir_name)
    parent_dir = os.path.join(field_dir, band)
    print('\nDefault parent dir: %s' % parent_dir)
    if not os.path.isdir(parent_dir):
        print('Generating parent directory!\n')
        os.makedirs(parent_dir)
    else:
        pass
    return parent_dir, field_dir

#%%


def datalistdownload(parentdir, target, band, date, header_filter=None):
    print('Verifying parent dir: %s' % parentdir)
    if not os.path.isdir(parentdir):
        os.makedirs(parentdir)
    else:
        pass
    if band == 'Z':
        try:
            full_ramp_list, m_list = get_data_files(date=date, objname=target, filter1=band, filter2='Open', header_filter=header_filter)
            if any(lst for lst in m_list):
                print('Missing files!: ', m_list)
        except:
            print('Error fetching data!')
            sys.exit(0)
    else:
        full_ramp_list, m_list = get_data_files(date=date, objname=target, filter1='Open', filter2=band, header_filter=header_filter)
        if any(lst for lst in m_list):
            print('Missing files!: ', m_list)
    return full_ramp_list, m_list


def datadownload(parentdir, target, band, date, chip, header_filter=None):
    print('Verifying parent dir: %s' % parentdir)
    if not os.path.isdir(parentdir):
        os.mkdir(parentdir)
    else:
        pass
    if band == 'Z':
        try:
            download_data(date=date, ip='192.168.212.22', save_dir=parentdir, objname=target, chip=chip,
                          filter1=band, filter2='Open')
        except:
            print('Error downloading data!')
    else:
        try:
            download_data(date=date, ip='192.168.212.22', save_dir=parentdir, objname=target, chip=chip,
                          filter1='Open', filter2=band)
        except:
            print('Error downloading data!')

#%%


def refineprocess(parentdir, chip, band, date, rot_val, sky_override_path, removal, fullramplist,
                  bulge, rampnum):
    if chip == 1 or chip == 2:
        sigma = defaults['sigma']
    elif chip == 3 or chip == 4:
        sigma = defaults['sigma'] + 2
    else:
        sigma = None
    master(parentdir=parentdir, chip=chip, band=band, sigma=sigma, date=date, rot_val=rot_val,
           sky_override=sky_override_path, removal=removal, fullramplist=fullramplist, bulge=bulge, rampnum=rampnum)


def shiftprocess(parentdir, chip, band, date, rot_val, sky_override_path, removal, fullramplist,
                 bulge, rampnum):
    if chip == 1 or chip == 2:
        sigma = defaults['sigma']
    elif chip == 3 or chip == 4:
        sigma = defaults['sigma'] + 2
    else:
        sigma = None
    master(parentdir=parentdir, chip=chip, band=band, sigma=sigma, date=date, rot_val=rot_val, sky_override=sky_override_path,
           removal=removal, fullramplist=fullramplist, bulge=bulge, rampnum=rampnum)


def baseprocess(parentdir, chip, band, date, rot_val, sky_override_path, removal, fullramplist,
                bulge, rampnum):
    if chip == 1 or chip == 2:
        sigma = defaults['sigma']
    elif chip == 3 or chip == 4:
        sigma = defaults['sigma'] + 2
    else:
        sigma = None
    master(parentdir=parentdir, chip=chip, band=band, sigma=sigma, date=date, rot_val=rot_val,
           sky_override=sky_override_path, removal=removal, fullramplist=fullramplist, bulge=bulge, rampnum=rampnum)


def processparallel(target, date, band, chips):
    commands = []
    for chip in chips:
        if chip == 1 or chip == 2:
            sigma = 4
        elif chip == 3 or chip == 4:
            sigma = 6
        # command = master(parentdir=parentdir, chip=chip, band=band, sigma=sigma, rot_val=rot_val, sky_override=sky_override_path)
        command = ('photometrus pipeline -target %s -date %s -band %s -sigma %s -no_download' % (target, date, band, sigma))
        commands.append(command)
    print(commands)
    procs = [Popen(i.split()) for i in commands]
    for p in procs:
        p.wait()

#%%


def multi_master(
        target=defaults['target'], date=defaults['date'], band=defaults['band'], chip=defaults['chip'],
        parentdir=defaults['parent'], rot_val=defaults['rot_val'], no_shift=defaults['no_shift'], astromnet=defaults['astromnet'],
        no_download=defaults['no_download'], sky_override_path=defaults['sky_override_path'], removal=defaults['removal'],
        no_get_files=defaults['no_get_files'], no_mflat=defaults['no_mflat'], rampnum=defaults['rampnum'],
        bulge=defaults['bulge'], auto_mode=defaults['automode'], input_ramp_lists=defaults['ramplist'], header_filter=None
):
    if not no_mflat:
        auto_mflat_gen(date, band, chip)

    if parentdir != defaults['parent']:
        chosen_parent = parentdir
    else:
        chosen_parent, field_dir = parentcreation(target, date, band)

    if not chip:
        chips = [1, 2, 3, 4]
    else:
        if type(chip) is int:
            chips = chip
        else:
            chips = chip.split(',')
            chips = [int(f) for f in chips]

    if no_get_files:
        print('Omitting usage of getfiles!')
        if not no_download:
            chips_str = ','.join(str(x) for x in [chips])
            datadownload(chosen_parent, target, band, date, chips_str)

        if chip:
            if type(chip) is int:
                chip_path = os.path.join(chosen_parent, 'C%s/' % chips)
            else:
                chip_path = os.path.join(chosen_parent, 'C%s/' % chips[0])
            print(chip_path)
        else:
            chip_path = os.path.join(chosen_parent, 'C1/')

        if not os.listdir(chip_path):
            print('Error downloading, No data! Or perhaps wrong date or target?')
            if parentdir:
                pass
            else:
                print('Removing default directory: %s' % field_dir)
                try:
                    os.chdir(defaults['parent'])
                    shutil.rmtree(field_dir, ignore_errors=True)
                except FileNotFoundError:
                    print('Directory already no longer exists.')
        if type(chips) is int:
            if astromnet:
                refineprocess(chosen_parent, chips, band, date, rot_val, sky_override_path, removal=removal, bulge=bulge,
                              fullramplist=None, rampnum=rampnum)
            elif not no_shift:
                shiftprocess(chosen_parent, chips, band, date, rot_val, sky_override_path, removal=removal, bulge=bulge,
                              fullramplist=None, rampnum=rampnum)
            else:
                baseprocess(chosen_parent, chips, band, date, rot_val, sky_override_path, removal=removal, bulge=bulge,
                              fullramplist=None, rampnum=rampnum)
        else:
            for f in chips:
                if astromnet:
                    refineprocess(chosen_parent, f, band, date, rot_val, sky_override_path, removal=removal, bulge=bulge,
                              fullramplist=None, rampnum=rampnum)
                elif not no_shift:
                    shiftprocess(chosen_parent, f, band, date, rot_val, sky_override_path, removal=removal, bulge=bulge,
                              fullramplist=None, rampnum=rampnum)
                else:
                    baseprocess(chosen_parent, f, band, date, rot_val, sky_override_path, removal=removal, bulge=bulge,
                              fullramplist=None, rampnum=rampnum)

    else:
        if not auto_mode:
            full_ramp_list, m_list = datalistdownload(chosen_parent, target, band, date, header_filter=header_filter)
        else:
            full_ramp_list = input_ramp_lists[0]
            m_list = input_ramp_lists[1]
            # removal = False
        if all(not lst for lst in full_ramp_list):
            print('Error finding files, No data! Or perhaps wrong date or target?')
            if parentdir:
                pass
            else:
                print('Removing default directory: %s' % field_dir)
                try:
                    os.chdir(defaults['parent'])
                    shutil.rmtree(field_dir, ignore_errors=True)
                except FileNotFoundError:
                    print('Directory already no longer exists.')
        # if parallel:
        #     print('Processing all chosen chips in parallel!')
        #     processparallel(target, date, band, chips)
        if type(chips) is int:
            if astromnet:
                refineprocess(chosen_parent, chips, band, date, rot_val, sky_override_path, removal=removal, bulge=bulge,
                              fullramplist=full_ramp_list, rampnum=rampnum)
            elif not no_shift:
                shiftprocess(chosen_parent, chips, band, date, rot_val, sky_override_path, removal=removal, bulge=bulge,
                             fullramplist=full_ramp_list, rampnum=rampnum)
            else:
                baseprocess(chosen_parent, chips, band, date, rot_val, sky_override_path, removal=removal, bulge=bulge,
                            fullramplist=full_ramp_list, rampnum=rampnum)
        else:
            for f in chips:
                if astromnet:
                    refineprocess(chosen_parent, f, band, date, rot_val, sky_override_path, removal=removal, bulge=bulge,
                                  fullramplist=full_ramp_list, rampnum=rampnum)
                elif not no_shift:
                    shiftprocess(chosen_parent, f, band, date, rot_val, sky_override_path, removal=removal, bulge=bulge,
                                 fullramplist=full_ramp_list, rampnum=rampnum)
                else:
                    baseprocess(chosen_parent, f, band, date, rot_val, sky_override_path, removal=removal, bulge=bulge,
                                fullramplist=full_ramp_list, rampnum=rampnum)


def main():
    parser = argparse.ArgumentParser(description='Use to process whole observations (all chips)')
    # parser.add_argument('-parallel', action='store_true', help='optional flag, process multiple chips simultaneously,'
    #                                                            ' only use on obs. w/ small amount of images!')
    parser.add_argument('-parent', type=str, help='[str] *NOW OPTIONAL* specify parent directory to '
                                                  'store all data products, otherwise it will automatically generate @ '
                                                  'the default directory w/ the format "/target_date/band/"',
                                            default=defaults['parent'])
    parser.add_argument('-target', type=str, help='[str] target field, objname in log, ex. "field1234"',
                        default=defaults['target'])
    parser.add_argument('-date', type=str, help='[str] date of observation, in yyyymmdd format',
                        default=defaults['date'])
    parser.add_argument('-band', type=str, help='[str] filter, ex. "J"', default=defaults['band'])
    parser.add_argument('-chip', type=str, help='[str] Optional, use to process specific chips, use "1,2,3,4"'
                                                ' format.', default=defaults['chip'])
    parser.add_argument('-header_filter', type=str, default=None,
        help='[str] Optional, use to require fits headers match the desired pairing. Keys and headers are ":" delimited and key:header pairs are "," delimited'
    )
    parser.add_argument('-no_get_files', action='store_true', help='optional flag, use if you want to download '
                                                                     'through old method (scp), new method passes file '
                                                                     'paths (new saves space and time)',
                        default=defaults['no_get_files'])
    parser.add_argument('-no_download', action='store_true', help='optional flag, use if you *ALREADY* have the data'
                                                                  'downloaded, *NOT* to use new file path method.  Use w/ '
                                                                  '-no_get_files!',
                        default=defaults['no_download'])
    parser.add_argument('-no_shift', action='store_true', help='optional flag, DO NOT use astrometric shift'
                                            ' script in place of astrom.net, will not use either (shift is default)',
                        default=defaults['no_shift'])
    parser.add_argument('-astromnet', action='store_true', help='optional flag, use astrom.net to reinforce astrometry',
                        default=defaults['astromnet'])
    parser.add_argument('-removal', action='store_true',
                        help='optional flag, used to remove intermediate subdirectories and data, leaving only the '
                             'stacks & skies; intended for space saving in large nights of observation',
                        default=defaults['removal'])
    parser.add_argument('-rot_val', type=float, help='[float] optional, put in your rot angle in deg,'
                                                     ' if you had a non-default rotation angle in your obs'
                                                     ' (default = 48 deg or 172800")', default=defaults['rot_val'])
    parser.add_argument('-sky_override', type=str, help='[str], Optional path to specify already generated '
                                                        'sky to use in sky sub, skipping sky gen. Input full file path.',
                        default=defaults['sky_override_path'])
    parser.add_argument('-no_mflat', action='store_true', help='optional flag, use if you *DO NOT* want to'
                                                               ' automatically generate mflats for this night if none'
                                                               ' exist', default=defaults['no_mflat'])
    parser.add_argument('-rampnum', type=int, help='[int], optional arg to specify how many ramp images from '
                                                   'your observation you want to include in the processing, helpful for '
                                                   'dodging bad individual images',
                        default=defaults["rampnum"])
    parser.add_argument('-bulge', action='store_true',
                        help='optional flag, utilize setup specifically designed for bulge fields.  Hopefully we can'
                             ' automate this in the future', default=defaults['bulge'])
    parser.add_argument('-ramplist', type=str,
                        help='[str], list of filepaths to ramp files to conduct processing on, in format: '
                             '"image1.fits,image2.fits,image3.fits,...".  For usage w/ '
                             'multi_master.py, better alternative to downloading ramp data', default=defaults['ramplist'])
    parser.add_argument('-auto', action='store_true',
                        help='use setup designed for full night automatic processing', default=defaults['automode'])
    args, unknown = parser.parse_known_args()  # TODO: get default arguments from defaults dict
    print('unknown args', unknown)
    multi_master(target=args.target, date=args.date, band=args.band, chip=args.chip, parentdir=args.parent,
                 rot_val=args.rot_val, no_shift=args.no_shift, astromnet=args.astromnet, no_download=args.no_download,
                 sky_override_path=args.sky_override, removal=args.removal, no_get_files=args.no_get_files,
                 no_mflat=args.no_mflat, rampnum=args.rampnum, bulge=args.bulge, auto_mode=args.auto, input_ramp_lists=args.ramplist,
                 header_filter=args.header_filter)


if __name__ == "__main__":
    main()
