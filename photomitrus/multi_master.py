import os
import sys
import subprocess
import argparse
from subprocess import Popen
import shutil

from photomitrus.getdata import download_data
from photomitrus.master import master
from photomitrus.settings import PIPELINE_DEFAULT_DIR

#%%


def parentcreation(target, date, band):
    field_dir_name = '%s_%s' % (target, date)
    field_dir = os.path.join(PIPELINE_DEFAULT_DIR, field_dir_name)
    parent_dir = os.path.join(field_dir, band)
    print('Default parent dir: %s' % parent_dir)
    if not os.path.isdir(parent_dir):
        print('Generating parent directory!\n')
        os.makedirs(parent_dir)
    else:
        pass
    return parent_dir, field_dir

#%%


def datadownload(parentdir, target, band, date, chip):
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


def refineprocess(parentdir, chip, band, rot_val=48, sky_override_path=None, removal=False):
    if chip == 1 or chip == 2:
        sigma = 4
    elif chip == 3 or chip == 4:
        sigma = 6
    else:
        sigma = None
    master(parentdir=parentdir, chip=chip, band=band, sigma=sigma, rot_val=rot_val,net_refine=True,
           sky_override=sky_override_path, removal=removal)


def shiftprocess(parentdir, chip, band, rot_val=48, sky_override_path=None, removal=False):
    if chip == 1 or chip == 2:
        sigma = 4
    elif chip == 3 or chip == 4:
        sigma = 6
    else:
        sigma = None
    master(parentdir=parentdir, chip=chip, band=band, sigma=sigma, rot_val=rot_val, sky_override=sky_override_path,
           removal=removal)


def baseprocess(parentdir, chip, band, rot_val=48, sky_override_path=None, removal=False):
    if chip == 1 or chip == 2:
        sigma = 4
    elif chip == 3 or chip == 4:
        sigma = 6
    else:
        sigma = None
    master(parentdir=parentdir, chip=chip, band=band, sigma=sigma, rot_val=rot_val, no_shift=True,
           sky_override=sky_override_path, removal=removal)


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
        target, date, band, chip=None, parentdir=None, rot_val=48, no_shift=False, astromnet=False, parallel=False,
        no_download=False, sky_override_path=None, removal=False
):

    if parentdir:
        chosen_parent = parentdir
    else:
        chosen_parent, field_dir = parentcreation(target, date, band)

    if not chip:
        chips = [1, 2, 3, 4]
    else:
        chips = chip.split(',')
        chips = [int(f) for f in chips]

    if no_download:
        pass
    else:
        chips_str = ','.join(str(x) for x in chips)
        datadownload(chosen_parent, target, band, date, chips_str)

    if chip:
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
                os.chdir(PIPELINE_DEFAULT_DIR)
                shutil.rmtree(field_dir, ignore_errors=True)
            except FileNotFoundError:
                print('Directory already no longer exists.')
    elif parallel:
        print('Processing all chosen chips in parallel!')
        processparallel(target, date, band, chips)
    else:
        for f in chips:
            if astromnet:
                refineprocess(chosen_parent, f, band, rot_val, sky_override_path, removal=removal)
            elif not no_shift:
                shiftprocess(chosen_parent, f, band, rot_val, sky_override_path, removal=removal)
            else:
                baseprocess(chosen_parent, f, band, rot_val, sky_override_path, removal=removal)


def main():
    parser = argparse.ArgumentParser(description='Use to process whole observations (all chips)')
    parser.add_argument('-parallel', action='store_true', help='optional flag, process multiple chips simultaneously,'
                                                               ' only use on obs. w/ small amount of images!')
    parser.add_argument('-no_download', action='store_true', help='optional flag, use if you already have the data')
    parser.add_argument('-no_shift', action='store_true', help='optional flag, DO NOT use astrometric shift'
                                                               ' script in place of astrom.net, will not use either (shift is default)')
    parser.add_argument('-astromnet', action='store_true', help='optional flag, use astrom.net to reinforce astrometry')
    parser.add_argument('-removal', action='store_true',
                        help='optional flag, used to remove intermediate subdirectories and data, leaving only the '
                             'stacks & skies; intended for space saving in large nights of observation')
    parser.add_argument('-parent', type=str, help='[str] *NOW OPTIONAL* specify parent directory to '
                                                  'store all data products, otherwise it will automatically generate w/'
                                                  'the format "/target_date/band/"', default=None)
    parser.add_argument('-target', type=str, help='[str] target field, objname in log, ex. "field1234"')
    parser.add_argument('-date', type=str, help='[str] date of observation, in yyyymmdd format')
    parser.add_argument('-band', type=str, help='[str] filter, ex. "J"')
    parser.add_argument('-chip', type=str, help='[str] Optional, use to process specific chips',default=None)
    parser.add_argument('-rot_val', type=float, help='[float] optional, put in your rot angle in deg,'
                                                     ' if you had a non-default rotation angle in your obs'
                                                     ' (default = 48 deg or 172800")', default=48)
    parser.add_argument('-sky_override', type=str, help='[str], Optional path to specify already generated '
                                                        'sky to use in sky sub, skipping sky gen. Input full file path.',
                        default=None)
    args, unknown = parser.parse_known_args()

    multi_master(args.target, args.date, args.band, args.chip, args.parent, args.rot_val, args.no_shift, args.astromnet,
                 args.parallel, args.no_download, args.sky_override, args.removal)


if __name__ == "__main__":
    main()
