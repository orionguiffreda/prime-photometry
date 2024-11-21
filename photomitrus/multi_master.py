import os
import sys
import subprocess
import argparse
from subprocess import Popen

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
    return parent_dir

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


def refineprocess(parentdir, chip, band, rot_val=48, sky_override_path=None):
    if chip == 1 or chip == 2:
        sigma = 4
    elif chip == 3 or chip == 4:
        sigma = 6
    else:
        sigma = None
    master(parentdir=parentdir, chip=chip, band=band, sigma=sigma, rot_val=rot_val, net_refine=True,
           sky_override=sky_override_path)


def shiftprocess(parentdir, chip, band, rot_val=48, sky_override_path=None):
    if chip == 1 or chip == 2:
        sigma = 4
    elif chip == 3 or chip == 4:
        sigma = 6
    else:
        sigma = None
    master(parentdir=parentdir, chip=chip, band=band, sigma=sigma, rot_val=rot_val, sky_override=sky_override_path)


def baseprocess(parentdir, chip, band, rot_val=48, sky_override_path=None):
    if chip == 1 or chip == 2:
        sigma = 4
    elif chip == 3 or chip == 4:
        sigma = 6
    else:
        sigma = None
    master(parentdir=parentdir, chip=chip, band=band, sigma=sigma, rot_val=rot_val, no_shift=True,
           sky_override=sky_override_path)

#%%


def refineprocessparallel(parentdir,chips,band):
    commands = []
    for f in chips:
        command = ('python ./master.py -fpack -refine -angle -FF -sex -parent %s '
                   '-chip %i -filter %s') % (parentdir,f,band)
        commands.append(command)
    print('Processing chips in parallel...')
    procs = [Popen(i.split(),stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) for i in commands]
    for p in procs:
        p.wait()

#%%


def checkprocess(parentdir,chip,band):
    try:
        command = ('python ./master.py -FF -skygen_start -fpack -parent %s '
                   '-chip %i -filter %s') % (parentdir, chip, band)
        print('Executing command: %s' % command)
        subprocess.run(command.split(), check=True)
    except subprocess.CalledProcessError as err:
        print('Could not run with exit error %s' % err)


def checkprocessparallel(parentdir,chips,band):
    commands = []
    for f in chips:
        command = ('python ./master.py -fpack -angle -FF -sex -parent %s '
                   '-chip %i -filter %s') % (parentdir,f,band)
        commands.append(command)
    print('Processing chips in parallel...')
    procs = [Popen(i.split(),stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) for i in commands]
    for p in procs:
        p.wait()

#%%


def multi_master(
        target, date, band, chip, parentdir=None, rot_val=48, no_shift=False, astromnet=False, parallel=False,
        no_download=False, sky_override_path=None
):

    if parentdir:
        chosen_parent = parentdir
    else:
        chosen_parent = parentcreation(target, date, band)

    if not chip:
        chips = [1, 2, 3, 4]
    else:
        chips = [chip]

    if no_download:
        pass
    else:
        chips_str = ','.join(str(x) for x in chips)
        datadownload(chosen_parent, target, band, date, chips_str)

    if chip:
        chip_path = os.path.join(chosen_parent, 'C%s/' % chip)
    else:
        chip_path = os.path.join(chosen_parent, 'C1/')

    if not os.listdir(chip_path):
        print('Error downloading! perhaps wrong date or target?')
    else:

        # if args.parallel:
        #    refineprocessparallel(args.parent,refinechips,args.filter)
        #    for f in checkchips:
        #        checkprocess(args.parent, f, args.filter)

        for f in chips:
            if astromnet:
                refineprocess(chosen_parent, f, band, rot_val, sky_override_path)
            elif not no_shift:
                shiftprocess(chosen_parent, f, band, rot_val, sky_override_path)
            else:
                baseprocess(chosen_parent, f, band, rot_val, sky_override_path)


def main():
    parser = argparse.ArgumentParser(description='Use to process whole observations (all chips)')
    parser.add_argument('-parallel', action='store_true', help='optional flag, process multiple chips simultaneously')
    parser.add_argument('-no_download', action='store_true', help='optional flag, use if you already have the data')
    parser.add_argument('-no_shift', action='store_true', help='optional flag, DO NOT use astrometric shift'
                                                               ' script in place of astrom.net, will not use either (shift is default)')
    parser.add_argument('-astromnet', action='store_true', help='optional flag, use astrom.net to reinforce astrometry')
    parser.add_argument('-parent', type=str, help='[str] *NOW OPTIONAL* specify parent directory to '
                                                  'store all data products, otherwise it will automatically generate w/'
                                                  'the format "/target_date/band/"', default=None)
    parser.add_argument('-target', type=str, help='[str] target field, objname in log, ex. "field1234"')
    parser.add_argument('-date', type=str, help='[str] date of observation, in yyyymmdd format')
    parser.add_argument('-band', type=str, help='[str] filter, ex. "J"')
    parser.add_argument('-chip', type=int, help='[int] Optional, use to process only 1 specific chip',default=None)
    parser.add_argument('-rot_val', type=float, help='[float] optional, put in your rot angle in deg,'
                                                     ' if you had a non-default rotation angle in your obs'
                                                     ' (default = 48 deg or 172800")', default=48)
    parser.add_argument('-sky_override', type=str, help='[str], Optional path to specify already generated '
                                                        'sky to use in sky sub, skipping sky gen. Input full file path.',
                        default=None)
    args, unknown = parser.parse_known_args()

    multi_master(args.target, args.date, args.band, args.chip, args.parent, args.rot_val, args.no_shift, args.astromnet,
                 args.parallel, args.no_download, args.sky_override)


if __name__ == "__main__":
    main()
