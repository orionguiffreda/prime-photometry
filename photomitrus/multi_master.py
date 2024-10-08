import os
import sys
import subprocess
import argparse
from subprocess import Popen

#%%
def datadownload(parentdir,target,filter,date,chip=None):
    print('Creating parent dir: %s' % parentdir)
    if not os.path.isdir(parentdir):
        os.mkdir(parentdir)
    else:
        pass
    if chip:
        if filter == 'Z':
            try:
                command = ('python ./getdata.py -d %s -i 192.168.212.22 -f ramp --objname %s '
                           '-c %s --filter1 %s --filter2 Open %s') % (parentdir, target, chip, filter, date)
                print('Executing command: %s' % command)
                rval = subprocess.run(command.split(), check=True)
            except subprocess.CalledProcessError as err:
                print('Could not run with exit error %s' % err)
        else:
            try:
                command = ('python ./getdata.py -d %s -i 192.168.212.22 -f ramp --objname %s '
                           '-c %s --filter1 Open --filter2 %s %s') % (parentdir,target,chip,filter,date)
                print('Executing command: %s' % command)
                rval = subprocess.run(command.split(), check=True)
            except subprocess.CalledProcessError as err:
                print('Could not run with exit error %s' % err)
    else:
        try:
            command = ('python ./getdata.py -d %s -i 192.168.212.22 -f ramp --objname %s '
                       '-c 1,2,3,4 --filter1 Open --filter2 %s %s') % (parentdir,target,filter,date)
            print('Executing command: %s' % command)
            rval = subprocess.run(command.split(), check=True)
        except subprocess.CalledProcessError as err:
            print('Could not run with exit error %s' % err)

#%%


def refineprocess(parentdir,chip,filter,rot_val=None):
    if rot_val:
        if chip == 1 or chip == 2:
            try:
                command = ('python ./master.py -FF -angle -refine -parent %s '
                           '-chip %i -filter %s -sigma 4 -rot_val %s') % (parentdir,chip,filter,rot_val)
                print('Executing command: %s' % command)
                rval = subprocess.run(command.split(), check=True)
            except subprocess.CalledProcessError as err:
                print('Could not run with exit error %s' % err)
        elif chip == 3 or chip == 4:
            try:
                command = ('python ./master.py -FF -angle -refine -parent %s '
                           '-chip %i -filter %s -sigma 6 -rot_val %s') % (parentdir,chip,filter,rot_val)
                print('Executing command: %s' % command)
                rval = subprocess.run(command.split(), check=True)
            except subprocess.CalledProcessError as err:
                print('Could not run with exit error %s' % err)
    else:
        if chip == 1 or chip == 2:
            try:
                command = ('python ./master.py -FF -angle -refine -parent %s '
                           '-chip %i -filter %s -sigma 4') % (parentdir, chip, filter)
                print('Executing command: %s' % command)
                rval = subprocess.run(command.split(), check=True)
            except subprocess.CalledProcessError as err:
                print('Could not run with exit error %s' % err)
        elif chip == 3 or chip == 4:
            try:
                command = ('python ./master.py -FF -angle -refine -parent %s '
                           '-chip %i -filter %s -sigma 6') % (parentdir, chip, filter)
                print('Executing command: %s' % command)
                rval = subprocess.run(command.split(), check=True)
            except subprocess.CalledProcessError as err:
                print('Could not run with exit error %s' % err)

def shiftprocess(parentdir,chip,filter,rot_val=None):
    if rot_val:
        if chip == 1 or chip == 2:
            try:
                command = ('python ./master.py -FF -angle -shift -parent %s '
                           '-chip %i -filter %s -sigma 4 -rot_val %s') % (parentdir,chip,filter,rot_val)
                print('Executing command: %s' % command)
                rval = subprocess.run(command.split(), check=True)
            except subprocess.CalledProcessError as err:
                print('Could not run with exit error %s' % err)
        elif chip == 3 or chip == 4:
            try:
                command = ('python ./master.py -FF -angle -shift -parent %s '
                           '-chip %i -filter %s -sigma 6 -rot_val %s') % (parentdir,chip,filter,rot_val)
                print('Executing command: %s' % command)
                rval = subprocess.run(command.split(), check=True)
            except subprocess.CalledProcessError as err:
                print('Could not run with exit error %s' % err)
    else:
        if chip == 1 or chip == 2:
            try:
                command = ('python ./master.py -FF -angle -shift -parent %s '
                           '-chip %i -filter %s -sigma 4') % (parentdir,chip,filter)
                print('Executing command: %s' % command)
                rval = subprocess.run(command.split(), check=True)
            except subprocess.CalledProcessError as err:
                print('Could not run with exit error %s' % err)
        elif chip == 3 or chip == 4:
            try:
                command = ('python ./master.py -FF -angle -shift -parent %s '
                           '-chip %i -filter %s -sigma 6') % (parentdir,chip,filter)
                print('Executing command: %s' % command)
                rval = subprocess.run(command.split(), check=True)
            except subprocess.CalledProcessError as err:
                print('Could not run with exit error %s' % err)

def baseprocess(parentdir,chip,filter,rot_val=None):
    if rot_val:
        if chip == 1 or chip == 2:
            try:
                command = ('python ./master.py -FF -angle -parent %s '
                           '-chip %i -filter %s -sigma 4 -rot_val %s') % (parentdir,chip,filter,rot_val)
                print('Executing command: %s' % command)
                rval = subprocess.run(command.split(), check=True)
            except subprocess.CalledProcessError as err:
                print('Could not run with exit error %s' % err)
        elif chip == 3 or chip == 4:
            try:
                command = ('python ./master.py -FF -angle -parent %s '
                           '-chip %i -filter %s -sigma 6 -rot_val %s') % (parentdir,chip,filter,rot_val)
                print('Executing command: %s' % command)
                rval = subprocess.run(command.split(), check=True)
            except subprocess.CalledProcessError as err:
                print('Could not run with exit error %s' % err)
    else:
        if chip == 1 or chip == 2:
            try:
                command = ('python ./master.py -FF -angle -parent %s '
                           '-chip %i -filter %s -sigma 4') % (parentdir, chip, filter)
                print('Executing command: %s' % command)
                rval = subprocess.run(command.split(), check=True)
            except subprocess.CalledProcessError as err:
                print('Could not run with exit error %s' % err)
        elif chip == 3 or chip == 4:
            try:
                command = ('python ./master.py -FF -angle -parent %s '
                           '-chip %i -filter %s -sigma 6') % (parentdir, chip, filter)
                print('Executing command: %s' % command)
                rval = subprocess.run(command.split(), check=True)
            except subprocess.CalledProcessError as err:
                print('Could not run with exit error %s' % err)

#%%

def refineprocessparallel(parentdir,chips,filter):
    commands = []
    for f in chips:
        command = ('python ./master.py -fpack -refine -angle -FF -sex -parent %s '
                   '-chip %i -filter %s') % (parentdir,f,filter)
        commands.append(command)
    print('Processing chips in parallel...')
    procs = [Popen(i.split(),stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) for i in commands]
    for p in procs:
        p.wait()

#%%
def checkprocess(parentdir,chip,filter):
    try:
        command = ('python ./master.py -FF -skygen_start -fpack -parent %s '
                   '-chip %i -filter %s') % (parentdir, chip, filter)
        print('Executing command: %s' % command)
        rval = subprocess.run(command.split(), check=True)
    except subprocess.CalledProcessError as err:
        print('Could not run with exit error %s' % err)

def checkprocessparallel(parentdir,chips,filter):
    commands = []
    for f in chips:
        command = ('python ./master.py -fpack -angle -FF -sex -parent %s '
                   '-chip %i -filter %s') % (parentdir,f,filter)
        commands.append(command)
    print('Processing chips in parallel...')
    procs = [Popen(i.split(),stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) for i in commands]
    for p in procs:
        p.wait()

#%%
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Use to process whole observations (all chips)')
    parser.add_argument('-parallel', action='store_true', help='optional flag, process multiple chips simultaneously')
    parser.add_argument('-no_download', action='store_true', help='optional flag, use if you already have the data')
    parser.add_argument('-shift', action='store_true', help='optional flag, use astrometric shift script in place of astrom.net')
    parser.add_argument('-astromnet', action='store_true', help='optional flag, use astrom.net to reinforce astrometry')
    parser.add_argument('-parent', type=str, help='[str] parent directory to store all data products')
    parser.add_argument('-target', type=str, help='[str] target field, objname in log, ex. "field1234"')
    parser.add_argument('-date', type=str, help='[str] date of observation, in yyyymmdd format')
    parser.add_argument('-filter', type=str, help='[str] filter, ex. "J"')
    parser.add_argument('-chip', type=int, help='[int] Optional, use to process only 1 specific chip',default=None)
    parser.add_argument('-rot_val', type=float, help='[float] optional, put in your rot angle in deg,'
                                                     ' if you had a non-default rotation angle in your obs'
                                                     ' (default = 48 deg or 172800")', default=None)
    args = parser.parse_args()

    if args.no_download:
        pass
    else:
        datadownload(args.parent,args.target,args.filter,args.date,args.chip)

    if args.chip:
        chip_path = args.parent+'C%s/' % args.chip
    else:
        chip_path = args.parent + 'C1/'

    if not os.listdir(chip_path):
        print('Error downloading! perhaps wrong date or target?')
    else:
        if not args.chip:
            chips = [1, 2, 3, 4]
        else:
            chips = [args.chip]
        #refinechips = [f for f in chips if f == 1 or f == 3]
        #checkchips = [f for f in chips if f == 2 or f == 4]

        #if args.parallel:
        #    refineprocessparallel(args.parent,refinechips,args.filter)
        #    for f in checkchips:
        #        checkprocess(args.parent, f, args.filter)

        for f in chips:
            if args.astromnet:
                refineprocess(args.parent,f,args.filter,args.rot_val)
            elif args.shift:
                shiftprocess(args.parent,f,args.filter,args.rot_val)
            else:
                baseprocess(args.parent,f,args.filter,args.rot_val)
            #for f in checkchips:
           #     checkprocess(args.parent,f,args.filter)