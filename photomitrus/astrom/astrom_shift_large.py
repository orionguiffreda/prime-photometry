import subprocess
import os
import argparse
from settings import gen_pipeline_file_name

#%%


def split_shift(subpath, num, filter):
    os.chdir(gen_pipeline_file_name())
    print('Splitting obs list for increased accuracy...')
    all_fits = [f for f in sorted(os.listdir(subpath)) if f.endswith('.flat.fits')]
    tot = len(all_fits)
    for i in range(0, tot, num):
        split_arr = all_fits[0 + i:num + i]
        print('\nApplying shift to %i files: %s - %s\n' % (len(split_arr), split_arr[0], split_arr[-1]))
        imgname = split_arr[0]
        split_arr_arg = ' '.join(split_arr)
        try:
            command = 'python ./astrom/astrom_shift.py -remove -segment -pipeline -split -dir %s -imagename %s -filter %s -arr %s' % (subpath, imgname, filter, split_arr_arg)
            print('Executing command: %s' % command)
            rval = subprocess.run(command.split(), check=True)
        except subprocess.CalledProcessError as err:
            print('Could not run with exit error %s' % err)

def double_run(subpath, filter):
    os.chdir(gen_pipeline_file_name())
    try:
        command = 'python ./astrom/astrom_shift.py -remove -segment -pipeline -dir %s -imagename %s -filter %s' % (subpath, imgname, filter)
        print('Executing command: %s' % command)
        rval = subprocess.run(command.split(), check=True)
    except subprocess.CalledProcessError as err:
        print('Could not run with exit error %s' % err)

#%%

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='For observations with many files, either gives option to run astrom_shift '
                                                 'again, or splits up the observation file list'
                                                 ' into sublists of specified length, runs astrom_shift on sublists')
    parser.add_argument('-dir', type=str, help='[str] path where input files are stored (should run on proc. image, '
                                               'so likely should be /C#_sub/)')
    parser.add_argument('-filter', type=str, help='[str] filter used, ex. "J"')
    parser.add_argument('-num', type=int, help='[int] # of files to split list by, ex. -num 50 would run on file '
                                               'index 0 and apply to file indices 0-49, then run on file index 50 and apply to '
                                               'indices 50-99, etc')

    args = parser.parse_args()

    split_shift(args.dir, args.num, args.filter)
