#!/usr/bin/python
"""
Author: Joe Durbak, durbak.3@gmail.com

CLI to download PRIME data

installation instructions (assumes Anaconda or Miniconda are installed):
  1) Install package dependencies:
        Option A) new environment
          conda create -n getdata pysftp pandas
          conda activate getdata
        Option B)
          conda install pysftp pandas
  2) Download data:
      python getdata.py <yyyymmdd>
  3) Check out command line options:
      python getdata.py -h
"""
import os

from six.moves.urllib import request
from argparse import ArgumentParser

import pysftp
from pandas import read_csv, to_datetime

from photometrus.utils.defaults import FILE_DEFAULTS

file_types = ['raw_fz', 'ramp', 'raw']
file_types_str = ','.join(file_types)

defaults = FILE_DEFAULTS
filter_1_options = ','.join(['NB', 'Open', 'Z', 'Dark'])
filter_2_options = ','.join(['Open', 'Y', 'J', 'H'])

log_link_format = 'http://www-ir.ess.sci.osaka-u.ac.jp/prime_staff/LOG/Ramp_LOG/ramp_fit_log_{}.dat'
log_folder_location = '/prime01/work/LOG/Ramp_LOG/ramp_fit_log_{}.dat'

cam_dirs = ['home', 'xion2', 'xion3', 'xion4']
cam_dirs_dict_all = {(cam_num+1): cam_dir for cam_num, cam_dir in enumerate(cam_dirs)}
download_directory_format = '/{}/prime/Data/{}'
file_prefix = '{0:08d}'
file_formats = {
    'ramp': 'C{0}.fits.ramp',
    'raw': 'C{0}.fits',
    'raw_fz': 'C{0}.fits.fz'
}

file_type_dirs = {
    'ramp': 'ramp',
    'raw': 'raw',
    'raw_fz': 'raw_fz'
}

raw_fz_subdir = '{1:08d}'

remote_file_formats = {
    'ramp': [
        '/'.join(
            (download_directory_format.format(_dir, 'ramp'), file_prefix+file_formats['ramp'].format(_i+1))
        ) for _i, _dir in enumerate(cam_dirs)
    ],
    'raw': [
        '/'.join(
            (download_directory_format.format(_dir, 'raw'), file_prefix+file_formats['raw'].format(_i+1))
        ) for _i, _dir in enumerate(cam_dirs)
    ],
    'raw_fz': [
        '/'.join(('/xion3/prime/raw_fz/C1', raw_fz_subdir, file_prefix+file_formats['raw_fz'].format(1))),
        '/'.join(('/xion4/prime/raw_fz/C2', raw_fz_subdir, file_prefix+file_formats['raw_fz'].format(2))),
        '/'.join(('/xion3/prime/raw_fz/C3', raw_fz_subdir, file_prefix+file_formats['raw_fz'].format(3))),
        '/'.join(('/xion4/prime/raw_fz/C4', raw_fz_subdir, file_prefix+file_formats['raw_fz'].format(4))),
    ]
}

replace_list = [
    ('all sky grid', 'all_sky_grid'),
    ('no_grid  ', 'no_grid '),
    ('no_grid_test  ', 'no_grid_test '),
    ('field10554_ test', 'field10554_test'),
    ('fucus test', 'focus_test'),
    ('galactic plane', 'galactic_plane'),
    ('standard star', 'standard_star'),
    ('rotator test', 'rotator_test'),
    ('SKY M7', 'SKY_M7'),
    ('LVC S230522n', 'LVC_S230522n'),
    ('Dec test', 'Dec_test'),
    ('Galctic plane', 'Galactic_plane'),
    ('Glactic plane', 'Galactic_plane'),
    ('Galactic plane', 'Galactic_plane'),
    ('Galactic Plane', 'Galactic_plane'),
    ('Vignetting Test', 'Vignetting_Test'),
    ('SPIS J0539-0059', 'SPIS_J0539-0059'),
    ('Frost Check', 'Frost_Check'),
    ('GainTest LED', 'GainTest_LED'),

    # ('OH FilterTest', 'OH_FilterTest'),
    # ('OY FilterTest', 'OY_FilterTest'),
    # ('OJ FilterTest', 'OJ_FilterTest'),
    # ('OO FilterTest', 'OO_FilterTest'),
    # ('ZH FilterTest', 'ZH_FilterTest'),
    # ('ZY FilterTest', 'ZY_FilterTest'),
    # ('ZJ FilterTest', 'ZJ_FilterTest'),
    # ('ZO FilterTest', 'ZO_FilterTest'),
    # ('NBH FilterTest', 'NBH_FilterTest'),
    # ('NBY FilterTest', 'NBY_FilterTest'),
    # ('NBJ FilterTest', 'NBJ_FilterTest'),
    # ('NBO FilterTest', 'NBO_FilterTest'),
    # ('DO FilterTest', 'DO_FilterTest'),
]


def replace_bad_string(lines):
    replaced_lines = []
    for line in lines:
        fixed_line = line
        for replacement in replace_list:
            fixed_line = fixed_line.replace(replacement[0], replacement[1])
        replaced_lines.append(fixed_line)
    return replaced_lines


def input_none(var, var_string):
    if var is None:
        return input(var_string+': ')
    else:
        return var


def get_log_file(
        date, ip=defaults['ip'], user=defaults['user'], password=defaults['password'], overwrite=defaults['redownload'],
        save_dir=defaults['save_dir'], n_retry=defaults['n_retry']
):
    # log_link = log_link_format.format(date)
    ip = input_none(ip, 'ip')
    user = input_none(user, 'user')
    password = input_none(password, 'password')
    log_link = log_folder_location.format(date)
    print(log_link)
    filename = log_link.split('/')[-1]
    clean_filename = filename.replace('.dat', '.clean.dat')
    if not os.path.isfile(clean_filename) or overwrite:
        # output_file, http_message = request.urlretrieve(log_link, filename)
        with pysftp.Connection(ip, username=user, password=password) as sftp:
            exists = sftp.exists(log_link)
            if not exists:
                raise FileNotFoundError('Log file, {}, not found.'.format(log_link))
            get_with_retry(sftp, log_link, filename, overwrite, save_dir, n_retry)
        with open(filename, 'r') as _f:
            read_lines = _f.readlines()
        cleaned_lines = [
            line for line in read_lines if not (
                    line.startswith('Warning') or line.startswith('ERROR') or line.startswith('#') or (line.find('!')!=-1)
            )
        ]
        first_line = ' '.join([field.strip() for field in read_lines[0][1:].split(',')]) + '\n'
        cleaned_lines = replace_bad_string(cleaned_lines)
        cleaned_lines = [first_line] + cleaned_lines
        print(first_line)
        with open(clean_filename, 'w') as _f:
            _f.writelines(cleaned_lines)
    df = read_csv(clean_filename, sep=' ', on_bad_lines='warn')
    return df


def filter_df(
    df, objname=defaults['objname'], objtype=defaults['objtype'], observer=defaults['observer'],
    filter1=defaults['filter1'], filter2=defaults['filter2']
):
    if objname is not None:
        df = df[df['OBJNAME'] == objname]
    if objtype is not None:
        df = df[df['OBJTYPE'] == objtype]
    if observer is not None:
        df = df[df['OBSERVER'] == observer]
    if filter1 is not None:
        df = df[df['FILTER1'] == filter1]
    if filter2 is not None:
        df = df[df['FILTER2'] == filter2]
    return df


def get_file_numbers(df, ftype=defaults['ftype']):
    f_numbers = [int(f_c1[:8]) for f_c1 in df['filename']]
    if 'raw' in ftype:
        new_f_numbers = []
        nframes = df['NFRAMES']
        for start, frames in zip(f_numbers, nframes):
            for frame in range(int(frames)):
                new_f_numbers.append(start+frame)
        f_numbers = new_f_numbers
    return f_numbers


def truncate_1000(number):
    return int(number/1000) * 1000


def get_file_names(file_numbers, camera, ftype=defaults['ftype']):
    remote_file_format = remote_file_formats[ftype][camera]
    if ftype.endswith('fz'):
        file_names = [remote_file_format.format(_n, truncate_1000(_n)) for _n in file_numbers]
    else:
        file_names = [remote_file_format.format(_n) for _n in file_numbers]
    print(file_names)
    return file_names


def get_with_retry(
    sftp, f_name, save_name, overwrite=defaults['overwrite'], save_dir=defaults['save_dir'], n_retry=defaults['n_retry']
):
    retries = 0
    # save_name = f_name.replace('C4', 'C1')
    # save_name = f_name.replace('.fits.ramp', '.ramp.fits')
    if overwrite or not os.path.isfile(save_name):
        while retries < n_retry:
            try:
                sftp.get(f_name, save_name)
                print('success {}'.format(f_name))
                retries = n_retry
            except OSError as e:
                os.remove(save_name)
                print(e)
                retries += 1
        if not os.path.isfile(save_name):
            print('failed to download {}'.format(save_name))
            logfile = os.path.join(save_dir, 'failed_downloads_list.txt')
            with open(logfile, 'a') as _f:
                _f.write(save_name)


def download_camera_data(
    camera, file_numbers, save_dir=defaults['save_dir'], overwrite=defaults['overwrite'], ftype=defaults['ftype'],
    user=defaults['user'], ip=defaults['ip'], password=defaults['password'],
    n_retry=defaults['n_retry']
):
    ip = input_none(ip, 'ip')
    password = input_none(password, 'password')
    user = int(input_none(user, 'user'))
    # file_format = file_formats[ftype]
    # cam_dir = cam_dirs[camera]
    # file_type_dir = file_type_dirs[ftype]
    # download_directory = download_directory_format.format(cam_dir, file_type_dir)
    remote_file_names = get_file_names(file_numbers, camera, ftype)
    cam_save_dir = os.path.join(save_dir, 'C{}'.format(camera+1))
    if not os.path.isdir(cam_save_dir):
        os.makedirs(cam_save_dir)
    # os.chdir(working_dir)
    with pysftp.Connection(ip, username=user, password=password) as sftp:
        for f in remote_file_names:
            # f = file_format.format(i, camera)
            exists = sftp.exists(f)
            save_f = os.path.join(cam_save_dir, os.path.basename(f).replace('.fits.ramp', '.ramp.fits'))
            if exists and (not os.path.isfile(save_f) or overwrite):
                print(f)
                get_with_retry(sftp, f, save_f, overwrite, save_dir, n_retry)
            elif not exists:
                print('{} not in server directory, skipping...'.format(f))
            elif os.path.isfile(save_f):
                print('{} already downloaded, skipping...'.format(f))


def download_data(
    date, save_dir='.', redownload=defaults['redownload'], overwrite=defaults['overwrite'], ftype=defaults['ftype'],
    ip=defaults['ip'], user=defaults['user'], password=defaults['password'],
    objname=defaults['objname'], objtype=defaults['objtype'], observer=defaults['observer'], chip=defaults['chip'],
    filter1=defaults['filter1'], filter2=defaults['filter2'], n_retry=defaults['n_retry']
):
    """
    Main interface function to download a data set

    Parameters
    ----------
    date: str
        date string, flexible formatting, example: '20230912' or '2023/09/12'
    save_dir: str, path
        directory where data will be saved
    redownload: bool

    overwrite
    ftype
    ip
    user
    password
    objname
    objtype
    observer
    chip
    filter1
    filter2
    n_retry

    Returns
    -------

    """
    datetime = to_datetime(date)
    date = datetime.strftime('%Y-%m-%d')

    if not os.path.isdir(save_dir):
        os.makedirs(save_dir)
    os.chdir(save_dir)

    log_file_df = get_log_file(date, ip, user, password, redownload, save_dir, n_retry)
    log_file_df = filter_df(log_file_df, objname, objtype, observer, filter1, filter2)
    file_numbers = get_file_numbers(log_file_df, ftype)
    cameras = [int(j)-1 for j in chip]
    for camera in cameras:
        download_camera_data(camera, file_numbers, save_dir, overwrite, ftype, user, ip, password, n_retry)


def main():
    parser = ArgumentParser()

    parser.add_argument('date', type=str, help='[str] log file date yyyymmdd')
    parser.add_argument(
        '-d', '--directory', type=str, help="[str], optional, save directory", default=defaults['save_dir']
    )
    parser.add_argument(
        '-r', '--redownload', action='store_true', help='overwrite the previously downloaded log file',
        default=defaults['redownload']
    )
    parser.add_argument(
        '-o', '--overwrite', action='store_true', help='overwrite existing fits files', default=defaults['overwrite']
    )
    parser.add_argument(
        '-f', '--ftype', type=str, help="[str], optional, file type desired: {}".format(file_types_str),
        default=defaults['ftype']
    )
    parser.add_argument(
        '-i', '--ip', type=str, help="[str], optional, data server ip address", default=defaults['ip']
    )
    parser.add_argument(
        '-u', '--user', type=str, help="[str], optional, connection user name", default=defaults['user']
    )
    parser.add_argument(
        '-p', '--password', type=str, help="[str], optional, connection password", default=defaults['password']
    )
    parser.add_argument(
        '--objname', type=str, help="[str], optional, object name, only downloads files with specified object name",
        default=defaults['objname']
    )
    parser.add_argument(
        '--objtype', type=str, help="[str], optional, object type, only downloads files with specified object type",
        default=defaults['objtype']
    )
    parser.add_argument(
        '--observer', type=str, help="[str], optional, observer, only downloads files with specified observer",
        default=defaults['observer']
    )
    parser.add_argument(
        '-c', '--chip', type=str,
        help="[str], optional, only downloads files for the specified detectors in format <chip_a>,<chip_b>,... e.g. 1,3,4",
        default=','.join(defaults['chip'])
    )
    parser.add_argument(
        '--filter1', type=str,
        help="[str], optional, filter 1, only downloads files with specified filter 1: {}".format(filter_1_options),
        default=None
    )
    parser.add_argument(
        '--filter2', type=str,
        help="[str], optional, filter 1, only downloads files with specified filter 2: {}".format(filter_2_options),
        default=None
    )
    parser.add_argument(
        '--n_retry', type=int,
        help="[int], optional, number of attempts to download a file, default: {}".format(defaults['n_retry']),
        default=defaults['n_retry']
    )
    args = parser.parse_args()

    download_data(
        args.date, args.directory, args.redownload, args.overwrite, args.ftype, args.ip, args.user, args.password,
        args.objname, args.objtype, args.observer, args.chip.split(','), args.filter1, args.filter2, args.n_retry
    )


if __name__ == '__main__':
    main()
