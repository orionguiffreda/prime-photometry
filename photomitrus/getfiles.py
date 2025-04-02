#!/usr/bin/python
"""
Author: Joe Durbak, durbak.3@gmail.com
"""
import os
import io

from pandas import read_csv, to_datetime

file_types = ['raw_fz', 'ramp', 'raw']
file_types_str = ','.join(file_types)

defaults = dict(
    save_dir='.', redownload=False, overwrite=False, ftype='ramp',
    objname=None, objtype=None, observer=None, chip=(1,2,3,4), filter1=None, filter2=None,
)
filter_1_options = ','.join(['NB', 'Open', 'Z', 'Dark'])
filter_2_options = ','.join(['Open', 'Y', 'J', 'H'])

log_folder_location = '/nfs/prime01/work/LOG/Ramp_LOG/ramp_fit_log_{}.dat'

cam_dirs = ['home', 'xion2', 'xion3', 'xion4']
cam_dirs_dict_all = {(cam_num+1): cam_dir for cam_num, cam_dir in enumerate(cam_dirs)}
download_directory_format = '/nfs/{}/prime/Data/{}'
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
        '/'.join(('/nfs/xion3/prime/raw_fz/C1', raw_fz_subdir, file_prefix+file_formats['raw_fz'].format(1))),
        '/'.join(('/nfs/xion4/prime/raw_fz/C2', raw_fz_subdir, file_prefix+file_formats['raw_fz'].format(2))),
        '/'.join(('/nfs/xion3/prime/raw_fz/C3', raw_fz_subdir, file_prefix+file_formats['raw_fz'].format(3))),
        '/'.join(('/nfs/xion4/prime/raw_fz/C4', raw_fz_subdir, file_prefix+file_formats['raw_fz'].format(4))),
    ]
}

replace_list = [
    'all sky grid'
]


def replace_bad_string(lines):
    replaced_lines = []
    for line in lines:
        fixed_line = line
        for replacement in replace_list:
            fixed_line = fixed_line.replace(replacement, replacement.replace(' ', '_'))
        replaced_lines.append(fixed_line)
    return replaced_lines


def get_log_file(date):
    # log_link = log_link_format.format(date)
    log_link = log_folder_location.format(date)
    print(log_link)
    # output_file, http_message = request.urlretrieve(log_link, filename)

    with open(log_link, 'r') as _f:
        read_lines = _f.readlines()
    cleaned_lines = [
        line for line in read_lines if not (
                line.startswith('Warning') or line.startswith('ERROR') or line.startswith('#') or (line.find('!')!=-1)
        )
    ]
    first_line = ' '.join([field.strip() for field in read_lines[0][1:].split(',')]) + '\n'
    cleaned_lines = replace_bad_string(cleaned_lines)
    cleaned_lines = [first_line] + cleaned_lines
    # print(first_line)
    df = read_csv(io.StringIO('\n'.join(cleaned_lines)), sep=' ', on_bad_lines='warn')
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


def get_file_names(file_numbers, cameras=defaults['chip'], ftype=defaults['ftype']):
    filenames = []
    missing_filenames = []
    cameras = [camera - 1 for camera in cameras]
    for camera in cameras:
        remote_file_format = remote_file_formats[ftype][camera]
        if ftype.endswith('fz'):
            camera_filenames = [remote_file_format.format(_n, truncate_1000(_n)) for _n in file_numbers]
        else:
            camera_filenames = [remote_file_format.format(_n) for _n in file_numbers]
        existing_camera_filenames = []
        missing_camera_filenames = []
        for f in camera_filenames:
            if os.path.exists(f):
                existing_camera_filenames.append(f)
            else:
                missing_camera_filenames.append(f)
        filenames.append(existing_camera_filenames)
        missing_filenames.append(missing_camera_filenames)
    # print(filenames)
    return filenames, missing_filenames


# def download_camera_data(
#     camera, file_numbers, save_dir=defaults['save_dir'], ftype=defaults['ftype'],
# ):
    # file_format = file_formats[ftype]
    # cam_dir = cam_dirs[camera]
    # file_type_dir = file_type_dirs[ftype]
    # download_directory = download_directory_format.format(cam_dir, file_type_dir)
    # remote_file_names = get_file_names(file_numbers, camera, ftype)
    # cam_save_dir = os.path.join(save_dir, 'C{}'.format(camera+1))
    # if not os.path.isdir(cam_save_dir):
    #     os.makedirs(cam_save_dir)
    # os.chdir(working_dir)
    # return


def get_data_files(
    date, ftype=defaults['ftype'],
    objname=defaults['objname'], objtype=defaults['objtype'], observer=defaults['observer'], cameras=defaults['chip'],
    filter1=defaults['filter1'], filter2=defaults['filter2']
):
    """
    Main interface function to download a data set

    Parameters
    ----------
    date: str
        date string, flexible formatting, example: '20230912' or '2023/09/12'
    ftype
    objname
    objtype
    observer
    cameras
    filter1
    filter2

    Returns filenames, missing_filenames
    -------

    """
    datetime = to_datetime(date)
    date = datetime.strftime('%Y-%m-%d')
    log_file_df = get_log_file(date)
    log_file_df = filter_df(log_file_df, objname, objtype, observer, filter1, filter2)
    file_numbers = get_file_numbers(log_file_df, ftype)
    file_names, missing_file_names = get_file_names(file_numbers, cameras, ftype)
    return file_names, missing_file_names
