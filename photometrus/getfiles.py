#!/usr/bin/python
"""
Author: Joe Durbak, durbak.3@gmail.com
"""
import os
import io

from pandas import read_csv, to_datetime

from photometrus.utils.defaults import FILE_DEFAULTS as defaults
from photometrus.settings import GET_DATA_SETTINGS

file_types = ['raw_fz', 'ramp', 'raw']
file_types_str = ','.join(file_types)

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
    'raw_fz': 'C{0}.fits.fz',
    'ramp_fz': 'C{0}.fits.ramp.fz',
}

file_type_dirs = {
    'ramp': 'ramp',
    'raw': 'raw',
    'raw_fz': 'raw_fz',
    'ramp_fz': 'ramp_fz',
}

backup_lists = {
    'ramp': ['ramp_fz', 'regen'],
    'raw': ['raw_fz',],
    'raw_fz': ['raw',],
    'ramp_fz': ['ramp', 'regen'],
}

raw_fz_subdir = '{1:08d}'
funpack_output_dir = '/mnt/photometry/unarchived/C{0}/'+raw_fz_subdir

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
        '/'.join((
            '/nfs/xion{}/prime/raw_fz/C{}'.format(_i%2+3, _i+1),
            raw_fz_subdir,
            file_prefix+file_formats['raw_fz'].format(_i+1)
        )) for _i in range(len(cam_dirs))],
    'ramp_fz': [
        '/'.join((
            '/nfs/xion{}/prime/ramp_fz/C{}'.format(_i % 2 + 3, _i + 1),
            raw_fz_subdir,
            file_prefix + file_formats['ramp_fz'].format(_i + 1)
        )) for _i in range(len(cam_dirs))],
}

replace_list = GET_DATA_SETTINGS['replace_list']


def funpack(fpacked_file, funpacked_file):
    outdir = os.path.dirname(funpacked_file)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    os.system('funpack -O {} {}'.format(funpacked_file, fpacked_file))

def replace_bad_string(lines):
    replaced_lines = []
    for line in lines:
        fixed_line = line
        for replacement in replace_list:
            fixed_line = fixed_line.replace(replacement[0], replacement[1])
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
    if not read_lines[0].startswith('#'):
        read_lines.insert(
            0,
            '# DATEBEG, filename,OBJNAME,OBJTYPE,OBSERVER,CHIP,FILTER1,FILTER2,EXPTIME,EXPTIMEE,EXPTIMEC,NFRAMES,DITHTYP,DITHRAD,DITHPH,DITH_REP,TEMPDET,FOCUS,NINT,INT'
         )
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


def get_file_name(file_number, ftype, camera, funpack_fz=defaults['funpack_fz']):
    remote_file_format = remote_file_formats[ftype][camera]
    if ftype.endswith('fz'):
        file_name = remote_file_format.format(file_number, truncate_1000(file_number))
        if funpack_fz:
            if os.path.exists(file_name):
                output_dir = funpack_output_dir.format(camera+1, truncate_1000(file_number))
                funpack_file_name = os.path.basename(file_name).replace('.fz', '')
                funpack_file_name = os.path.join(output_dir, funpack_file_name)
                funpack(file_name, funpack_file_name)
                file_name = funpack_file_name
        return file_name
    else:
        return remote_file_format.format(file_number)


def regen_ramp(file_number, camera):
    return 'does not exist'  # TODO: incorporate Ryusei's data reduction tool


def get_backup_file(file_number, backup_types, camera, funpack_fz=defaults['funpack_fz']):
    for backup_type in backup_types:
        if backup_type == 'regen':
            backup_filename = regen_ramp(file_number, camera)
        else:
            backup_filename = get_file_name(file_number, backup_type, camera, funpack_fz=funpack_fz)
        if os.path.exists(backup_filename):
            return backup_filename
    raise FileNotFoundError


def get_file_names(
    file_numbers, backup_file_types, cameras=defaults['chip'], ftype=defaults['ftype'],
    funpack_fz=defaults['funpack_fz'],
):

    filenames = []
    missing_filenames = []
    cameras = [camera - 1 for camera in cameras]
    for camera in cameras:
        camera_filenames = [get_file_name(_n, ftype, camera, funpack_fz=funpack_fz) for _n in file_numbers]
        existing_camera_filenames = []
        missing_camera_filenames = []
        for f in camera_filenames:
            if os.path.exists(f):
                existing_camera_filenames.append(f)
            else:
                try:
                    file_number = int(os.path.basename(f)[:8])
                    existing_camera_filenames.append(
                        get_backup_file(file_number, backup_file_types, camera, funpack_fz=funpack_fz)
                    )
                except FileNotFoundError:
                    missing_camera_filenames.append(f)
        filenames.append(existing_camera_filenames)
        missing_filenames.append(missing_camera_filenames)
    return filenames, missing_filenames


def get_data_files(
    date, ftype=defaults['ftype'],
    objname=defaults['objname'], objtype=defaults['objtype'], observer=defaults['observer'], cameras=defaults['chip'],
    filter1=defaults['filter1'], filter2=defaults['filter2'], use_backups=True, backup_file_types=None,
    funpack_fz=defaults['funpack_fz'],
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
    use_backups
    backup_file_types
    funpack_fz

    Returns filenames, missing_filenames
    -------

    """
    if backup_file_types is None:
        backup_file_types = backup_lists[ftype]
    if not use_backups:
        backup_file_types = tuple()
    datetime = to_datetime(date)
    date = datetime.strftime('%Y-%m-%d')
    log_file_df = get_log_file(date)
    log_file_df = filter_df(log_file_df, objname, objtype, observer, filter1, filter2)
    file_numbers = get_file_numbers(log_file_df, ftype)
    file_names, missing_file_names = get_file_names(file_numbers, backup_file_types, cameras, ftype, funpack_fz)
    return file_names, missing_file_names
