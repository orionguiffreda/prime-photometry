#!/usr/bin/python
"""
Author: Joe Durbak, durbak.3@gmail.com
"""
import os
import io

from astropy.io.fits import PrimaryHDU
from pandas import read_csv, to_datetime
from astropy.io import fits
import numpy as np

from photometrus.utils.defaults import FILE_DEFAULTS
from photometrus.settings import GET_DATA_SETTINGS
from photometrus.ramp import do_ramp, calc_darklim, reduce_image_from_file_list

defaults = FILE_DEFAULTS.copy()


def path_replace(path):
    path_replace_format = (
        ('%SERIAL%', '{0:08d}'),
        ('%SERIAL_TRUNC%', '{1:08d}'),
        ('%CHIP%', '{0}'),
    )
    for replace_format in path_replace_format:
        path = path.replace(replace_format[0], replace_format[1])
    return path


file_types = ['raw_fz', 'ramp', 'raw']
file_types_str = ','.join(file_types)

filter_1_options = ','.join(['NB', 'Open', 'Z', 'Dark'])
filter_2_options = ','.join(['Open', 'Y', 'J', 'H'])

log_folder_location = os.path.join(GET_DATA_SETTINGS['log_folder_location'], 'ramp_fit_log_{}.dat')

file_prefix = '{0:08d}'

backup_lists = {
    'ramp': ['ramp_fz', 'regen'],
    'raw': ['raw_fz',],
    'raw_fz': ['raw',],
    'ramp_fz': ['ramp', 'regen'],
}

funpack_output_dir = path_replace(GET_DATA_SETTINGS['funpack_output_dir'])

remote_file_formats = GET_DATA_SETTINGS['remote_file_formats']
for k,v in remote_file_formats.items():
    for i, fmt in enumerate(v):
        remote_file_formats[k][i] = path_replace(fmt)

replace_list = GET_DATA_SETTINGS['replace_list']


class NoCalError(Exception):
    pass


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
    nframes = df['NFRAMES'].tolist()
    if 'raw' in ftype:
        new_f_numbers = []
        for start, frames in zip(f_numbers, nframes):
            for frame in range(int(frames)):
                new_f_numbers.append(start+frame)
        f_numbers = new_f_numbers
    return f_numbers, nframes


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


def get_temperature_date(datetime, temperature):
    """
    picks calibration date based on datetime and temperature
    """
    temperature_date_dict = GET_DATA_SETTINGS['temperature_date_dict']
    cutoff_dates = sorted(temperature_date_dict.keys())
    temp_dates = {}
    for cutoff_date in cutoff_dates:  # getting era with identical load files
        td = to_datetime(cutoff_date) - datetime
        if td.total_seconds() > 0:
            temp_dates = temperature_date_dict[cutoff_date]
            break
    if not temp_dates:  # raise an error if no dataset exists
        raise NoCalError
    temperature_keys = temp_dates.keys()
    temperatures = np.asarray(temperature_keys, dtype=float)
    nearest_temp = np.abs(temperatures-temperature).argmin()  # getting closest detector temperature with cal data
    date_list = temp_dates[temperature_keys[nearest_temp]]
    if len(date_list) == 1:
        return date_list[0]
    else:
        date_array = np.asarray([to_datetime(d) for d in date_list])
        nearest_date = np.abs(date_array-date_array).argmin()  # getting closest date with cal data
        return date_list[nearest_date]


def get_ramp_cal_files(header, base_dir=GET_DATA_SETTINGS['ramp_cal_directory']):
    """
    picks the appropriate calibration files based on the chip, temperature and date
    """
    chip = header['CHIP']
    temperature = header['TEMPDET']
    date = to_datetime(header['DATE'])
    date = get_temperature_date(date, temperature)
    coe_R = f"coe_R/coe_R.C{chip}.fits.cube.{date}"
    coe_D = f"coe_D/coe_D.C{chip}.fits.cube.{date}"
    darklim = f"darklim/darklim.C{chip}.fits.{date}"
    superbias = f"bias/superbias_minRMS.C{chip}.fits.{date}"
    satulim = f"satulim/satulim.C{chip}.fits.{date}"
    return [os.path.join(base_dir, f) for f in (superbias, satulim, coe_R, coe_D, darklim)]


def get_ramp_cal(header):
    cal_files = get_ramp_cal_files(header)
    superbias, satulim, coe_R, coe_D, darklim = [fits.getdata(f) for f in cal_files]
    coe_R_tr = np.transpose(coe_R, (1, 2, 0)).copy(order='C')
    coe_D_tr = np.transpose(coe_D, (1, 2, 0)).copy(order='C')
    Adarklim, Fdarklim = calc_darklim(coe_D_tr, darklim)
    satulim_flat = satulim[4:-4, 4:-4].flatten()
    satulim_thre = np.nanmean(satulim_flat) - 5 * np.nanstd(satulim_flat)
    mask = np.where(satulim > satulim_thre, 1, 0).astype(np.uint8)
    return superbias, satulim, mask, coe_R_tr, coe_D_tr, darklim, Adarklim, Fdarklim


def update_header(header, f0_num, fl_num):
    superbias, satulim, coe_R, coe_D, darklim =\
        [os.path.basename(f) for f in get_ramp_cal_files(header)]
    header["BZERO"] = 0
    header["BSCALE"] = 1
    header["COER"] = coe_R
    header["COED"] = coe_D
    header["DARKLIM"] = darklim
    header["SPBIAS"] = superbias
    header["SATULIM"] = satulim
    header["USEDNUMS"] = f0_num
    header["USEDNUML"] = fl_num
    return header


def make_ramp_fits(ramp, header, out_file, f0_num, fl_num):
    header = update_header(header, f0_num, fl_num)
    fits.HDUList([PrimaryHDU(ramp.astype(np.float32, header))]).writeto(out_file, overwrite=True)
    print(f"make {out_file}")


def regen_ramp(file_number, nframe, camera):
    output_dir = funpack_output_dir.format(camera+1, truncate_1000(file_number), '{:08d}.ramp.fits'.format(file_number))
    file_numbers = [file_number+i for i in range(nframe)]
    nframes = [nframe for i in range(nframe)]
    raw_files, missing_raw_files = get_file_names(
        file_numbers, nframes, backup_file_types=backup_lists['raw'], ftype='raw', cameras=(camera,)
    )
    raw_files = raw_files[0]
    raw_files.sort()
    missing_raw_files = missing_raw_files[0]
    if missing_raw_files:
        return 'does not exist'
    ext_dict = {'.fz': 1, '.fits': 0, '.ramp': 0}
    extension = ext_dict[os.path.splitext(raw_files[0])[1]]
    try:
        superbias, satulim, mask, coe_R_tr, coe_D_tr, darklim, Adarklim, Fdarklim = get_ramp_cal_files(fits.getheader(raw_files[0], ext=extension))
        header, ramp = do_ramp(raw_files, superbias, satulim, mask, coe_R_tr, coe_D_tr, darklim, Adarklim, Fdarklim, extension)
    except NoCalError:
        header, ramp = reduce_image_from_file_list(raw_files, hdu_ext=extension)
    make_ramp_fits(ramp, header, output_dir, file_numbers[0], file_numbers[-1])
    return output_dir


def get_backup_file(file_number, nframe, backup_types, camera, funpack_fz=defaults['funpack_fz']):
    for backup_type in backup_types:
        if backup_type == 'regen':
            backup_filename = regen_ramp(file_number, nframe, camera)
        else:
            backup_filename = get_file_name(file_number, backup_type, camera, funpack_fz=funpack_fz)
        if os.path.exists(backup_filename):
            return backup_filename
    raise FileNotFoundError


def get_file_names(
    file_numbers, nframes, backup_file_types, cameras=defaults['chip'], ftype=defaults['ftype'],
    funpack_fz=defaults['funpack_fz'],
):

    filenames = []
    missing_filenames = []
    cameras = [camera - 1 for camera in cameras]
    for camera in cameras:
        camera_filenames = [get_file_name(_n, ftype, camera, funpack_fz=funpack_fz) for _n in file_numbers]
        existing_camera_filenames = []
        missing_camera_filenames = []
        for f, nframe in zip(camera_filenames, nframes):
            if os.path.exists(f):
                existing_camera_filenames.append(f)
            else:
                try:
                    file_number = int(os.path.basename(f)[:8])
                    existing_camera_filenames.append(
                        get_backup_file(file_number, nframe, backup_file_types, camera, funpack_fz=funpack_fz)
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
    file_numbers, nframes = get_file_numbers(log_file_df, ftype)
    file_names, missing_file_names = get_file_names(file_numbers, nframes, backup_file_types, cameras, ftype, funpack_fz)
    return file_names, missing_file_names
