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
import fitsio

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
backup_log_folder_location = os.path.join(GET_DATA_SETTINGS['backup_log_folder_location'], 'ramp_log_C1.{}.dat')  # /home/prime/hamada/log

file_prefix = '{0:08d}'

backup_lists = {
    'ramp': ['real_time_ramp', 'ramp_fz', 'regen'],
    'raw': ['raw_fz',],
    'raw_fz': ['raw',],
    'ramp_fz': ['ramp', 'real_time_ramp', 'regen'],
    'real_time_ramp': ['ramp', 'ramp_fz', 'regen'],
    'regen': ['regen']
}

funpack_output_dir = path_replace(GET_DATA_SETTINGS['funpack_output_dir'])

remote_file_formats = GET_DATA_SETTINGS['remote_file_formats']
for _k, _v in remote_file_formats.items():
    for _i, fmt in enumerate(_v):
        remote_file_formats[_k][_i] = path_replace(fmt)

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
    if not os.path.exists(log_link):
        log_link = backup_log_folder_location.format(date.replace('-', ''))
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
    if ftype.endswith('fz') or ftype.startswith('regen'):
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
        td = datetime - to_datetime(cutoff_date)
        if td.total_seconds() > 0:
            temp_dates = temperature_date_dict[cutoff_date]
            break
    if not temp_dates:  # raise an error if no dataset exists
        raise NoCalError
    temperature_keys = temp_dates.keys()
    temperatures = np.asarray(list(temperature_keys), dtype=float)
    nearest_temp = np.abs(temperatures-temperature).argmin()  # getting closest detector temperature with cal data
    date_list = temp_dates['{:.1f}'.format(temperatures[nearest_temp])]
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
    date = to_datetime(header['DATE'])
    try:
        temperature = header['TEMPDET']
    except KeyError:
        try:
            date_str = date.strftime('%Y%m%d')
            temperature = GET_DATA_SETTINGS['date_temperature_dict'][date_str]
        except KeyError:
            print(GET_DATA_SETTINGS)
            raise KeyError('TEMPDET header is missing from file, and no backup temperature is available. An update to photometrus.json5, "date_temperature_dict" is required for date {} to rereduce this data'.format(date_str))
    date = get_temperature_date(date, temperature)
    coe_R = f"coe_R/coe_R.C{chip}.fits.cube.{date}"
    coe_D = f"coe_D/coe_D.C{chip}.fits.cube.{date}"
    darklim = f"darklim/darklim.C{chip}.fits.{date}"
    superbias = f"bias/superbias_minRMS.C{chip}.fits.{date}"
    satulim = f"satulim/satulim.C{chip}.fits.{date}"
    cal_files = [os.path.join(base_dir, f) for f in (superbias, satulim, coe_R, coe_D, darklim)]
    sb = cal_files[0]
    if not os.path.isfile(sb):
        cal_files[0] = sb.replace('bias/superbias_minRMS.', 'bias/super_bias')
    return [os.path.join(base_dir, f) for f in (superbias, satulim, coe_R, coe_D, darklim)]


def get_ramp_cal(header):
    cal_files = get_ramp_cal_files(header)
    superbias, satulim, coe_R, coe_D, darklim = [fits.getdata(f).astype(np.float64) for f in cal_files]
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
    out_dir = os.path.dirname(out_file)
    if not os.path.isdir(out_dir) and out_dir != '':
        os.makedirs(out_dir)
    if not isinstance(header, fits.Header):  # nlc ramp fit code uses fitsio instead of astropy.io.fits
        fitsio.write(out_file, ramp.astype(np.float32), header=header, clobber=True)
    else:
        fits.HDUList([PrimaryHDU(ramp.astype(np.float32), header)]).writeto(out_file, overwrite=True)
    print(f"make {out_file}")


def regen_ramp(file_number, nframe, camera):
    output_dir = funpack_output_dir.format(camera+1, truncate_1000(file_number), '{:08d}.ramp.fits'.format(file_number))
    file_numbers = [file_number+i for i in range(nframe)]
    nframes = [nframe for i in range(nframe)]
    print('regenerating_ramp', file_number, nframe, camera)
    raw_files, missing_raw_files, filtered = get_file_names(
        file_numbers, nframes, backup_file_types=backup_lists['raw'], ftype='raw', cameras=(camera+1,)
    )
    print('raw_files')
    print(raw_files)
    raw_files = raw_files[0]
    print('selected raw_files')
    print(raw_files)
    raw_files.sort()
    missing_raw_files = missing_raw_files[0]
    if missing_raw_files:
        print('missing files: {}'.format(missing_raw_files))
        return 'does not exist'
    ext_dict = {'.fz': 1, '.fits': 0, '.ramp': 0}
    extension = ext_dict[os.path.splitext(raw_files[0])[1]]
    output_file = os.path.join(output_dir, os.path.basename(raw_files[0]).replace('.fz', '').replace('.fits', '.ramp.fits'))
    if os.path.isfile(output_file):
        try:
            fits.getdata(output_file)
            return output_file
        except OSError:
            print('{} exists, but seems to be corrupted. Regenerating file...')
    try:
        print('regenerating ramp: {}'.format(output_file))
        superbias, satulim, mask, coe_R_tr, coe_D_tr, darklim, Adarklim, Fdarklim = get_ramp_cal(fits.getheader(raw_files[0], ext=extension))
        header, ramp = do_ramp(raw_files, superbias, satulim, mask, coe_R_tr, coe_D_tr, darklim, Adarklim, Fdarklim, extension)
    except NoCalError:
        header, ramp = reduce_image_from_file_list(raw_files, hdu_ext=extension)
    make_ramp_fits(ramp, header, output_file, file_numbers[0], file_numbers[-1])
    return output_file


def get_backup_file(file_number, nframe, backup_types, camera, funpack_fz=defaults['funpack_fz']):
    for backup_type in backup_types:
        if backup_type == 'regen':
            backup_filename = regen_ramp(file_number, nframe, camera)
        else:
            backup_filename = get_file_name(file_number, backup_type, camera, funpack_fz=funpack_fz)
        if os.path.exists(backup_filename):
            return backup_filename
    raise FileNotFoundError


def header_filter_file(filename, header_filter):
    hdr = fits.getheader(filename)
    for k, v in header_filter.items():
        if str(hdr[k]).strip() == v:
            continue
        else:
            return True
    return False


def get_file_names(
    file_numbers, nframes, backup_file_types=None, cameras=defaults['chip'], ftype=defaults['ftype'],
    funpack_fz=defaults['funpack_fz'], header_filter=None
):
    if backup_file_types is None:
        backup_file_types = backup_lists[ftype]
    filenames = []
    missing_filenames = []
    filtered_filenames = []
    cameras = [camera - 1 for camera in cameras]
    for camera in cameras:
        camera_filenames = [get_file_name(_n, ftype, camera, funpack_fz=funpack_fz) for _n in file_numbers]
        existing_camera_filenames = []
        missing_camera_filenames = []
        for f, nframe in zip(camera_filenames, nframes):
            if os.path.exists(f):
                if header_filter is None or not header_filter_file(f, header_filter):
                    existing_camera_filenames.append(f)
                else:
                    filtered_filenames.append(f)
            else:
                try:
                    file_number = int(os.path.basename(f)[:8])
                    backup_file = get_backup_file(file_number, nframe, backup_file_types, camera, funpack_fz=funpack_fz)
                    if header_filter is None or not header_filter_file(backup_file, header_filter):
                        existing_camera_filenames.append(backup_file)
                    else:
                        filtered_filenames.append(backup_file)
                except FileNotFoundError:
                    missing_camera_filenames.append(f)
        filenames.append(existing_camera_filenames)
        missing_filenames.append(missing_camera_filenames)
    return filenames, missing_filenames, filtered_filenames


def get_data_files(
    date, ftype=defaults['ftype'],
    objname=defaults['objname'], objtype=defaults['objtype'], observer=defaults['observer'], cameras=defaults['chip'],
    filter1=defaults['filter1'], filter2=defaults['filter2'], use_backups=True, backup_file_types=None,
    funpack_fz=defaults['funpack_fz'], header_filter=defaults['header_filter'],
    # skip_run_numbers=None
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
    header_filter: expects string of comma deliminated
    skip_run_numbers: expects string of comma deliminated integers

    Returns filenames, missing_filenames
    -------

    """
    if backup_file_types is None:
        backup_file_types = backup_lists[ftype]
    if not use_backups:
        backup_file_types = tuple()
    # if skip_run_numbers is not None:
    #     skip_run_numbers = np.asarray([int(i) for i in skip_run_numbers.split(',')])
    if header_filter is not None:
        tuples = [i.split(':') for i in header_filter.split(',')]
        header_filter = {k.upper().strip(): v.strip() for k, v in tuples}
    datetime = to_datetime(date)
    date = datetime.strftime('%Y-%m-%d')
    log_file_df = get_log_file(date)
    log_file_df = filter_df(log_file_df, objname, objtype, observer, filter1, filter2)
    file_numbers, nframes = get_file_numbers(
        log_file_df, ftype,
        # skip_run_numbers
    )
    file_names, missing_file_names, filtered_file_name = get_file_names(
        file_numbers, nframes, backup_file_types, cameras, ftype, funpack_fz, header_filter
    )
    print('header_filtered', filtered_file_name)
    return file_names, missing_file_names
