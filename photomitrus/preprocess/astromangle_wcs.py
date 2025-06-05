import os
from fnmatch import fnmatch
import shutil
from datetime import datetime as dt

import argparse
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits
from astropy.wcs.utils import fit_wcs_from_points
from astropy.table import Table

from photomitrus.settings import gen_config_file_name

_mesh_file_dir = gen_config_file_name('wcs')
_wcs_matrix_tranlation = {
    'PC1_1': 'CD1_1',
    'PC1_2': 'CD1_2',
    'PC2_1': 'CD2_1',
    'PC2_2': 'CD2_2',
}
_default_sip_degree = 4
_default_downsample = 16


def get_sep_rot_table(chip_number, mesh_file_dir=_mesh_file_dir):
    sep_rot_files = [os.path.join(mesh_file_dir, 'C{}_{}.fits'.format(chip_number, s)) for s in ('sep', 'rot')]
    sep_rot_hdus = [fits.open(f)[0] for f in sep_rot_files]
    sep_array, rot_array = [hdu.data for hdu in sep_rot_hdus]
    header = sep_rot_hdus[0].header
    mesh_step = header['MESHSTEP']
    x_indices = np.ones(sep_array.shape, dtype=int)
    x_indices[:] = np.arange(sep_array.shape[0]) * mesh_step
    y_indices = x_indices.transpose()
    table = {
        'sep': sep_array.ravel(),
        'rot': rot_array.ravel(),
        'x': x_indices.ravel(),
        'y': y_indices.ravel(),
    }
    return Table(table)


def gen_ra_dec(ra_tel, dec_tel, rot, chip_number, mesh_file_dir=_mesh_file_dir):
    tel_coords = SkyCoord(ra_tel, dec_tel, unit=(u.degree, u.degree))
    table = get_sep_rot_table(chip_number, mesh_file_dir)
    table['rot'] += rot
    table['coord'] = tel_coords.directional_offset_by(table['rot']*u.deg, table['sep']*u.deg)
    return table


def calculate_wcs(
    ra_tel, dec_tel, rot, chip_number, mesh_file_dir=_mesh_file_dir, sip_degree=_default_sip_degree,
    downsample=_default_downsample,
):
    table = gen_ra_dec(ra_tel, dec_tel, rot, chip_number, mesh_file_dir)
    table = table[::downsample]
    wcs = fit_wcs_from_points((table['x'], table['y']), table['coord'], projection='TAN', sip_degree=sip_degree)
    return wcs


def calculate_wcs_header(
    file_header, rot_val, mesh_file_dir=_mesh_file_dir, sip_degree=_default_sip_degree, downsample=_default_downsample
):
    ra_tel = file_header['RA-D']
    dec_tel = file_header['DEC-D']
    rot = file_header['ROTOFF']  # - 48
    try:
        rot = int(rot)
    except ValueError:
        rot = rot_val
    chip_number = file_header['CHIP']
    wcs = calculate_wcs(ra_tel, dec_tel, rot, chip_number, mesh_file_dir, sip_degree, downsample=downsample)
    return wcs


def calculate_wcs_file(
    filename, rot_val, mesh_file_dir=_mesh_file_dir, sip_degree=_default_sip_degree, downsample=_default_downsample
):
    file_header = fits.getheader(filename)
    return calculate_wcs_header(file_header, rot_val, mesh_file_dir, sip_degree, downsample=downsample)


def update_ra_dec(
    fits_file, rot_val, mesh_file_dir=_mesh_file_dir, sip_degree=_default_sip_degree, downsample=_default_downsample
):
    print(fits_file)
    start_time = dt.now()
    hdr_check = fits.getheader(fits_file)
    if 'CHECKSUM' in hdr_check and len(hdr_check) <= 10:
        print('File is fpack compressed...')
        os.system('funpack -F %s' % fits_file)
        print('funpacked!')
    with fits.open(fits_file, 'update') as f:
        wcs = calculate_wcs_header(f[0].header, rot_val, mesh_file_dir, sip_degree, downsample=downsample)
        wcs_header = wcs.to_header(relax=True)
        for k, v in _wcs_matrix_tranlation.items():
            wcs_header.set(v, wcs_header[k], comment=wcs_header.comments[k])  # , before='CDELT1')
            del wcs_header[k]
        f[0].header.update(wcs_header)
        # print(fitsheader)
        # f[0].data = f[0].data[4:-4, 4:-4]
    end_time = dt.now()
    print('wcs gen time time:', (end_time - start_time).total_seconds())


def get_files(directory):
    ls = os.listdir(directory)
    print(ls)
    fits_files = [f for f in ls if fnmatch(f, '????????C?.fits') or fnmatch(f, '????????C?.fits.ramp') or
                  fnmatch(f, '????????C?.ramp.new') or fnmatch(f, '????????C?.ramp.fits')
                  or fnmatch(f, '????????C?.sky.flat.fits')]
    print(fits_files)
    return fits_files


def update_ra_dec_directory(directory, rot_val, downsample=_default_downsample):
    fits_files = [os.path.join(directory, f) for f in get_files(directory)]
    for fits_file in fits_files:
        update_ra_dec(fits_file, rot_val, downsample=downsample)


def update_ra_dec_move_directory(input_dir, output_dir, rot_val, downsample=_default_downsample):
    for f in sorted(os.listdir(input_dir)):
        if f.endswith('.ramp.fits'):
            origpath = os.path.join(input_dir, f)
            if f.endswith('ramp.fits.fz'):
                fnewname = f.replace('.ramp.fits.fz', '.ramp.new')
            else:
                fnewname = f.replace('.ramp.fits', '.ramp.new')
            newpath = os.path.join(output_dir, fnewname)
            shutil.copyfile(origpath, newpath)
            update_ra_dec(newpath, rot_val, downsample=downsample)
            print('%s updated, renamed, and moved!' % fnewname)


def update_ra_dec_list(input_list, output_dir, rot_val, downsample=_default_downsample):
    for filepath in input_list:
        if filepath.endswith('.fz'):
            filename = filepath[-23:]
            filenewname = filename.replace('.fits.ramp.fz', '.ramp.new')
        else:
            filename = filepath[-20:]
            filenewname = filename.replace('.fits.ramp', '.ramp.new')
        filenewpath = os.path.join(output_dir, filenewname)
        shutil.copyfile(filepath, filenewpath)
        update_ra_dec(filenewpath, rot_val, downsample=downsample)


def astrom_angle(input_field, output_dir, rot_val=48, downsample=_default_downsample):
    if type(input_field) is list:
        update_ra_dec_list(input_field, output_dir, rot_val, downsample)
    elif os.path.isdir(input_field):
        if output_dir:
            update_ra_dec_move_directory(input_field, output_dir, rot_val, downsample)
        if not output_dir:
            update_ra_dec_directory(input_field, rot_val, downsample)
    elif os.path.isfile(input_field):
        update_ra_dec(input_field, rot_val, downsample=downsample)


def main():
    parser = argparse.ArgumentParser(description='alternative to astrometry.net, generates initial astrometry on imgs'
                                                 ' using telescope pointing and corner positions - new')
    parser.add_argument('-input', type=str, help='[str] input path or single file for ramp images, put only this field if you '
                                                 'want to generate astrometry w/o changing names or path')
    parser.add_argument('-output', type=str, help='[str] output path for images w/ astrometry (ramp.new), to '
                                                  'be used in pipeline',default=None)
    parser.add_argument('-rot_val', type=float, help='[float] Use if you want to input a custom rotation value '
                                                     '(the default is 48 deg)', default=48)
    parser.add_argument('-downsample', type=int, help='[float] downsample for wcs grid', default=_default_downsample)
    args, unknown = parser.parse_known_args()
    astrom_angle(args.input, args.output, args.rot_val, args.downsample)


if __name__ == "__main__":
    main()
