#!/usr/bin/env python
# Calculates RA and DEC given separation and angle for each detector
#
# Usage:
#    update single file:
#        python angle.py 00087801C3.fits
#        python angle.py 00087801C3.fits.ramp
#    update all ????????C?.fits or ????????C?.fits.ramp files in a directory:
#        python angle.py /path/to/data/directory  # note: this skips files that have already been updated
#    get ra,dec strings back based on input ra,dec,rot:
#        python angle.py <ra> <dec> <rot>
#        python angle.py 08:20:09.400 -46:40:34.300 268.6

import os
import sys
import shutil
from fnmatch import fnmatch
import warnings

import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits
from astropy.wcs.utils import fit_wcs_from_points
import argparse
from astropy.utils.exceptions import AstropyWarning
#%%
warnings.simplefilter('ignore', category=AstropyWarning)

tel_angle_corr = -48
# rot_val = 48       # default is 48, change if the rot of the telescope is not default (172800) for your obs
# rot val is default at 48, if you want to use a custom val, use the -rot_val flag (run this script w/ the -h flag to see how it works)

offset_dict = {
    '2': {
        'separation': 0.44850763816602646,
        'angle': tel_angle_corr + 45.66331831019727,
        'rotation': 270.25 + tel_angle_corr,
        'corners': {
            'separation': (0.8440546940520994, 0.5992982059603985, 0.04713389920612984, 0.6021396275144256),
            'angle': (45.80165572567937+tel_angle_corr, 87.63450103257266+tel_angle_corr, 43.379142144408796+tel_angle_corr, 3.939995837173418+tel_angle_corr),
        }
    },
    '4': {
        'separation': 0.45007521998044664,
        'angle': tel_angle_corr + 135.1541758426945,
        'rotation': 90.10 + tel_angle_corr,
        'corners': {
            'separation': (0.601016428757574, 0.04901034640323484, 0.6029977744583211, 0.8456369866921561),
            'angle': (176.96558950424802+tel_angle_corr, 133.66527108136722+tel_angle_corr, 93.4803481359835+tel_angle_corr, 135.2329686892295+tel_angle_corr),
        }
    },
    '3': {
        'separation': 0.44169840857795223,
        'angle': tel_angle_corr + 225.56641620402425,
        'rotation': 91.41 + tel_angle_corr,
        'corners': {
            'separation': (0.837484328429415, 0.5967910467757283, 0.040129482544494374, 0.5947483907048631),
            'angle': (225.4819246084007+tel_angle_corr, 267.77509433798167+tel_angle_corr, 227.13371609154987+tel_angle_corr, 183.21310183951354+tel_angle_corr),
        }
    },
    '1': {
        'separation': 0.4457960233925488,
        'angle': tel_angle_corr + 315.7002683305685,
        'rotation': 270.56 + tel_angle_corr,
        'corners': {
            'separation': (0.5989415235089534, 0.044278559406327196, 0.5985614992422466, 0.8416094636073621),
            'angle': (357.70292702168376+tel_angle_corr, 315.8881546490197+tel_angle_corr, 273.6751062684781+tel_angle_corr, 315.6959537354861+tel_angle_corr),
        }
    },
    'corner_x_y': [(0, 0), (0, 4096), (4096, 4096), (4096, 0)]
}


def calc_offset(ra, dec, angle, sep):
    radec = SkyCoord(ra, dec, frame="fk5", unit=(u.hourangle, u.deg))
    # print (radec.ra.hms, radec.dec.dms)
    coord = radec.directional_offset_by(angle, sep)
    return coord


def calc_offset_detector(ra, dec, rot, chip):
    angle_off = offset_dict[chip]['angle'] + float(rot)
    separation_off = offset_dict[chip]['separation']
    chip_coords = calc_offset(ra, dec, angle_off * u.deg, separation_off * u.deg)
    new_rotation = (rot + offset_dict[chip]['rotation']) % 360
    return chip_coords, new_rotation


def calc_corners_detector(ra, dec, rot, chip):
    corner_coords = []
    for separation, angle in zip(offset_dict[chip]['corners']['separation'], offset_dict[chip]['corners']['angle']):
        angle_off = angle + rot
        corner_coords.append(calc_offset(ra, dec, angle_off * u.deg, separation * u.deg))
    return corner_coords


def gen_wcs(ra, dec, rot, chip):
    # wcs = WCS(naxis=2)
    # wcs.wcs.crpix = [image_width/2, image_height/2]
    # wcs.wcs.crval = [ra, dec]
    # wcs.wcs.cdelt = [-pixel_scale, -pixel_scale]
    # wcs.wcs.ctype = ["RA---AIR", "DEC--AIR"]
    corner_coords = calc_corners_detector(ra, dec, rot, chip)
    xy_array = np.asarray(offset_dict['corner_x_y']).transpose()
    wcs = fit_wcs_from_points((xy_array[0], xy_array[1]), SkyCoord(corner_coords),projection='TAN')
    return wcs


def get_files(directory):
    ls = os.listdir(directory)
    print(ls)
    fits_files = [f for f in ls if fnmatch(f, '????????C?.fits') or fnmatch(f, '????????C?.fits.ramp') or
                  fnmatch(f, '????????C?.ramp.new') or fnmatch(f, '????????C?.ramp.fits')
                  or fnmatch(f, '????????C?.sky.flat.fits')]
    print(fits_files)
    return fits_files


def update_ra_dec(fits_file, rot_val):
    print(fits_file)
    with fits.open(fits_file, 'update') as f:
        header = f[0].header
        if not header.get('COORD_UP', False):
            ra = header['RA']
            dec = header['DEC']
            rad = header['RA-D']
            decd = header['DEC-D']
            rotation = header['ROTOFF']
            try:
                float(rotation)
            except ValueError:
                rotation = rot_val  #48  # manually adjust value for now! until issue is fixed!
            # rotation = 48
            rot = header['ROT']
            chip = str(header['CHIP'])
            chip_coords, new_rotation = calc_offset_detector(ra, dec, rotation, chip)

            header['RA-D'] = chip_coords.ra.deg
            header['DEC-D'] = chip_coords.dec.deg
            header['RA'] = '{:02d}:{:02d}:{:.03f}'.format(
                int(chip_coords.ra.hms[0]), int(chip_coords.ra.hms[1]), chip_coords.ra.hms[2]
            )
            header['DEC'] = '{:02d}:{:02d}:{:02.03f}'.format(
                int(chip_coords.dec.dms[0]), abs(int(chip_coords.dec.dms[1])), abs(chip_coords.dec.dms[2])
            )
            header['ROT'] = new_rotation
            header['RA-TEL'] = ra
            header['DEC-TEL'] = dec
            header['ROTOFTEL'] = rot
            header['RA-D-TEL'] = rad
            header['DECD-TEL'] = decd
            header['COORD_UP'] = True

            for i, corner in enumerate(calc_corners_detector(ra, dec, rotation, chip)):
                header['RA-CNR{}D'.format(i)] = corner.ra.deg
                header['DECCNR{}D'.format(i)] = corner.dec.deg
                header['RA-CNR{}'.format(i)] = '{:02d}:{:02d}:{:.03f}'.format(
                    int(corner.ra.hms[0]), int(corner.ra.hms[1]), corner.ra.hms[2]
                )
                header['DECCNR{}'.format(i)] = '{:02d}:{:02d}:{:02.03f}'.format(
                    int(corner.dec.dms[0]), abs(int(corner.dec.dms[1])), abs(corner.dec.dms[2])
                )
            print(
                'Updated ra,dec,rotation: {},{},{}'.format(
                    str(chip_coords.ra), str(chip_coords.dec), new_rotation)
            )

            header.update(gen_wcs(ra, dec, rotation, chip).to_header(relax=True))
        else:
            print('Already updated')

    #update fits file pc matrix to cd matrix for reading by scamp

    #fits.delval(fits_file, 'CD1_1')
    #fits.delval(fits_file, 'CD1_2')
    #fits.delval(fits_file, 'CD2_2')
    #fits.delval(fits_file, 'CD2_1')
    pc1_1 = fits.getval(fits_file, 'PC1_1')
    pc1_2 = fits.getval(fits_file, 'PC1_2')
    pc2_1 = fits.getval(fits_file, 'PC2_1')
    pc2_2 = fits.getval(fits_file, 'PC2_2')
    fits.setval(fits_file, 'CD2_2', value=pc2_2, after='CRPIX2')
    fits.setval(fits_file, 'CD2_1', value=pc2_1, after='CRPIX2')
    fits.setval(fits_file, 'CD1_2', value=pc1_2, after='CRPIX2')
    fits.setval(fits_file, 'CD1_1', value=pc1_1, after='CRPIX2')
    fits.delval(fits_file, 'PC1_1')
    fits.delval(fits_file, 'PC1_2')
    fits.delval(fits_file, 'PC2_1')
    fits.delval(fits_file, 'PC2_2')

    #CRPIX alteration
    if tel_angle_corr == -48 and rot_val == 48:
        if "C3" in fits_file:
            print('CRPIX fine shift for C3 implemented!')
            cr1 = fits.getval(fits_file, 'CRPIX1')
            cr2 = fits.getval(fits_file, 'CRPIX2')
            fits.setval(fits_file, 'CRPIX1', value=cr1+33)
            fits.setval(fits_file, 'CRPIX2', value=cr2-25)
        elif "C2" in fits_file:
            print('CRPIX fine shift for C2 implemented!')
            cr1 = fits.getval(fits_file, 'CRPIX1')
            cr2 = fits.getval(fits_file, 'CRPIX2')
            fits.setval(fits_file, 'CRPIX1', value=cr1+15)
            fits.setval(fits_file, 'CRPIX2', value=cr2+16)
        elif "C1" in fits_file:
            print('CRPIX fine shift for C1 implemented!')
            cr1 = fits.getval(fits_file, 'CRPIX1')
            cr2 = fits.getval(fits_file, 'CRPIX2')
            fits.setval(fits_file, 'CRPIX1', value=cr1+-5)
            fits.setval(fits_file, 'CRPIX2', value=cr2+46)
        elif "C4" in fits_file:
            print('CRPIX fine shift for C4 implemented!')
            cr1 = fits.getval(fits_file, 'CRPIX1')
            cr2 = fits.getval(fits_file, 'CRPIX2')
            fits.setval(fits_file, 'CRPIX1', value=cr1+12)
            fits.setval(fits_file, 'CRPIX2', value=cr2+4)
    # elif tel_angle_corr == 0:
        # place found CRPIX values here
           # c1=(0,-44)              #Swift_1246989   #c3 and c4 values are not available
           # c2=(-21,-19)
def update_ra_dec_directory(directory, rot_val):
    fits_files = [os.path.join(directory, f) for f in get_files(directory)]
    for fits_file in fits_files:
        update_ra_dec(fits_file, rot_val)

def update_ra_dec_move_directory(input, output, rot_val):
    for f in sorted(os.listdir(input)):
        if f.endswith('.ramp.fits'):
            origpath = os.path.join(input,f)
            fnewname = f.replace('.ramp.fits', '.ramp.new')
            newpath = os.path.join(output,fnewname)
            shutil.copyfile(origpath, newpath)
            update_ra_dec(newpath, rot_val)
            print('%s updated, renamed, and moved!' % fnewname)


def old_main():
    (ra, dec, rot) = sys.argv[1:4]
    C1_sep = offset_dict['1']['separation']
    C1_ang = float(rot) + offset_dict['1']['angle']
    C3_sep = offset_dict['3']['separation']
    C3_ang = float(rot) + offset_dict['3']['angle']
    C4_sep = offset_dict['4']['separation']
    C4_ang = float(rot) + offset_dict['4']['angle']
    coord_C1 = calc_offset(ra, dec, C1_ang*u.deg, C1_sep*u.deg)
    coord_C3 = calc_offset(ra, dec, C3_ang*u.deg, C3_sep*u.deg)
    coord_C4 = calc_offset(ra, dec, C4_ang*u.deg, C4_sep*u.deg)
    print("C1: " + str(coord_C1.to_string("hmsdms")))
    print("C3: " + str(coord_C3.to_string("hmsdms")))
    print("C4: " + str(coord_C4.to_string("hmsdms")))
#%%

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='alternative to astrometry.net, generates initial astrometry on imgs'
                                                 ' using telescope pointing and corner positions - new')
    parser.add_argument('-input', type=str, help='[str] input path or single file for ramp images, put only this field if you '
                                                 'want to generate astrometry w/o changing names or path')
    parser.add_argument('-output', type=str, help='[str] output path for images w/ astrometry (ramp.new), to '
                                                  'be used in pipeline',default=None)
    parser.add_argument('-rot_val', type=float, help='[float] Use if you want to input a custom rotation value '
                                                     '(the default is 48 deg)', default=48)

    args = parser.parse_args()

    if os.path.isdir(args.input):
        if args.output:
            update_ra_dec_move_directory(args.input, args.output, args.rot_val)
        if not args.output:
            update_ra_dec_directory(args.input, args.rot_val)
    elif os.path.isfile(args.input):
        update_ra_dec(args.input, args.rot_val)
