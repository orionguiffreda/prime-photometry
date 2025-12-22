"""
Runs improved astrom on a given file using sextractor and scamp
"""

import os
import sys
import argparse
import subprocess
import sys
import numpy as np
from astropy.io import fits
import threading
from astroquery.vizier import Vizier
import astropy.units as u
from astropy.coordinates import Angle, SkyCoord
from datetime import datetime as dt
from astropy.table import Table, Column

# sys.path.insert(0,'C:\PycharmProjects\prime-photometry\photometrus')
from photometrus.settings import (gen_config_file_name, auto_bulge_detect, ASTROM_QUERY_CATALOGS)
from photometrus.utils.utils import combine_header_and_fits_list
from photometrus.utils.defaults import PROCESSING_DEFAULTS as defaults

#%%


def sextract(imgdir, f, sx, ap):
    pre = os.path.splitext(f)[0]
    ext = os.path.splitext(f)[1]
    com = f'sex {os.path.join(imgdir, pre + ext)} -c {sx} -CATALOG_NAME {pre}.cat -PARAMETERS_NAME {ap}'
    out = subprocess.Popen([com], shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    out.wait()
    # print(pre + ext + ' sextracted!')


def sex(imgdir):
    os.chdir(str(imgdir))
    bulge = auto_bulge_detect(imgdir)
    if bulge:
        print('Using bulge astromatic configs...')
        sx = gen_config_file_name('bulge_new.config')
        ap = gen_config_file_name('tempsource.param')
    else:
        sx = gen_config_file_name('sex.config')
        ap = gen_config_file_name('astrom.param')
    threads = []
    print('Sextracting all images in %s' % imgdir)
    for f in sorted(os.listdir(imgdir)):
        if f.endswith('flat.fits') or f.endswith('.flat.new'):
            thread = threading.Thread(target=sextract, args=(imgdir, f, sx, ap), name=f)
            thread.start()
            threads.append(thread)
    for thread in threads:
        thread.join()
    print('Sextraction of all images complete!')


def sexback(imgdir):
    os.chdir(str(imgdir))
    sx = gen_config_file_name('sex.config')
    ap = gen_config_file_name('astrom.param')
    for f in sorted(os.listdir(str(imgdir))):
        if f.endswith('ramp.fits'):
            pre = os.path.splitext(f)[0]
            ext = os.path.splitext(f)[1]
            com = ["sex ", imgdir + pre + ext, ' -c '+sx, " -CATALOG_NAME " + pre + '.cat', ' -PARAMETERS_NAME '+ap,
                   ' -CHECKIMAGE_NAME '+pre+'.back.fits']
            s0 = ''
            com = s0.join(com)
            out = subprocess.Popen([com], shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            out.wait()
            print(pre + '.back.fits back subbed!')


def scamp(imgdir, distortdeg=None, swarpcat=None, band=None):
    os.chdir(imgdir)
    sc = gen_config_file_name('scamp.conf')
    if band:
        if band == 'Z' or band == 'Y':
            addition = f' -ASTREF_BAND J'
        else:
            addition = f' -ASTREF_BAND {band}'
    else:
        addition = ''
    if swarpcat:
        command = ('scamp %s -c %s' % (swarpcat, sc))
        command = command + addition
        # print('Executing command: %s' % command)
        subprocess.run(command.split(), check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    else:
        img_list = [f for f in sorted(os.listdir(imgdir)) if f.endswith('.cat')]
        img_list = [os.path.join(imgdir, f) for f in img_list]
        img_list = ','.join(img_list)
        if distortdeg:
            command = ('scamp %s -c %s -DISTORT_DEGREES %s' % (img_list, sc, distortdeg))
            command = command + addition
            print('Executing command: %s' % command)
        else:
            command = ('scamp %s -c %s' % (img_list, sc))
            command = command + addition
            print('Executing command: %s' % command)
        subprocess.run(command.split(), check=True)
    # print(pre + ext + ' scamped!')


def improved_scamp(imgdir, band, distortdeg=None, swarpcat=None, bulge=None):
    os.chdir(imgdir)

    fits_list = [f for f in sorted(os.listdir(imgdir)) if f.endswith('.new')]
    firsthdr = fits.getheader(fits_list[0])

    ldac_check = [f for f in sorted(os.listdir(imgdir)) if f.endswith('.ldac')]
    if len(ldac_check) < 1:
        raImage = firsthdr['CRVAL1']
        decImage = firsthdr['CRVAL2']

        width = 48
        if bulge:  # if bulge field, change to galactic coords for query
            print('Galactic bulge field detected!  Adjusting query parameters accordingly...')
            coords = SkyCoord(ra=raImage * u.degree, dec=decImage * u.degree, frame='fk5')
            coords = coords.galactic  # galactic conversion for bulge fields
            chosen_frame = 'galactic'
            frame_long = coords.l.deg
            frame_long_str = 'l = %.4f' % frame_long
            frame_lat = coords.b.deg
            frame_lat_str = 'b = %.4f' % frame_lat
            print('Converting coords to galactic: %s, %s' % (frame_long_str, frame_lat_str))

        else:
            coords = SkyCoord(ra=raImage * u.degree, dec=decImage * u.degree, frame='fk5')
            chosen_frame = 'fk5'
            frame_long = raImage
            frame_long_str = 'RA: %.4f' % frame_long
            frame_lat = decImage
            frame_lat_str = 'DEC: %.4f' % frame_lat

        catalog_dict = ASTROM_QUERY_CATALOGS

        catalogs = []
        if band == 'J' or band == 'H':
            for k, v in catalog_dict.items():
                if 'J' in v[0] or 'H' in v[0]:
                    catalogs.append((k, v[1]))
        elif band == 'Z':
            for k, v in catalog_dict.items():
                if 'Z' in v[0]:
                    catalogs.append((k, v[1]))
        elif band == 'Y':
            for k, v in catalog_dict.items():
                if 'Y' in v[0]:
                    catalogs.append((k, v[1]))
        else:
            print('Only J, H, Y, and Z band are supported!')

        checkwidth = 28
        # current columns
        v = Vizier(columns=['RAJ2000', 'DEJ2000', 'e_RAJ2000', 'e_DEJ2000', 'RPmag', 'e_RPmag',
                            f'{band}mag', f'e_{band}mag', 'Epoch'])
        try:
            result = v.query_region(coords, width=str(checkwidth) + 'm', catalog=[f[1] for f in catalogs])
            # print(result[0])
        except IndexError:
            raise IndexError('Sadly, no current surveys available in current area in %s band' % band)

        keys = result.format_table_list()

        print('%s, %s, Box Width: %s arcmin... '
              '\n%s band initial query resulting in: \n%s' % (
                  frame_long_str, frame_lat_str, checkwidth, band, keys))

        keycheck = result.keys()

        for chosen_survey, catNum in catalogs:
            for k in keycheck:
                if catNum in k:
                    print('%s catalog found!' % k)
                    vhs_table = result[''.join(k)]
                    cols = vhs_table.colnames
                    vhs_band_col = vhs_table[cols[4]]

                    if not np.all(vhs_band_col.mask):
                        print('Survey = %s' % k)

                        print('\nQuerying Vizier %s around %s, %s, boxwidth %.2f arcmin'
                              % (catNum, frame_long_str, frame_lat_str, width))
                        try:
                            v = Vizier(columns=['%s' % cols[0], '%s' % cols[1], '%s' % cols[2], '%s' % cols[3],
                                                '%s' % cols[4], '%s' % cols[5], '%s' % cols[6]],
                                       column_filters={
                                                       "%sFlag" % band.lower(): "<4",
                                                       "%sflags1" % band.lower(): "<16"
                                                       # "%sperrbits" % band: '<=16',
                                                       }, row_limit=-1)
                            astrom_query = v.query_region(SkyCoord(coords, unit=(u.deg, u.deg)),
                                               width=str(width) + 'm'
                                               , catalog=k, cache=False, frame=chosen_frame)
                            if astrom_query and len(astrom_query[0]) > 0:
                                print('Queried source total = ', len(astrom_query[0]))
                                break  # success
                            else:
                                # in case survey provides no sources for some reason, try next available
                                print(f"No sources found in {catNum}, trying fallback if available...")
                        except Exception as e:
                            print('Error in Vizier query.')
                            print(f"Error details: {e}")
                            continue
                    else:
                        print('Query unsuccessful, defaulting to next fallback catalog...')
                else:
                    continue
            else:
                continue
            break

        if swarpcat:
            img_list = swarpcat
            doner_cat = swarpcat
        else:
            img_list = [f for f in sorted(os.listdir(imgdir)) if f.endswith('.cat')]
            img_list_list = [os.path.join(imgdir, f) for f in img_list]
            # for img in img_list_list:
            #     data = Table.read(img, hdu=2)
            #     data = data[(data['FLAGS'] <= 1) &
            #                 (data['FLUX_RADIUS'] * 0.498 > 2) &
            #                 (data['FLUX_MAX'] / data['FLUX_AUTO'] < 0.25)
            #                 ]
            #     hdu_objects = fits.BinTableHDU(data)
            #     hdu_objects.header['EXTNAME'] = 'LDAC_OBJECTS'
            #
            #     primary = fits.PrimaryHDU()
            #     hdul = fits.HDUList([primary, hdu_objects])
            #     hdul.writeto(img, overwrite=True)
            img_list = ','.join(img_list_list)

            doner_cat = img_list_list[0]

        # scamp running
        t = astrom_query[0]
        t[cols[2]] = t[cols[2]].to(u.deg)
        t[cols[3]] = t[cols[3]].to(u.deg)

        # open the working SExtractor cat and get its LDAC_IMHEAD
        doner = fits.open(doner_cat)
        imhead = doner[1]

        # ensure flags exist
        if 'FLAGS' not in t.colnames:
            n = len(t[cols[0]])
            flags_col = Column(np.zeros(n, dtype=np.int32), name='FLAGS', dtype=np.int32)
            gaia_tbl = Table(t).copy()
            gaia_tbl.add_column(flags_col)
        else:
            gaia_tbl = Table(t)

        # Ensure double precision
        for col in gaia_tbl.colnames:
            if col != 'FLAGS':
                gaia_tbl[col] = gaia_tbl[col].astype(np.float32)

        # Add dummy columns
        # dummy_float_cols = ['X_IMAGE', 'Y_IMAGE', 'MAG_AUTO', 'MAGERR_AUTO', 'FLUX_AUTO', 'FLUXERR_AUTO',
        #                     'FLUX_RADIUS', 'SNR_WIN', 'ELONGATION', 'ELLIPTICITY', 'FWHM_IMAGE',
        #                     'XWIN_IMAGE', 'YWIN_IMAGE', 'ERRAWIN_IMAGE', 'ERRBWIN_IMAGE', 'ERRTHETAWIN_IMAGE']
        #
        # for cname in dummy_float_cols:
        #     if cname not in t.colnames:
        #         val = 1.0 if cname in ['ELONGATION', 'FWHM_IMAGE'] else 0.0
        #         gaia_tbl[cname] = Column(np.full(len(gaia_tbl), val, dtype=np.float32))

        # Dummy EPOCH if missing or bad
        if np.all((gaia_tbl['Epoch'] == 0) | np.isnan(gaia_tbl['Epoch'])):
            gaia_tbl['Epoch'] = np.full(len(gaia_tbl), 2016.0)

        # LDAC structure
        hdu_imhead = fits.BinTableHDU(imhead.data, header=imhead.header)
        hdu_imhead.header['EXTNAME'] = 'LDAC_IMHEAD'

        hdu_objects = fits.BinTableHDU(gaia_tbl)
        hdu_objects.header['EXTNAME'] = 'LDAC_OBJECTS'

        primary = fits.PrimaryHDU()
        hdul = fits.HDUList([primary, hdu_imhead, hdu_objects])
        astromcatname = 'astrom_%s_query.ldac' % (band)
        hdul.writeto(astromcatname, overwrite=True)
        print('astrom_%s_query.ldac written!' % (band))

    else:
        print('.ldac catalog found already in directory, taking it and moving on...')

        astromcatname = ldac_check[0]
        cols = Table.read(astromcatname, hdu=2).colnames

        if swarpcat:
            img_list = swarpcat
        else:
            img_list = [os.path.join(imgdir, f) for f in sorted(os.listdir(imgdir)) if f.endswith('.cat')][0]
            img_list = ''.join(img_list)

    sc = gen_config_file_name('scamp_testing.conf')
    command = (f'scamp {img_list} -c {sc} -ASTREF_CATALOG FILE -ASTREFCAT_NAME {astromcatname} '
               f'-ASTREFMAG_LIMITS 8.5,20 -DISTORT_DEGREES {distortdeg} '
               f'-ASTREFCENT_KEYS {cols[0]},{cols[1]} -ASTREFERR_KEYS {cols[2]},{cols[3]} '
               f'-ASTREFMAG_KEY {cols[4]} -ASTREFMAGERR_KEY {cols[5]} -ASTREFOBSDATE_KEY {cols[6]}')
    print(command)
    os.system(command)
    # subprocess.run(command.split(), check=True)


def missfits(imgdir):
    os.chdir(imgdir)
    mc = gen_config_file_name('default.missfits')
    img_list = [f for f in sorted(os.listdir(imgdir)) if f.endswith('flat.fits') or f.endswith('flat.new')]
    img_list = [os.path.join(imgdir, f) for f in img_list]
    # img_list = ' '.join(img_list)
    combine_header_and_fits_list(img_list)
    command = ('missfits -c %s %s' % (mc, img_list))
    print('Executing command: %s' % command)
    # subprocess.run(command.split(), check=True)
    print('Complete!')


#%%

def remove_head(directory):
    fnames = ['.head']
    for f in os.listdir(directory):
        for name in fnames:
            if f.endswith(name):
                path = os.path.join(directory, f)
                try:
                    os.remove(path)
                except Exception as e:
                    print(f"Error removing file: {path} - {e}")


def double_astrom(imgdir, band=None):
    # beginning from where astrom_shift_bulge solved
    start_time = dt.now()
    remove_head(imgdir)
    print('\nSextracting shift-corrected fits files!')
    sex(imgdir)                     # sextract shift-solved fits files
    print('Running SCAMP w/ 2nd order distortion polynomial...')
    scamp(imgdir, distortdeg=2, band=band)                 # scamp shift-solved cat files w/ 2d poly solve
    print('Adding .head files directly to fits hdrs...')
    missfits(imgdir)                            # add 2d-solved scamp hdrs to shifted fits files
    print('Removing 2nd order .head files...')
    remove_head(imgdir)                         # renames .head files to .2d.head to differentiate betw. later scamp run

    print('\nSextracting 2nd order scamp-corrected fits files!')
    sex(imgdir)                     # sextract 2d-solved fits files
    print('Running SCAMP w/ 4th order distortion polynomial...')
    scamp(imgdir, band=band)                               # scamp 2d-solved cat files w/ 4d poly solve
    print('Adding 4th order .head files directly to fits hdrs...')
    missfits(imgdir)                            # replace 2d-solved fits file scamp hdrs w/ 4d soln scamp hdrs
    end_time = dt.now()
    print('\ndouble scamp astrometry time:', (end_time - start_time).total_seconds())


#%%


def astrometry(path, band=None, run_sex=False, run_scamp=False, run_miss=False, double_solve=False, bulge=False,
               improved=False):
    if run_sex:
        sex(path)
    elif run_scamp:
        scamp(path, band=band)
    elif run_miss:
        missfits(path)
    elif double_solve:
        double_astrom(path, band=band)
    elif improved:
        remove_head(path)
        improved_scamp(imgdir=path, band=band, distortdeg=2, bulge=bulge)
        # missfits(path)
    else:
        start_time = dt.now()
        remove_head(path)
        sex(path)
        scamp(path, band=band)
        missfits(path)
        end_time = dt.now()
        print('\nnormal astrometry time:', (end_time - start_time).total_seconds())


def main():
    parser = argparse.ArgumentParser(description='runs sextractor and scamp on input imgs; generates LDAC .cat and .head files')
    parser.add_argument('-sex', action='store_true', help='if you want to run JUST sextractor')
    parser.add_argument('-scamp', action='store_true', help='if you want to run JUST scamp')
    parser.add_argument('-improved', action='store_true', help='if you want to run improved version')
    parser.add_argument('-missfits', action='store_true', help='if you want to run JUST missfits')
    parser.add_argument('-double_solve', action='store_true', help='runs a 2nd order'
                                                            ' then a 4th order poly scamp solution for max accuracy')
    parser.add_argument('-path', type=str, help='[str] Images path (currently just dumps .cat '
                                                                       '& .head files in same path)')
    parser.add_argument('-band', type=str, help='[str] Optional to perhaps increase scamp accuracy, specify'
                                                'band')
    parser.add_argument('-bulge', action='store_true', help='if you want to run on bulge fields',
                        default=defaults['bulge'])
    args, unknown = parser.parse_known_args()

    astrometry(args.path, args.band, args.sex, args.scamp, args.missfits, args.double_solve, args.bulge, args.improved)


if __name__ == "__main__":
    main()
