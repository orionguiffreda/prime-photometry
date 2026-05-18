import os
from datetime import datetime
import numpy as np
import pandas as pd
from astropy import units as u
from astropy_healpix import HEALPix
from astropy.coordinates import Angle
from filelock import FileLock, Timeout
from datetime import datetime as dt
import argparse

#%%

# Settings
hpx_level = 6


def open_bulk_file(filename, timeout=30):

    lock_path = filename + ".lock"
    lock = FileLock(lock_path, timeout=timeout)

    try:
        with lock:
            bulk_file = open(filename, "w")
            return bulk_file
    except Timeout:
        raise Timeout(f"Could not acquire lock for '{filename}' within {timeout}s")


def ref_file_setup(output_dir):

    DR3            = False            # Set it to False to select EDR3
    target_table   = 'gaia_source'   # Alternative values: 'Astrophysical_parameters/astrophysical_parameters/', etc

    # Download basic parameters ============

    if not os.path.isdir(f'{output_dir}'):
        os.system(f'mkdir {output_dir}')

    output_file    = os.path.join(output_dir, 'bulk_download_files.txt')

    print(f'Files will be downloaded to: {output_dir}')

    if DR3:
        gaia_dr_flag = 'DR3'
    else:
        gaia_dr_flag = 'EDR3'

    print('='*120)
    print(f'Preparing selection of Gaia {gaia_dr_flag}: ""{target_table}" files')
    print('='*120)
    print('')


    url_prefix      = f'http://cdn.gea.esac.esa.int/Gaia/g{gaia_dr_flag.lower()}/{target_table}/'
    md5sum_file_url = url_prefix + '_MD5SUM.txt'
    md5sum_file     = pd.read_csv(md5sum_file_url, header=None, delim_whitespace=True, names=['md5Sum', 'file'])

    if DR3:
        md5sum_file.drop(md5sum_file.tail(1).index,inplace=True) # The last row in the "_MD5SUM.txt" file in the DR3 directories includes the md5Sum value of the _MD5SUM.txt file

    # Extract HEALPix level-8 from file name ======================================
    healpix_8_min  = [int(file[file.find('_')+1:file.rfind('-')])     for file in md5sum_file['file']]
    healpix_8_max  = [int(file[file.rfind('-')+1:file.rfind('.csv')]) for file in md5sum_file['file']]
    reference_file = pd.DataFrame({'file':md5sum_file['file'], 'healpix8_min':healpix_8_min, 'healpix8_max':healpix_8_max}).reset_index(drop=True)

    # Compute HEALPix levels 6,7, and 9 ===========================================
    reference_file['healpix7_min'] = [inp >> 2 for inp in reference_file['healpix8_min']]
    reference_file['healpix7_max'] = [inp >> 2 for inp in reference_file['healpix8_max']]

    reference_file['healpix6_min'] = [inp >> 2 for inp in reference_file['healpix7_min']]
    reference_file['healpix6_max'] = [inp >> 2 for inp in reference_file['healpix7_max']]

    reference_file['healpix9_min'] = [inp << 2       for inp in reference_file['healpix8_min']]
    reference_file['healpix9_max'] = [(inp << 2) + 3 for inp in reference_file['healpix8_max']]

    # Generate reference file =====================================================
    ncols          = ['file', 'healpix6_min', 'healpix6_max', 'healpix7_min', 'healpix7_max', 'healpix8_min', 'healpix8_max', 'healpix9_min', 'healpix9_max']
    reference_file = reference_file[ncols]

    return reference_file, output_file, url_prefix


def write_dl_file_list(lon, lat, radius, reference_file, output_dir, output_file, url_prefix):

    # Compute Healpix indexes associated to the selected circular region =====================================================
    print(f'Computing HEALPix Level {hpx_level} encompasing a Cone Search (Radius, longitude, latitude): '
          f'{round(radius.value,2)} {radius.unit},  {round(lon.value,2)} {lon.unit}, {round(lat.value,2)} {lat.unit}')

    hp             = HEALPix(nside=2**hpx_level, order='nested')
    hp_cone_search = hp.cone_search_lonlat(lon, lat, radius = radius)

    # Write applicable files to output for download =====================================================
    # f = open(output_file, "w")
    f = open_bulk_file(output_file, timeout=30)

    dl_list = []
    for index in reference_file.index:
        row = reference_file.iloc[index]
        hp_min, hp_max = row[f'healpix{hpx_level}_min'], row[f'healpix{hpx_level}_max']
        if np.any(np.logical_and(hp_min <= hp_cone_search, hp_cone_search <= hp_max)):
            bulk_file = url_prefix + row['file'] + '\n'
            cat_file = row['file']
            if cat_file not in sorted(os.listdir(output_dir)):
                dl_list.append(cat_file)
                f.write(bulk_file)
    f.close()

    print(f'A total of {len(dl_list)} files for download were written in {output_file}')
    print(dl_list)


def gaia_file_downloader(output_dir, output_file):
    os.system(
        f'wget -i {output_file} -P {output_dir}/ -q  --show-progress --progress=bar:force 2>&1')

    print('='*120)


def mass_gaia_downloader(
        output_dir,
        grid_file='/home/alex/PycharmProjects/prime-photometry/photometrus/configs/obsable_all_sky_grid.csv',
        field_range=None
):
    start_time = dt.now()

    df = pd.read_csv(grid_file, comment="#")

    if field_range is not None:
        range_start, range_end = field_range
        print(f'Downloading for PRIME field range of: {range_start} - {range_end}\n')
        range_start = range_start-1
        range_end = range_end-1
        df = df[range_start:range_end]

    Field_arr = df['ObjectName']
    RA_arr = df['RA']
    DEC_arr = df['DEC']

    reference_file, output_file, url_prefix = ref_file_setup(output_dir)

    output_dir = os.path.join(output_dir, 'downloads')
    if not os.path.isdir(f'{output_dir}'):
        os.system(f'mkdir {output_dir}')

    for id, pointing_ra, pointing_dec in zip(Field_arr, RA_arr, DEC_arr):
        pointing_ra = Angle(pointing_ra, unit=u.hourangle)
        pointing_ra = pointing_ra.deg * u.deg
        pointing_dec = Angle(pointing_dec, unit=u.deg)
        pointing_dec = pointing_dec.deg * u.deg
        pointing_rad = 55 * u.arcmin

        print('')
        print('=' * 120)
        print(f'Downloading for {id}')
        write_dl_file_list(pointing_ra, pointing_dec, pointing_rad, reference_file, output_dir, output_file, url_prefix)

        gaia_file_downloader(output_dir=output_dir, output_file=output_file)

    end_time = dt.now()
    print('Full time to download range of files:', (end_time - start_time).total_seconds())

#%%
# mass_gaia_downloader(output_dir='/mnt/photometry/local_catalog_db_test/gaia_test/')


def main():
    parser = argparse.ArgumentParser(
        description='Downloads csv.gz files for GAIA EDR3 for all PRIME fields'
    )
    parser.add_argument('--storage_directory', type=str, metavar='SD', help='Specify'
                                                                            ' storage directory where downloaded .csv.gz'
                                                                            ' files are stored.')
    parser.add_argument('--range', type=int, nargs=2, metavar=('START', 'END'))
    args = parser.parse_args()

    field_range = tuple(args.range) if args.range else None

    mass_gaia_downloader(output_dir=args.storage_directory, field_range=field_range)


if __name__ == "__main__":
    main()
