import numpy as np
import argparse
from astropy.table import Table
import os
import pandas as pd


def split_coordinates(data, n_segs, img_size=5000):
    """
    Splits the x and y coordinate data into equal areas based on array_length.

    Parameters:
        data (pd.DataFrame): DataFrame with 'X_IMAGE' and 'Y_IMAGE' columns.
        array_length (int): Length of the array to split the data into segments.
        img_size (int): Size of the image (default is 5000).

    Returns:
        dict: A dictionary where keys are quadrant identifiers and values are DataFrames with coordinates.
    """
    # Calculate the number of segments along one dimension
    segs = int(np.sqrt(n_segs))
    if segs * segs != n_segs:
        segs += 1

    # Calculate segment boundaries for x and y dimensions
    x_bounds = np.linspace(0, img_size, segs + 1, endpoint=True)
    y_bounds = np.linspace(0, img_size, segs + 1, endpoint=True)

    # Create a dictionary to hold the coordinates for each segment
    segments = {}

    for i in range(segs):
        for j in range(segs):
            # Define the boundaries for the current segment
            x_min = x_bounds[i]
            x_max = x_bounds[i + 1]
            y_min = y_bounds[j]
            y_max = y_bounds[j + 1]

            # Filter the coordinates within the current segment
            segment_data = data[
                (data['X_IMAGE'] >= x_min) & (data['X_IMAGE'] < x_max) &
                (data['Y_IMAGE'] >= y_min) & (data['Y_IMAGE'] < y_max)
                ]

            #print('x vals = %s, %s' % (x_min, x_max))
            #print('y vals = %s, %s' % (y_min, y_max))

            # Store the segment data in the dictionary
            segment_key = (i, j)
            segments[segment_key] = segment_data

    return segments


def images_list(directory):
    proc_list = [f for f in sorted(os.listdir(directory)) if f.endswith('.cat')]
    check_imgs = [proc_list[0],proc_list[int(len(proc_list)/2)],proc_list[-1]]
    check_imgs_paths = [os.path.join(directory,f) for f in check_imgs]
    return check_imgs_paths


def ellipticity_log(directory,catname,parentdir=None):
    catpath = os.path.join(directory,catname)
    if catpath.endswith('.ecsv'):
        cat = Table.read(catpath)
    else:
        cat = Table.read(catpath, hdu=2)
    print('\nChecking ellipticity for %s...' % catname)
    stack_split = split_coordinates(cat, n_segs=9)
    full_mean = np.median(cat['ELLIPTICITY'])
    full_psf = np.median(cat['FLUX_RADIUS'])
    print('Median Ellipticity for whole image = %.3f' % full_mean)
    print('Median PSF Size for whole image = %.3f' % full_psf)

    sec_means = []
    sec_psfs = []
    for quad in stack_split:
        ell_arr = stack_split[quad]['ELLIPTICITY']
        psf_arr = stack_split[quad]['FLUX_RADIUS']
        quad_mean = np.median(ell_arr)
        quad_info = (quad, quad_mean)
        quad_mean_psf = np.median(psf_arr)
        quad_info_psf = (quad, quad_mean_psf)
        sec_means.append(quad_info)
        sec_psfs.append(quad_info_psf)

    all_means = [means[1] for means in sec_means]
    all_psfs = [psfs[1] for psfs in sec_psfs]
    std = np.std(all_means)
    print('Ellipticity Stdev = %.4f' % std)
    sig_2 = full_mean + 1 * std
    print('Ellipticity Sigma Limit = %.3f' % sig_2)
    lim_means_sec = []
    lim_means = []

    for mean in sec_means:
        if mean[1] >= sig_2:
            print(f'Nonant {mean[0]} is above ellipticity sigma limit!: Median = {mean[1]:.3f}')
        lim_means.append(mean[1])
        lim_means_sec.append(mean[0])

    std_psf = np.std(all_psfs)
    print('PSF Stdev = %.4f' % std_psf)
    sig_2_psf = full_psf + 1 * std_psf
    print('PSF Sigma Limit = %.3f' % sig_2_psf)
    lim_psfs_sec = []
    lim_psfs = []

    for psf in sec_psfs:
        if psf[1] >= sig_2_psf:
            print(f'Nonant {psf[0]} is above PSF sigma limit!: Median = {psf[1]:.3f}')
        lim_psfs.append(psf[1])
        lim_psfs_sec.append(psf[0])

    ell_reshape = np.array(lim_means).reshape(3, 3)
    psf_reshape = np.array(lim_psfs).reshape(3, 3)
    return ell_reshape, psf_reshape

"""    if parentdir:
        directory = parentdir
    dirname = directory.split('/')[-3]
    dirlogname = directory + dirname + '_log.txt'
    log = open(dirlogname,'a')
    log.write('\nEllipticity info on %s\n' % catname)
    log.write('Mean Ellipticity for whole image = %.3f\n' % full_mean)
    for sec,mean in zip(lim_means_sec,lim_means):
        log.write(f'Section {sec} is above sigma limit!: Mean = {mean:.3f}\n')

    log.write('\nPSF info on %s\n' % catname)
    log.write('Mean PSF Size for whole image = %.3f\n' % full_psf)
    for sec,psf in zip(lim_psfs_sec,lim_psfs):
        log.write(f'Section {sec} is above sigma limit!: Mean = {psf:.3f}\n')
    log.close()"""


def table_construction(master_ells, master_psfs, parentdir):
    upper_ells = np.hstack((master_ells[0], master_ells[1]))
    lower_ells = np.hstack((master_ells[0], master_ells[1]))
    final_ells = np.vstack((upper_ells, lower_ells))

    upper_psfs = np.hstack((master_psfs[0], master_psfs[1]))
    lower_psfs = np.hstack((master_psfs[0], master_psfs[1]))
    final_psfs = np.vstack((upper_psfs, lower_psfs))

    ell_df = pd.DataFrame(final_ells)
    psf_df = pd.DataFrame(final_psfs)

    dirname = parentdir.split('/')[-3]
    elltablename = parentdir + dirname + '_ellip_table.csv'
    psftablename = parentdir + dirname + '_psf_table.csv'

    ell_df.to_csv(elltablename)
    psf_df.to_csv(psftablename)
    print('Ellipticity and PSF table created in parent dir!: %s, %s' % (elltablename, psftablename))

#%%


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Logs ellipticity for a single observation (checks start, middle, and '
                                                 'end of observation along with stack')
    parser.add_argument('-dir', type=str, help='[str] parent path where observation is stored (C#_sub/ & stack/'
                                               ' should be in here), or path where specific image is held')
    parser.add_argument('-imgname', type=str, help='[str] optional manual field to specify image filename',
                        default=None)
    parser.add_argument('-include_sub', action='store_true', help='if you want to run logger on ALL '
                                                                  'processed images instead of just stacks, use this')
    args = parser.parse_args()

    if args.imgname:
        ellipticity_log(args.dir,args.imgname)
    else:
        chips = [1, 2, 3, 4]
        stackdir = args.dir + 'stack/'

        print('Checking ellipticity of stacked images...')
        stacklist = [f for f in sorted(os.listdir(stackdir)) if f.endswith('.ecsv') and f.startswith('coadd')]
        print(stacklist)
        master_ell = []
        master_psfs = []
        for stackpath in stacklist:
            ells, psfs = ellipticity_log(stackdir,stackpath,args.dir)
            master_ell.append(ells)
            master_psfs.append(psfs)

        table_construction(master_ell, master_psfs, args.dir)


