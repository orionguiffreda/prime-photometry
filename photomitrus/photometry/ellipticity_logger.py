import numpy as np
import argparse
from astropy.table import Table
import os
import pandas as pd
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.colors as mcolors



def split_coordinates(data, n_segs, img_size_x=5000, img_size_y=5000):
    """
    Splits the x and y coordinate data into equal areas based on n_segs.

    Parameters:
        data (pd.DataFrame): DataFrame with 'X_IMAGE' and 'Y_IMAGE' columns.
        n_segs (int): Total number of segments to divide the data into.
        img_size_x (int): Width of the image.
        img_size_y (int): Height of the image.

    Returns:
        dict: A dictionary where keys are segment identifiers and values are DataFrames with coordinates.
    """
    # Determine the number of segments along each axis
    segs_x = int(np.sqrt(n_segs * (img_size_x / img_size_y)))
    segs_y = int(np.sqrt(n_segs * (img_size_y / img_size_x)))

    # Ensure at least one segment per axis
    segs_x = max(segs_x, 1)
    segs_y = max(segs_y, 1)

    # Calculate segment boundaries
    x_bounds = np.linspace(0, img_size_x, segs_x + 1, endpoint=True)
    y_bounds = np.linspace(0, img_size_y, segs_y + 1, endpoint=True)

    # Dictionary to store segmented data
    segments = {}

    for i in range(segs_x):
        for j in range(segs_y):
            # Define boundaries for the current segment
            x_min, x_max = x_bounds[i], x_bounds[i + 1]
            y_min, y_max = y_bounds[j], y_bounds[j + 1]

            # Filter coordinates within the current segment
            segment_data = data[
                (data['X_IMAGE'] >= x_min) & (data['X_IMAGE'] < x_max) &
                (data['Y_IMAGE'] >= y_min) & (data['Y_IMAGE'] < y_max)
            ]

            # Store the segment data in the dictionary
            segment_key = (i, j)
            segments[segment_key] = segment_data

    # resort order by x-axis
    sorted_segments = dict(sorted(segments.items(), key=lambda item: item[0][1]))

    return sorted_segments


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
        cat = Table.read(catname, hdu=2)
    # img = fits.open(imgname)
    # hdr = img[0].header
    # img_x = hdr['NAXIS1']
    # img_y = hdr['NAXIS2']

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


def table_construction(master_ells, master_psfs, stackdir, single=False):
    upper_ells = np.hstack((master_ells[0], master_ells[1]))
    lower_ells = np.hstack((master_ells[2], master_ells[3]))
    final_ells = np.vstack((upper_ells, lower_ells))

    upper_psfs = np.hstack((master_psfs[0], master_psfs[1]))
    lower_psfs = np.hstack((master_psfs[2], master_psfs[3]))
    final_psfs = np.vstack((upper_psfs, lower_psfs))

    ell_df = pd.DataFrame(final_ells)
    psf_df = pd.DataFrame(final_psfs)

    split = stackdir.split('/')
    if single:
        fieldname = [f for f in split if 'field' in f or 'GRB' in f]
        if not fieldname:
            if stackdir[-1] == '/':
                fieldname = split[-3]
            else:
                fieldname = split[-2]
        else:
            fieldname = ' '.join(fieldname)

        elltablename = fieldname + '_single_ellip_table.csv'
        psftablename = fieldname + '_single_psf_table.csv'
    else:
        fieldname = [f for f in split if 'field' in f or 'GRB' in f]
        if not fieldname:
            if stackdir[-1] == '/':
                fieldname = split[-4]
            else:
                fieldname = split[-3]
        else:
            fieldname = ' '.join(fieldname)

        elltablename = fieldname + '_ellip_table.csv'
        psftablename = fieldname + '_psf_table.csv'

    elltablepath = os.path.join(stackdir, elltablename)
    psftablepath = os.path.join(stackdir, psftablename)

    ell_df.to_csv(elltablepath)
    psf_df.to_csv(psftablepath)
    print('Ellipticity and PSF table created in stack dir!: %s, %s' % (elltablepath, psftablepath))
    return final_ells, final_psfs, fieldname


def plots(final_ells, final_psfs, fieldname, stackdir, single=False):

    print('Generating 3d bar plots!')

    plt.close('all')
    plt.clf()
    ar_list = [final_ells, final_psfs]
    titles = ['Median Ellipticity across all chips', 'Median PSF Size across all chips']
    z_labels = ['Ellipticity Value', 'PSF Size Value']

    plot_mins = [0, 1.25]
    plot_maxes = [0.4, 2.5]
    z_lims_norm = [0.65, 4.0]
    z_lims_high = z_lims_norm

    fig, axes = plt.subplots(1, 2, figsize=(10, 6), subplot_kw={'projection': '3d'})
    cmap = plt.get_cmap('plasma')

    dx = dy = 1

    for i, (ar, ax) in enumerate(zip(ar_list, axes)):
        xpos, ypos = np.meshgrid(np.arange(ar.shape[1]), np.arange(ar.shape[0]))
        xpos, ypos = xpos.flatten(), ypos.flatten()
        zpos = np.zeros_like(xpos)
        dz = ar.flatten()

        zlim_max = np.max(ar)
        if zlim_max >= plot_maxes[i]:
            zlim = z_lims_high[i]
        else:
            zlim = z_lims_norm[i]

        norm = mcolors.Normalize(vmin=plot_mins[i], vmax=plot_maxes[i])
        colors = cmap(norm(dz))
        ax.bar3d(xpos, ypos, zpos, dx, dy, dz, shade=True, color=colors, edgecolor='black')
        ax.set_zlim(0, zlim)
        ax.xaxis.set_inverted(True)

        m = cm.ScalarMappable(cmap=cmap, norm=norm)
        m.set_array([])
        plt.colorbar(m, ax=ax, fraction=0.035, pad=0.1)

        ax.view_init(elev=50, azim=45)
        ax.set_xticks(np.arange(0, 7))
        ax.set_yticks(np.arange(0, 7))
        ax.set_xlabel('X Axis Nonants')
        ax.set_ylabel('Y Axis Nonants')
        ax.set_zlabel(z_labels[i])
        ax.set_title(titles[i])

    if single:
        fig.suptitle('Ellipticity and PSF Size 3D plots for %s - Single Image' % fieldname)
        savename = os.path.join(stackdir, '%s_single_detector_plots.png' % fieldname)
    else:
        fig.suptitle('Ellipticity and PSF Size 3D plots for %s - Stacked Image' % fieldname)
        savename = os.path.join(stackdir, '%s_detector_plots.png' % fieldname)
    plt.subplots_adjust(wspace=0.5)
    plt.savefig(savename, dpi=300)
    print('%s generated.' % savename)

    plt.close('all')

#%%


def logger(directory, single=False):
    # stackdir = os.path.join(directory, 'stack')
    if single:
        os.chdir(directory)
        chips = ['1','2','3','4']

        imglist = []
        for c in chips:
            subdir = os.path.join(directory, 'C%s_sub' % c)
            sublist = [f for f in sorted(os.listdir(subdir)) if f.endswith('.flat.cat')]
            imgpath = os.path.join(subdir,sublist[0])
            imglist.append(imgpath)

        print('Checking ellipticity & psfs of single processed images in dir...')
        print(imglist)
        master_ell = []
        master_psfs = []
        for catpath in imglist:
            ells, psfs = ellipticity_log(directory,catpath,directory)
            master_ell.append(ells)
            master_psfs.append(psfs)

        if not master_ell:
            raise FileNotFoundError('No matching catalogs found!')
        else:
            final_ells, final_psfs, fieldname = table_construction(master_ell, master_psfs, directory, single)
            plots(final_ells, final_psfs, fieldname, directory, single)
    else:
        stackdir = directory
        os.chdir(stackdir)

        print('Checking ellipticity & psfs of stacked images...')
        stacklist = [f for f in sorted(os.listdir(stackdir)) if f.endswith('.ecsv') and f.startswith('coadd')]
        # imglist = [f for f in sorted(os.listdir(stackdir)) if f.endswith('.fits') and f.startswith('coadd')]
        print(stacklist)
        master_ell = []
        master_psfs = []
        for catpath in stacklist:
            ells, psfs = ellipticity_log(stackdir,catpath,directory)
            master_ell.append(ells)
            master_psfs.append(psfs)

        if not master_ell:
            raise FileNotFoundError('No matching ecsvs found!')
        else:
            final_ells, final_psfs, fieldname = table_construction(master_ell, master_psfs, stackdir, single)
            plots(final_ells, final_psfs, fieldname, stackdir, single)


def main():
    parser = argparse.ArgumentParser(description='Logs ellipticity for a single observation (checks start, middle, and '
                                                 'end of observation along with stack')
    parser.add_argument('-dir', type=str, help='[str] parent path where observation is stored '
                                               '(likely "/stack"), or parent directory of observation (if you want to '
                                               'run on single images)')
    parser.add_argument('-single', action='store_true', help='if you want to run logger on a set of'
                                                             'processed images (from C#_sub dirs) '
                                                             'instead of just stacks, use this')
    args, unknown = parser.parse_known_args()

    logger(args.dir, args.single)


if __name__ == "__main__":
    main()

