import numpy as np
import argparse
from astropy.table import Table
import os


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
    full_mean = np.mean(cat['ELLIPTICITY'])
    print('Mean Ellipticity for whole image = %.3f' % full_mean)

    sec_means = []
    for quad in stack_split:
        ell_arr = stack_split[quad]['ELLIPTICITY']
        quad_mean = np.mean(ell_arr)
        quad_info = (quad, quad_mean)
        sec_means.append(quad_info)

    all_means = [means[1] for means in sec_means]
    std = np.std(all_means)
    print('Stdev = %.4f' % std)
    sig_2 = full_mean + 1 * std
    print('Sigma Limit = %.3f' % sig_2)
    lim_means_sec = []
    lim_means = []
    for mean in sec_means:
        if mean[1] >= sig_2:
            print(f'Section {mean[0]} is above sigma limit!: Mean = {mean[1]:.3f}')
            lim_means.append(mean[1])
            lim_means_sec.append(mean[0])

    if parentdir:
        directory = parentdir
    dirname = directory.split('/')[-2]
    dirlogname = directory + dirname + '_log.txt'
    log = open(dirlogname,'a')
    log.write('\nEllipticity info on %s\n' % catname)
    log.write('Mean Ellipticity for whole image = %.3f\n' % full_mean)
    for sec,mean in zip(lim_means_sec,lim_means):
        log.write(f'Section {sec} is above sigma limit!: Mean = {mean:.3f}\n')
    log.close()

#%%


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Logs ellipticity for a single observation (checks start, middle, and '
                                                 'end of observation along with stack')
    parser.add_argument('-manual', action='store_true', help='Use along w/ -imgname to manually check specific image')
    parser.add_argument('-dir', type=str, help='[str] parent path where observation is stored (C#_sub/ & stack/'
                                               ' should be in here), or path where specific image is held')
    parser.add_argument('-chip', type=int, help='[int] Chip number to check')
    parser.add_argument('-imgname', type=str, help='[str] optional manual field to specify image filename',
                        default=None)
    args = parser.parse_args()

    if args.manual:
        ellipticity_log(args.dir,args.imgname)
    else:
        subdir = args.dir + 'C%i_sub/' % args.chip
        stackdir = args.dir + 'stack/'

        img_paths = images_list(subdir)
        dirname = args.dir.split('/')[-2]
        dirlogname = args.dir + dirname + '_log.txt'
        if os.path.exists(dirlogname):
            os.remove(dirlogname)

        print('Checking ellipticity of processed images...')
        paths = images_list(subdir)
        for f in paths:
            ellipticity_log(subdir,f,args.dir)

        print('Checking ellipticity of stacked image...')
        for f in os.listdir(stackdir):
            if 'C%i' % args.chip in f and f.endswith('.ecsv'):
                try:
                    stackpath = f
                except NameError:
                    break
        ellipticity_log(stackdir,stackpath,args.dir)

