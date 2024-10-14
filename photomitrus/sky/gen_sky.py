"""
Generates sky for a given detector and a given set of data
"""
#%%
import os

import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clip
from scipy import ndimage
from skimage.morphology import disk
from skimage.filters import rank
import argparse
import warnings
from astropy.utils.exceptions import AstropyWarning
import sys
from datetime import datetime as dt
import matplotlib.pyplot as plt
import warnings

warnings.filterwarnings("ignore", category=RuntimeWarning)

# TODO pipeline test
#%%


def sigma_clipped(image, sigma, sky=0):
    if sigma is None:
        return image
    masked_array = sigma_clip(image-sky, sigma, maxiters=5)
    image[masked_array.mask] = np.nan
    return image


def median_filter_masking(image, size=50):
    nan_percent = 100 * np.count_nonzero(np.isnan(image)) / (image.shape[0] * image.shape[1])
    print('start_nan_percentage: {}'.format(nan_percent))
    median_filter_image = ndimage.median_filter(image, size=size)
    image[np.isnan(image)] = median_filter_image[np.isnan(image)]
    nan_percent = 100 * np.count_nonzero(np.isnan(image)) / (image.shape[0] * image.shape[1])
    print('end_nan_percentage: {}'.format(nan_percent))
    return image


def mean_filter_masking(image, size=30):
    footprint = disk(size)
    nan_percent = 100 * np.count_nonzero(np.isnan(image)) / (image.shape[0] * image.shape[1])
    print('start_nan_percentage: {}'.format(nan_percent))
    filter_image = image.copy()
    filter_image[np.isnan(filter_image)] = np.nanmedian(image)
    plt.imsave('filter_img.png',filter_image,vmin=-1,vmax=1)
    # filter_image_norm = np.max(np.abs(filter_image))
    nan_percent = 100 * np.count_nonzero(np.isnan(filter_image)) / (image.shape[0] * image.shape[1])
    print('filter_nan_percentage: {}'.format(nan_percent))
    # mean_filter_image = rank.mean(filter_image / filter_image_norm, footprint=footprint) * filter_image_norm
    mean_filter_image = ndimage.uniform_filter(filter_image, size=size,mode='constant')
    nan_percent = 100 * np.count_nonzero(np.isnan(mean_filter_image)) / (mean_filter_image.shape[0] * mean_filter_image.shape[1])
    print('mean_filter_nan_percentage: {}'.format(nan_percent))
    image[np.isnan(image)] = mean_filter_image[np.isnan(image)]
    nan_percent = 100 * np.count_nonzero(np.isnan(image)) / (image.shape[0] * image.shape[1])
    print('end_nan_percentage: {}'.format(nan_percent))
    return image


def gen_sky_image(science_data_directory,output_directory, sky_group_size=None,sigma=None,nan_thresh=3):
    warnings.simplefilter('ignore', category=AstropyWarning)
    image_fnames = [os.path.join(science_data_directory, f) for f in os.listdir(science_data_directory) if
                    f.endswith('.ramp.new') or f.endswith('.flat.fits')]
    nfiles = len(image_fnames)
    if sky_group_size is None:
        sky_group_size = nfiles
    image_fnames.sort()
    header = fits.getheader(image_fnames[-1])
    filter1 = header.get('FILTER1', 'unknown')
    filter2 = header.get('FILTER2', 'unknown')
    #save_name = 'sky.Open-J.00747455-00747767.C4.fits'
    save_name = 'sky.{}-{}.{}-{}.C{}.fits'.format(filter1, filter2, image_fnames[0][-19:-11],
                                                    image_fnames[-1][-19:-11], image_fnames[0][-10])
    #save_name = os.path.join(output_directory, save_name)
    print(save_name)
    # breaking the files into group to avoid using too much memory
    ngroups = int(nfiles / sky_group_size)
    if ngroups*sky_group_size > nfiles:
        ngroups -= 1
    median_array = []
    file_counter = 0
    for i in range(ngroups):
        # calculating sky frame for group
        group_files = image_fnames[file_counter:file_counter+sky_group_size]
        print(group_files)
        images = [sigma_clipped(fits.getdata(f), sigma) for f in group_files]
        images = [img / np.nanmedian(img) for img in images]
        # images = [sigma_clipped(image) for image in images]
        median_array.append(np.nanmedian(images, axis=0))
        file_counter += sky_group_size
    sky = np.nanmedian(median_array, axis=0)  # generating median image
    nan_percent = 100 * np.count_nonzero(np.isnan(sky)) / (sky.shape[0] * sky.shape[1])
    if nan_percent >= nan_thresh:
        sky = median_filter_masking(sky)  # filling in the all nan slices
    median = np.nanmedian(sky)
    print('filling remaining nan with median value: {}'.format(median))
    sky[np.isnan(sky)] = median
    fits.HDUList([fits.PrimaryHDU(header=header, data=sky)]).writeto(output_directory+save_name, overwrite=True)


def gen_flat_sky_image(science_data_directory,output_directory, sky_group_size=None,sigma=None,nan_thresh = 0):
    warnings.simplefilter('ignore', category=AstropyWarning)
    image_fnames = [os.path.join(science_data_directory, f) for f in os.listdir(science_data_directory) if
                    f.endswith('.ramp.new') or f.endswith('.flat.fits')]
    nfiles = len(image_fnames)
    if sky_group_size is None:
        sky_group_size = nfiles
    image_fnames.sort()
    header = fits.getheader(image_fnames[-1])
    filter1 = header.get('FILTER1', 'unknown')
    filter2 = header.get('FILTER2', 'unknown')
    save_name = 'sky.{}-{}.{}-{}.C{}.fits'.format(filter1, filter2, image_fnames[0][-20:-12],
                                                    image_fnames[-1][-20:-12], image_fnames[0][-11])
    #save_name = 'sky.test.fits'
    #save_name = os.path.join(output_directory, save_name)
    print(save_name)
    # breaking the files into group to avoid using too much memory
    ngroups = int(nfiles / sky_group_size)
    if ngroups*sky_group_size > nfiles:
        ngroups -= 1
    median_array = []
    file_counter = 0
    group_start = dt.now()
    for i in range(ngroups):
        # calculating sky frame for group
        group_files = image_fnames[file_counter:file_counter+sky_group_size]
        print(group_files)
        images = [sigma_clipped(fits.getdata(f), sigma) for f in group_files]
        images = [img / np.nanmedian(img) for img in images]
        # images = [sigma_clipped(image) for image in images]
        median_array.append(np.nanmedian(images, axis=0))
        file_counter += sky_group_size
    sky = np.nanmedian(median_array, axis=0)  # generating median image
    nan_percent = 100 * np.count_nonzero(np.isnan(sky)) / (sky.shape[0] * sky.shape[1])
    if nan_percent >= nan_thresh:
        sky = mean_filter_masking(sky)  # filling in the all nan slices
    median = np.nanmedian(sky)
    print('filling remaining nan with median value: {}'.format(median))
    sky[np.isnan(sky)] = median
    fits.HDUList([fits.PrimaryHDU(header=header, data=sky)]).writeto(output_directory+save_name, overwrite=True)


def gen_mean_flat_sky_image(science_data_directory,output_directory, sky_group_size=None,sigma=None):
    warnings.simplefilter('ignore', category=AstropyWarning)
    image_fnames = [os.path.join(science_data_directory, f) for f in os.listdir(science_data_directory) if
                    f.endswith('.ramp.new') or f.endswith('.flat.fits')]
    nfiles = len(image_fnames)
    if sky_group_size is None:
        sky_group_size = nfiles
    image_fnames.sort()
    header = fits.getheader(image_fnames[-1])
    filter1 = header.get('FILTER1', 'unknown')
    filter2 = header.get('FILTER2', 'unknown')
    save_name = 'sky.{}-{}.{}-{}.C{}.fits'.format(filter1, filter2, image_fnames[0][-20:-12],
                                                    image_fnames[-1][-20:-12], image_fnames[0][-11])
    #save_name = 'sky.test.fits'
    #save_name = os.path.join(output_directory, save_name)
    print(save_name)
    # breaking the files into group to avoid using too much memory
    ngroups = int(nfiles / sky_group_size)
    if ngroups*sky_group_size > nfiles:
        ngroups -= 1
    median_array = []
    file_counter = 0
    for i in range(ngroups):
        # calculating sky frame for group
        group_files = image_fnames[file_counter:file_counter+sky_group_size]
        print(group_files)
        images = [sigma_clipped(fits.getdata(f), sigma) for f in group_files]
        images = [img / np.nanmedian(img) for img in images]
        # images = [sigma_clipped(image) for image in images]
        median_array.append(np.nanmedian(images, axis=0))
        file_counter += sky_group_size
    sky = np.nanmedian(median_array, axis=0)  # generating median image
    sky = mean_filter_masking(sky)  # filling in the all nan slices
    fits.HDUList([fits.PrimaryHDU(header=header, data=sky)]).writeto(output_directory+save_name, overwrite=True)


def sky_gen(in_path, sky_path, sigma, no_flat=False):
    if no_flat:
        gen_sky_image(science_data_directory=in_path,output_directory=sky_path, sky_group_size=None,sigma=sigma)
    else:
        gen_flat_sky_image(science_data_directory=in_path, output_directory=sky_path, sky_group_size=None, sigma=sigma)

#%%


def main():
    parser = argparse.ArgumentParser(description='Generates sky for given filter and dataset')
    parser.add_argument('-no_flat',  action='store_true', help='if you had NOT flat-fielded, use this')
    # parser.add_argument('filter', nargs=1, type=str, metavar='f', help='Filter being utilized (put first)')
    parser.add_argument('-in_path', type=str, help='[str] Input imgs path (usually ramps w/ astrometry)')
    parser.add_argument('-sky_path', type=str, help='[str] output sky path')
    parser.add_argument('-sigma', type=int, help='[int] Sigma value for sigma clipping',default=None)
    args, unknown = parser.parse_known_args()

    sky_gen(args.in_path, args.sky_path, args.sigma, args.no_flat)


if __name__ == "__main__":
    main()