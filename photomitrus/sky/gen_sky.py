"""
Generates sky for a given detector and a given set of data
"""
#%%
import os

import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clip
from scipy import ndimage
import argparse
import sys

sys.path.insert(0,'/Users/orion/Desktop/PRIME/prime-photometry/photomitrus')
from settings import filter
from settings import flist

#a,b,c,d = filter('H')
#%%
def sigma_clipped(image, sigma, sky=0):
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
    med = np.nanmedian(image)
    print('filling remaining nan with median value: {}'.format(med))
    image[np.isnan(image)] = med
    return image

def mean_filter_masking(image, size=50):
    nan_percent = 100 * np.count_nonzero(np.isnan(image)) / (image.shape[0] * image.shape[1])
    print('start_nan_percentage: {}'.format(nan_percent))
    mean_filter_image = ndimage.uniform_filter(image, size=size)
    image[np.isnan(image)] = mean_filter_image[np.isnan(image)]
    nan_percent = 100 * np.count_nonzero(np.isnan(image)) / (image.shape[0] * image.shape[1])
    print('end_nan_percentage: {}'.format(nan_percent))
    med = np.nanmean(image)
    print('filling remaining nan with median value: {}'.format(med))
    image[np.isnan(image)] = med
    return image

def gen_sky_image(science_data_directory,output_directory, sky_group_size=None,sigma=None):
    image_fnames = [os.path.join(science_data_directory, f) for f in os.listdir(science_data_directory) if
                    f.endswith('.ramp.new') or f.endswith('.ramp.fits')]
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
    save_name = os.path.join(output_directory, save_name)
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
        # images = [sigma_clipped(image) for image in images]
        median_array.append(np.nanmedian(images, axis=0))
        file_counter += sky_group_size
    sky = np.nanmedian(median_array, axis=0)  # generating median image
    sky = median_filter_masking(sky)  # filling in the all nan slices
    fits.HDUList([fits.PrimaryHDU(header=header, data=sky)]).writeto(save_name, overwrite=True)

def gen_prev_sky_image(science_data_directory,output_directory, previous_sky=None, sky_group_size=None, sigma=None):
    image_fnames = [os.path.join(science_data_directory, f) for f in os.listdir(science_data_directory) if
                    f.endswith('.cds.flat.fits')]
    nfiles = len(image_fnames)
    if sky_group_size is None:
        sky_group_size = nfiles
    image_fnames.sort()
    header = fits.getheader(image_fnames[-1])
    filter1 = header.get('FILTER1', 'unknown')
    filter2 = header.get('FILTER2', 'unknown')
    save_name = 'sky.{}-{}.{}-{}.C{}.fits'.format(filter1, filter2, image_fnames[0][-24:-16], image_fnames[-1][-24:-16],
                                                  image_fnames[0][-19])
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
        images = [sigma_clipped(fits.getdata(f), sigma, fits.getdata(previous_sky)) for f in group_files]
        # images = [sigma_clipped(image) for image in images]
        median_array.append(np.nanmedian(images, axis=0))
        file_counter += sky_group_size
    sky = np.nanmedian(median_array, axis=0)  # generating median image
    sky = median_filter_masking(sky)  # filling in the all nan slices
    fits.HDUList([fits.PrimaryHDU(header=header, data=sky)]).writeto(output_directory+save_name, overwrite=True)

def gen_prev_meansky_image(science_data_directory, output_directory, previous_sky=0, sky_group_size=None, sigma=2.5):
    image_fnames = [os.path.join(science_data_directory, f) for f in os.listdir(science_data_directory) if f.endswith('.cds.flat.fits')]
    nfiles = len(image_fnames)
    if sky_group_size is None:
        sky_group_size = nfiles
    image_fnames.sort()
    header = fits.getheader(image_fnames[-1])
    filter1 = header.get('FILTER1', 'unknown')
    filter2 = header.get('FILTER2', 'unknown')
    save_name = 'sky.{}-{}.{}-{}.C{}.fits'.format(filter1, filter2, image_fnames[0][-24:-16], image_fnames[-1][-24:-16],
                                                  image_fnames[0][-19])
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
        images = [sigma_clipped(fits.getdata(f), sigma, fits.getdata(previous_sky)) for f in group_files]
        # images = [sigma_clipped(image) for image in images]
        median_array.append(np.nanmean(images, axis=0))
        file_counter += sky_group_size
    sky = np.nanmean(median_array, axis=0)  # generating median image
    sky = mean_filter_masking(sky)  # filling in the all nan slices
    fits.HDUList([fits.PrimaryHDU(header=header, data=sky)]).writeto(output_directory+save_name, overwrite=True)

#%%
#NEED TO ADD DIFF DETECTOR SETTINGS, ONLY HAVE FILTER SETTINGS
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generates sky for given filter and dataset')
    parser.add_argument('-prev_sky',  action='store_true', help='if including prev sky in sky gen, include'
                                                                ' this optional argument')
    #parser.add_argument('filter', nargs=1, type=str, metavar='f', help='Filter being utilized (put first)')
    parser.add_argument('paths', nargs='*', type=str, metavar='p', help='Put input path, output path,'
                                                                        'then (if applicable) put prev sky path')
    parser.add_argument('sigma', nargs=1, type=float, metavar='s', help='Sigma value for sigma clipping'
                                                                      ', put after paths')
    args = parser.parse_args()
    #nint,nframes,sky_size,start_index = filter(args.filter[0])
    if args.prev_sky:
        gen_prev_sky_image(science_data_directory=args.paths[0],output_directory=args.paths[1],previous_sky=args.paths[2],
                      sky_group_size=None,sigma=args.sigma[0])
    else:
        gen_sky_image(science_data_directory=args.paths[0],output_directory=args.paths[1], sky_group_size=None,sigma=args.sigma[0])

    #%%

"""    gen_prev_sky_image(science_data_directory='/Users/orion/Desktop/PRIME/GRB/H_band/H_astrom_flat/',output_directory='/Users/orion/Desktop/PRIME/GRB/H_band/H_sky/',
                       previous_sky='/Users/orion/Desktop/PRIME/GRB/H_band/H_sky/sky.2.5-3.fits', sigma=2.5)
"""