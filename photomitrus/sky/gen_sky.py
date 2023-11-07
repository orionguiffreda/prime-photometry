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

os.chdir('/Users/orion/Desktop/PRIME/prime-photometry')
from photomitrus.settings import filter
#%%
def sigma_clipped(image, sigma=4):
    masked_array = sigma_clip(image, sigma, maxiters=5)
    return masked_array.filled(np.nan)
#%%
def median_filter_masking(image, size=40):
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
#%%
def gen_sky_image(science_data_directory, output_directory, sky_group_size=None):
    image_fnames = [os.path.join(science_data_directory, f) for f in os.listdir(science_data_directory) if f.endswith('.cds.flat.fits')]
    nfiles = len(image_fnames)
    if sky_group_size is None:
        sky_group_size = nfiles
    image_fnames.sort()
    header = fits.getheader(image_fnames[-1])
    filter1 = header.get('FILTER1', 'unknown')
    filter2 = header.get('FILTER2', 'unknown')
    save_name = 'sky.{}-{}.{}-{}.C{}.fits'.format(filter1, filter2, image_fnames[0][-19:-11], image_fnames[-1][-19:-11], image_fnames[0][-10:-9])
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
        images = [sigma_clipped(fits.getdata(f)) for f in group_files]
        # images = [sigma_clipped(image) for image in images]
        median_array.append(np.nanmedian(images, axis=0))
        file_counter += sky_group_size
    sky = np.nanmedian(median_array, axis=0)  # generating median image
    sky = median_filter_masking(sky)  # filling in the all nan slices
    fits.HDUList([fits.PrimaryHDU(header=header, data=sky)]).writeto(save_name, overwrite=True)
#%%
gen_sky_image('/Users/orion/Desktop/PRIME/GRB/H_band/GRB-Open-H2/','/Users/orion/Desktop/PRIME/GRB/H_band/H_sky', 30)

#%%
#NEED TO ADD DIFF DETECTOR SETTINGS, ONLY HAVE FILTER SETTINGS
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generates sky for given filter and dataset')
    parser.add_argument('filter', nargs=1, type=str, metavar='f', help='Filter being utilized')
    parser.add_argument('paths', nargs='*', type=str, metavar='p', help='Put input path / file first and output path last')
    args = parser.parse_args()
    nint,nframes,sky_size,start_index = filter(args.filter[0])
    gen_sky_image(args.files[0],args.files[1],sky_size)