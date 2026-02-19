"""
Generates sky for a given detector and a given set of data
"""
#%%
import os

import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clip
from astropy.modeling import models, fitting
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

from photometrus.utils.defaults import PROCESSING_DEFAULTS as defaults

warnings.filterwarnings("ignore", category=RuntimeWarning)

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
    return save_name


def gen_flat_sky_image(science_data_directory, output_directory, sky_group_size=None, sigma=None, nan_thresh=0):
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
    fits.HDUList(
        [fits.PrimaryHDU(header=header, data=sky)]
    ).writeto(os.path.join(output_directory, save_name), overwrite=True)
    return save_name


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
    fits.HDUList(
        [fits.PrimaryHDU(header=header, data=sky)]).writeto(os.path.join(output_directory, save_name), overwrite=True)


def gen_poly_fit(sky_img_path, poly_deg=defaults['poly_deg']):
    """
    Generate polynomial sky image.

    Parameters
    ----------
    imagedata: 2D numpy array
        FITS image data to construct fit
    poly_deg: int
        Degree of polynomial for sky fit

    Returns
    ----------
    model: 2D numpy array
        Fitted polynomial data
    coeffs: numpy array
        Polynomial coeffs
    """

    try:
        skydata = fits.getdata(sky_img_path)
    except FileNotFoundError or IndexError:
        raise Exception('Cannot get initial sky data for poly sky, is the path correct and or is the fits image intact?')
    skyhdr = fits.getheader(sky_img_path)

    print('Generating polynomial sky model!')
    imagedata_arr = skydata
    ny, nx = imagedata_arr.shape

    x_max = skyhdr['NAXIS1']
    y_max = skyhdr['NAXIS2']

    # make coord grids / flatten
    x = np.arange(nx)
    y = np.arange(ny)
    X, Y = np.meshgrid(x, y)

    # use Chebyshev2D model
    cheb_init = models.Chebyshev2D(x_degree=poly_deg, y_degree=poly_deg, x_domain=[0, x_max], y_domain=[0, y_max])

    # LinearLSQFitter for fitting
    fitter = fitting.LinearLSQFitter()

    # gen fitted model & get params
    print(f'Fitting Chebyshev2D polynomial of degree {poly_deg}...')
    cheb_model = fitter(cheb_init, X, Y, imagedata_arr)
    model = cheb_model(X, Y)

    coeffs = cheb_model.parameters
    print(f" Chebyshev coefficients: {coeffs}")

    # writing new sky image to sky path
    sky_dir = os.path.split(sky_img_path)[0]
    poly_sky_name = os.path.split(sky_img_path)[1].replace('sky.Open', f'cheb.{poly_deg}.sky.Open')
    poly_sky_path = os.path.join(sky_dir, poly_sky_name)
    fits.HDUList(fits.PrimaryHDU(header=skyhdr, data=model)).writeto(poly_sky_path, overwrite=True)

    return model, coeffs, poly_sky_path


def checkplot(output_directory, save_name, poly_deg=None):
    try:
        print('Generating histogram check plot!\n')
        if os.path.isfile(save_name):
            skypath = save_name
        else:
            skypath = os.path.join(output_directory, save_name)
        skyimg = fits.getdata(skypath)
        flat_sky = skyimg.flatten()

        # basic histogram check plot
        plt.figure(figsize=(10, 8))
        plt.hist(flat_sky, bins=100, density=True, edgecolor='black')
        plt.xlabel('Pixel Value')
        plt.ylabel('Normalized Frequency')
        plt.title('Sky Image Histogram')
        plt.savefig('%s.check_plot.png' % skypath, dpi=300)

        if not save_name.startswith('sky.'):
            print(f'Generating self sky - cheb2d sky residual data!')
            skyhdr = fits.getheader(skypath)
            orig_skypath = skypath.replace(f'cheb.{poly_deg}.','')
            orig_skyimg = fits.getdata(orig_skypath)
            sky_resid_img = orig_skyimg - skyimg

            fits.HDUList(fits.PrimaryHDU(header=skyhdr, data=sky_resid_img)).writeto(skypath.replace('cheb.','resid.'),
                                                                                     overwrite=True)

            n_pixels = orig_skyimg.size
            n_params = 5

            resid_std_err = np.nanstd(sky_resid_img) / np.sqrt(n_pixels)

            print('Resid Std Err:', resid_std_err)

            variance = np.nanstd(sky_resid_img)**2
            chi_squared = np.sum(sky_resid_img ** 2 / variance)

            # calculate reduced chi-squared
            dof = n_pixels - n_params
            reduced_chi_squared = chi_squared / dof
            print('red_chi_sq:', reduced_chi_squared)

            flat_sky_resid = sky_resid_img.flatten()
            # residual hist check plot
            plt.figure(figsize=(10, 8))
            plt.hist(flat_sky_resid, bins=100, edgecolor='black')
            plt.yscale('log')
            plt.xlim(-0.3, 0.3)
            plt.xlabel('Pixel Value')
            plt.ylabel('Frequency')
            plt.title(f'Self Sky - Cheb2d Sky ({poly_deg}: Residual Histogram')
            plt.savefig('%s.resid_plot.png' % skypath, dpi=300)

    except ValueError as e:
        raise Exception(f'*WARNING* Issue with sky check plot generation!: {e}'
                        f'\nPerhaps a critical issue with the flat?')


def sky_gen(in_path, sky_path, sigma, check_hist=False, poly=False):
    if poly:
        model, coeffs, poly_sky_path = gen_poly_fit(sky_img_path=sky_path, poly_deg=defaults['poly_deg'])
        checkplot(output_directory=sky_path, save_name=poly_sky_path, poly_deg=defaults['poly_deg'])
    else:
        save_name = gen_flat_sky_image(science_data_directory=in_path, output_directory=sky_path, sky_group_size=None,
                                       sigma=sigma)
        if check_hist:
            checkplot(output_directory=sky_path, save_name=save_name)

#%%


def main():
    parser = argparse.ArgumentParser(description='Generates sky for given filter and dataset')
    parser.add_argument('-poly', action='store_true', help='Optional arg, use to generate polynomial sky '
                                                           'image from initial sky.')
    parser.add_argument('-check_hist', action='store_true', help='Optional arg, use to generate histogram check plot '
                                                                 'for sky image.')
    # parser.add_argument('filter', nargs=1, type=str, metavar='f', help='Filter being utilized (put first)')
    parser.add_argument('-in_path', type=str, help='[str] Input imgs path (usually ramps w/ astrometry), not necessary '
                                                   'w/ -poly')
    parser.add_argument('-sky_path', type=str, help='[str] output sky path')
    parser.add_argument('-sigma', type=int, help='[int] Sigma value for sigma clipping',default=defaults['sigma'])
    args, unknown = parser.parse_known_args()

    sky_gen(args.in_path, args.sky_path, args.sigma, args.check_hist, args.poly)


if __name__ == "__main__":
    main()
