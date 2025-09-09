from astropy.io import fits
import numpy as np


def ramp_frame_reduce(images, saturation=50000):
    initial_shape = images.shape
    new_shape = (initial_shape[0], initial_shape[1]*initial_shape[2])
    final_shape = (initial_shape[1], initial_shape[2])
    images = np.reshape(images, new_shape)
    images = images.astype(float)
    # saturation = load_settings()['SATURATE']
    slope_flat = np.nan * np.ones(new_shape[-1])
    times = np.arange(initial_shape[0])
    for time in np.flip(times[1:]):
        # print(time)
        not_saturated = images[time] < saturation
        complete = np.isfinite(slope_flat)
        # print(np.sum(complete)/np.prod(complete.shape))
        calc_slopes = not_saturated & np.invert(complete)
        slope_images = np.asarray([images[t][calc_slopes] for t in range(time+1)])
        fit = np.polyfit(times[:time+1].astype(float), slope_images, 1)
        slope_flat[calc_slopes] = fit[-2]
    slope = np.reshape(slope_flat, final_shape)
    slope = slope.astype(np.float32)
    times = times.astype(np.float32)
    return slope


def reduce_image_from_file_list(file_list, hdu_ext):
    images = [fits.getdata(f, ext=hdu_ext) for f in file_list]
    header = fits.getheader(file_list[-1], ext=hdu_ext)
    return header, ramp_frame_reduce(np.asarray(images))

