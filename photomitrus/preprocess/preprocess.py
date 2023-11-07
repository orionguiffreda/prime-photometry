"""
Dark subtraction, flat fielding and rough astrometry
"""
#%%
import os

import numpy as np
from astropy.io import fits

# common settings
cameras = range(2, 3)
root_folder = ['', '', '', '']
data_hdu = 0
readout_channels = 33
saturation = 50000
active_pixel_slice_start_y = active_pixel_slice_start_x = 4
active_pixel_slice_end_y = active_pixel_slice_end_x = 4092
file_format = "{:08d}C{}.fits"

#%%
dark_ramp_fname = r"G:\My Drive\PRIME\prime_data\GRB230815A\mramp.Dark-Y.408423-408815.C2.fits"
with fits.open(dark_ramp_fname) as dark_ramp:
    master_dark_frame = dark_ramp[nframes-1].data  # choosing the corresponding dark frame to the exposure frame

#%%
def get_ref_pixel_corrected_frame_and_header(filename):
    hdu = fits.open(filename)[data_hdu]
    n_saturated_pixels = np.sum(hdu.data > saturation)  # calculating the number of saturated pixels
    percent_saturated_pixels = 100 * n_saturated_pixels / (hdu.data.shape[0] * hdu.data.shape[1])  # getting percentage
    print("{} saturation in {}".format(percent_saturated_pixels, os.path.basename(filename)))
    frame = ref_pixel_correct_frame(hdu.data)
    return hdu.header, frame

#%%
def ref_pixel_correct_frame(frame):
    ref_pixels_0 = frame[2:4, :]  # acquiring the desired reference pixels
    ref_pixels_1 = frame[-4:-2, :]  # acquiring the desired reference pixels
    ref_pixels = np.vstack((ref_pixels_0, ref_pixels_1))  # stacking reference pixels
    channel_width = int(frame.shape[1] / readout_channels)  # getting the width of each readout channel
    channel_starts = range(0, frame.shape[1], channel_width)  # getting the starting pixel value of each channel
    subtraction_row = np.zeros(frame.shape[1])  # initializing row array
    for channel_start in channel_starts:
        # finding the median reference pixel value for each channel and setting the subtraction reference subtraction for each channel based on that
        subtraction_row[channel_start:channel_start+channel_width] = \
            np.median(ref_pixels[:,channel_start:channel_start+channel_width])
    subtraction_array = np.zeros(frame.shape)
    subtraction_array[:] = subtraction_row  # creating subtraction array
    frame = frame.astype('int32') - subtraction_array
    return frame[active_pixel_slice_start_y:active_pixel_slice_end_y, active_pixel_slice_start_x:active_pixel_slice_end_x]  # removing reference pixels from final array

#%%
def gen_CDS_image(camera, ramp_start_index, integrations, frames_per_ramp, dark_frame=None):
    raw_directory = raw_directory_format.format(root_folder[camera-1])
    # raw_directory = r'G:\My Drive\PRIME\prime_data\GRB230815A\raw\C2\grb-H'
    ramp_starts = range(ramp_start_index, ramp_start_index + integrations * frames_per_ramp, frames_per_ramp)
    # filter1 = 'unknown'
    # filter2 = 'unknown'
    for ramp_start in ramp_starts:
        hdulist = []
        pedestal_fname = os.path.join(raw_directory, file_format.format(ramp_start, camera))
        pedestal_header, pedestal_frame = get_ref_pixel_corrected_frame_and_header(pedestal_fname)
        exposure_fname = os.path.join(raw_directory, file_format.format(ramp_start+frames_per_ramp-1, camera))
        exposure_header, exposure_frame = get_ref_pixel_corrected_frame_and_header(exposure_fname)
        cds_frame = exposure_frame-pedestal_frame
        if dark_frame is not None:
            cds_frame = cds_frame - dark_frame
        hdulist.append(fits.PrimaryHDU(header=exposure_header, data=cds_frame))
        saved_name = os.path.join(raw_directory, 'cds', file_format.format(ramp_start, camera).replace('.fits', '.cds.fits'))
        if not os.path.isdir(os.path.dirname(saved_name)):
            os.makedirs(os.path.dirname(saved_name))
        fits.HDUList(hdulist).writeto(saved_name,overwrite=True)

#%% flat fielding
def flatfield(directory,flat):
    flatdata = fits.getdata(flat)
    flat = flatdata[4:4 + 4088, 4:4 + 4088]
    for filename in sorted(os.listdir(directory)):
        if filename.endswith('.fits'):
            f = filename.replace('.fits','.flat.fits')
            with fits.open(directory+filename) as hdu:
                data = hdu[0].data
                hdr = hdu[0].header
                flatfield = data / flat
            fits.writeto(directory+f,flatfield,hdr,overwrite=True)
            print(directory+f+' flat fielded!')

#%%
flatfield('/Users/orion/Desktop/PRIME/GRB/H_band/H/','/Users/orion/Desktop/PRIME/GRB/H_band/H_flat/H.flat.C2.fits')