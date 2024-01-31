"""
Creates master flat
"""
import pandas as pd
from astropy.io import fits
import os
import numpy as np
import argparse

#%% getting lists of flats for beginning and end of night
"""
chip = 1
direct = '/mnt/d/PRIME_photometry_test_files/flat_testing/twilight_flats/'
directory = direct + 'C%i/' % chip
log = pd.read_csv('/mnt/d/PRIME_photometry_test_files/flat_testing/twilight_flats/ramp_fit_log_2024-01-14.clean.dat', delimiter=' ')
"""

def flatlists(path, log, chip):
    print('reading log and splitting flats into appropriate groups...')
    log = pd.read_csv(log,  delimiter=' ')
    direct = path
    directory = direct + 'C%i/' % chip
    log_start = log.iloc[:int(len(log)/2)]
    log_end = log.iloc[int(len(log)/2):]
    log_start_names = list(log_start['filename'][log_start['OBJNAME']=='FLAT'])
    log_end_names = list(log_end['filename'][log_end['OBJNAME']=='FLAT'])

    log_start_flats = []
    for f in log_start_names:
        new = f.replace('C1.fits.ramp', 'C%i.ramp.fits' % chip)
        full = directory + new
        log_start_flats.append(full)

    log_end_flats = []
    for f in log_end_names:
        new = f.replace('C1.fits.ramp', 'C%i.ramp.fits' % chip)
        full = directory + new
        log_end_flats.append(full)

    if log_start_flats:
        if log_end_flats:
            log_filter = list(log_start['FILTER2'][log_start['OBJNAME'] == 'FLAT'])
            print('On this night, flats were taken in %s band, both at the start and end of the night' % log_filter[0])
        if not log_end_flats:
            log_filter = list(log_start['FILTER2'][log_start['OBJNAME'] == 'FLAT'])
            print('On this night, flats were taken in %s band, just at the start of the night' % log_filter[0])
    elif log_end_flats:
        log_filter = list(log_end['FILTER2'][log_end['OBJNAME'] == 'FLAT'])
        print('On this night, flats were taken in %s band, just at the end of the night' % log_filter[0])
    else:
        print('No flats found... ')

    start_images_names_1 = log_start_flats[:int(len(log_start_flats)/2)]
    start_images_names_2 = log_start_flats[int(len(log_start_flats)/2):]

    end_images_names_1 = log_end_flats[:int(len(log_end_flats)/2)]
    end_images_names_2 = log_end_flats[int(len(log_end_flats)/2):]
    return start_images_names_1,start_images_names_2,end_images_names_1,end_images_names_2
#%% getting data for all groups and stacking along 3rd dim
def flatprocessing(start_images_names_1,start_images_names_2,end_images_names_1,end_images_names_2,direct):
    print('getting data for groups of flats and stacking...')
    image_list_1 = []
    for f in start_images_names_1:
        img = fits.getdata(f)
        image_list_1.append(img)
    start_images_1 = np.stack(image_list_1)
    start_images_1 = start_images_1[:,4:4092,10:]

    image_list_2 = []
    for f in start_images_names_2:
        img = fits.getdata(f)
        image_list_2.append(img)
    start_images_2 = np.stack(image_list_2)
    start_images_2 = start_images_2[:,4:4092,10:]

    image_list_3 = []
    for f in end_images_names_1:
        img = fits.getdata(f)
        image_list_3.append(img)
    end_images_1 = np.stack(image_list_3)
    end_images_1 = end_images_1[:,4:4092,10:]

    image_list_4 = []
    for f in end_images_names_2:
        img = fits.getdata(f)
        image_list_4.append(img)
    end_images_2 = np.stack(image_list_4)
    end_images_2 = end_images_2[:,4:4092,10:]

    #subtraction, averaging, then normalizing
    print('subtracting, averaging, and normalizing...')
    start_sub = start_images_1 - start_images_2
    start_mean = np.nanmean(start_sub,axis=0)
    start_mean_norm = start_mean/(np.nanmean(start_mean))

    end_sub = end_images_1 - end_images_2
    end_mean = np.nanmean(end_sub,axis=0)
    end_mean_norm = end_mean/(np.nanmean(end_mean))

    # writing master flat to file
    header_start = fits.getheader(start_images_names_1[-1])
    filter1_start = header_start.get('FILTER1', 'unknown')
    filter2_start = header_start.get('FILTER2', 'unknown')
    save_name_start = 'mflat.{}-{}.{}-{}.C{}.fits'.format(filter1_start, filter2_start, start_images_names_1[0][-20:-12],
                                                    start_images_names_2[-1][-20:-12], start_images_names_1[0][-11])
    output_fname_start = os.path.join(direct, save_name_start)
    print(output_fname_start + ' created!')

    fits.HDUList(fits.PrimaryHDU(header=header_start, data=start_mean_norm)).writeto(output_fname_start, overwrite=True)

    header_end = fits.getheader(end_images_names_1[-1])
    filter1_end = header_end.get('FILTER1', 'unknown')
    filter2_end = header_end.get('FILTER2', 'unknown')
    save_name_end = 'mflat.{}-{}.{}-{}.C{}.fits'.format(filter1_end, filter2_end, end_images_names_1[0][-20:-12],
                                                    end_images_names_2[-1][-20:-12], end_images_names_1[0][-11])
    output_fname_end = os.path.join(direct, save_name_end)
    print(output_fname_end + ' created!')

    fits.HDUList(fits.PrimaryHDU(header=header_end, data=end_mean_norm)).writeto(output_fname_end, overwrite=True)

#%%
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Automation of the backbone of pipeline, currently processes 1 chip at a time')
    parser.add_argument('-dir', type=str, help='[str], directory where folder of flats and appropriate log is stored')
    parser.add_argument('-log', type=str, help='[str], path to appropriate log file')
    parser.add_argument('-chip', type=int, help='[int], number of detector')
    args = parser.parse_args()

    start_images_names_1, start_images_names_2, end_images_names_1, end_images_names_2 = flatlists(args.dir,args.log,args.chip)
    flatprocessing(start_images_names_1, start_images_names_2, end_images_names_1, end_images_names_2)


