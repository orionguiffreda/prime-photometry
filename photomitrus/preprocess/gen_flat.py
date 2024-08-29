"""
Creates master flat
"""
import pandas as pd
from astropy.io import fits
import os
import numpy as np
import argparse
from pathlib import Path

#%% getting lists of flats for beginning and end of night
"""
chip = 1
direct = '/mnt/d/PRIME_photometry_test_files/flat_testing/twilight_flats/'
directory = direct + 'C%i/' % chip
log = pd.read_csv('/mnt/d/PRIME_photometry_test_files/flat_testing/twilight_flats/ramp_fit_log_2024-01-14.clean.dat', delimiter=' ')
"""

def flatlists(path, chip):
    import fnmatch
    for file in os.listdir(path):
        if fnmatch.fnmatch(file, '*.clean.dat'):
            logname = file
            print('taken log = ',logname)
    log = pd.read_csv(path+logname, delimiter=' ', on_bad_lines='warn')
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
            if log_start['FILTER1'][0] == 'Z':
                log_filter = list(log_start['FILTER1'][log_start['OBJNAME'] == 'FLAT'])
            else:
                log_filter = list(log_start['FILTER2'][log_start['OBJNAME'] == 'FLAT'])
            print('On this night, flats were taken in %s band, both at the start and end of the night' % log_filter[0])
        if not log_end_flats:
            if log_start['FILTER1'][0] == 'Z':
                log_filter = list(log_start['FILTER1'][log_start['OBJNAME'] == 'FLAT'])
            else:
                log_filter = list(log_start['FILTER2'][log_start['OBJNAME'] == 'FLAT'])
            print('On this night, flats were taken in %s band, just at the start of the night' % log_filter[0])
    elif log_end_flats:
        if log_start['FILTER1'][0] == 'Z':
            log_filter = list(log_start['FILTER1'][log_start['OBJNAME'] == 'FLAT'])
        else:
            log_filter = list(log_end['FILTER2'][log_end['OBJNAME'] == 'FLAT'])
        print('On this night, flats were taken in %s band, just at the end of the night' % log_filter[0])
    else:
        print('No flats found... ')

    if log_start_flats:
        if not log_end_flats:
            start_images_names_1 = log_start_flats[:int(len(log_start_flats)/2)]
            start_images_names_2 = log_start_flats[int(len(log_start_flats)/2):]
            end_images_names_1 = 0
            end_images_names_2 = 0
            #flag = 'start'

    if log_end_flats:
        start_images_names_1 = 0
        start_images_names_2 = 0
        end_images_names_1 = log_end_flats[:int(len(log_end_flats)/2)]
        end_images_names_2 = log_end_flats[int(len(log_end_flats)/2):]
        #flag = 'end'

    if log_start_flats:
        if log_end_flats:
            start_images_names_1 = log_start_flats[:int(len(log_start_flats) / 2)]
            start_images_names_2 = log_start_flats[int(len(log_start_flats) / 2):]
            end_images_names_1 = log_end_flats[:int(len(log_end_flats) / 2)]
            end_images_names_2 = log_end_flats[int(len(log_end_flats) / 2):]
            #flag = 'both'

    return start_images_names_1, start_images_names_2, end_images_names_1, end_images_names_2

#%% getting data for all groups and stacking along 3rd dim
def flatprocessing(direct,start_images_names_1=None,start_images_names_2=None,end_images_names_1=None,end_images_names_2=None):
    print('getting data for groups of flats and stacking...')

    if start_images_names_1:
        image_list_1 = []
        for f in start_images_names_1:
            img = fits.getdata(f)
            image_list_1.append(img)
        start_images_1 = np.stack(image_list_1)
        start_images_1 = start_images_1[:,4:4092,4:4092]

        image_list_2 = []
        for f in start_images_names_2:
            img = fits.getdata(f)
            image_list_2.append(img)
        start_images_2 = np.stack(image_list_2)
        start_images_2 = start_images_2[:,4:4092,4:4092]

        # subtraction, averaging, then normalizing
        print('start of night: subtracting, normalizing, then median combining...')

        start_median = []
        for i, j in zip(start_images_1, start_images_2):
            sub = i - j
            normsub = sub / np.nanmedian(sub)
            start_median.append(normsub)
        start_median_norm = np.nanmedian(np.stack(start_median), axis=0)

    if end_images_names_1:
        image_list_3 = []
        for f in end_images_names_1:
            img = fits.getdata(f)
            image_list_3.append(img)
        end_images_1 = np.stack(image_list_3)
        end_images_1 = end_images_1[:,4:4092,4:4092]   #remove if cropping issue is fixed

        image_list_4 = []
        for f in end_images_names_2:
            img = fits.getdata(f)
            image_list_4.append(img)
        end_images_2 = np.stack(image_list_4)
        end_images_2 = end_images_2[:,4:4092,4:4092]   #remove if cropping issue is fixed

        # subtraction, averaging, then normalizing
        print('end of night: subtracting, normalizing, then median combining...')

        end_median = []
        for i, j in zip(end_images_1, end_images_2):
            sub = i - j
            normsub = sub / np.nanmedian(sub)
            end_median.append(normsub)
        end_median_norm = np.nanmedian(np.stack(end_median), axis=0)

    """
    #    old algo
    start_sub = start_images_1 - start_images_2
    start_mean = np.nanmean(start_sub,axis=0)
    start_mean_norm = start_mean/(np.nanmean(start_mean))

    end_sub = end_images_1 - end_images_2
    end_mean = np.nanmean(end_sub,axis=0)
    end_mean_norm = end_mean/(np.nanmean(end_mean))
    """
    #low end cut off

    #cut_off = 0.45
    #start_median_norm[start_median_norm < cut_off] = np.nanmedian(start_median_norm)
    #end_median_norm[end_median_norm < cut_off] = np.nanmedian(end_median_norm)


    # writing master flat to file
    print('creating dir for mflats and writing...')
    os.chdir(direct)
    isExist = os.path.exists('mflats')
    if not isExist:
        os.mkdir('mflats')

    if start_images_names_1:
        header_start = fits.getheader(start_images_names_1[-1])
        filter1_start = header_start.get('FILTER1', 'unknown')
        filter2_start = header_start.get('FILTER2', 'unknown')
        save_name_start = 'mflat.{}-{}.{}-{}.C{}.fits'.format(filter1_start, filter2_start, start_images_names_1[0][-20:-12],
                                                        start_images_names_2[-1][-20:-12], start_images_names_1[0][-11])
        output_fname_start = os.path.join(direct+'mflats/', save_name_start)
        print(output_fname_start + ' created!')

        fits.HDUList(fits.PrimaryHDU(header=header_start, data=start_median_norm)).writeto(output_fname_start, overwrite=True)

    if end_images_names_1:
        header_end = fits.getheader(end_images_names_1[-1])
        filter1_end = header_end.get('FILTER1', 'unknown')
        filter2_end = header_end.get('FILTER2', 'unknown')
        save_name_end = 'mflat.{}-{}.{}-{}.C{}.fits'.format(filter1_end, filter2_end, end_images_names_1[0][-20:-12],
                                                        end_images_names_2[-1][-20:-12], end_images_names_1[0][-11])
        output_fname_end = os.path.join(direct+'mflats/', save_name_end)
        print(output_fname_end + ' created!')

        fits.HDUList(fits.PrimaryHDU(header=header_end, data=end_median_norm)).writeto(output_fname_end, overwrite=True)
#%%

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generates 2 master flats (1 from start of night & 1 from end) from given twilight flat data')
    parser.add_argument('-dir', type=str, help='[str], directory where folder of flats and appropriate log is stored')
    #parser.add_argument('-log', type=str, help='[str], path to appropriate log file')
    parser.add_argument('-chip', type=int, help='[int], number of detector')
    args = parser.parse_args()

    start_images_names_1, start_images_names_2,end_images_names_1, end_images_names_2 = flatlists(args.dir,args.chip)
    flatprocessing(args.dir,start_images_names_1, start_images_names_2, end_images_names_1, end_images_names_2)




#%% halves
"""
log = pd.read_csv('/mnt/d/PRIME_photometry_test_files/flat_testing/twilight_flats/ramp_fit_log_2024-01-14.clean.dat', delimiter=' ')
log_start = log.iloc[int(len(log)/2):]
log_start_names = list(log_start['filename'][log_start['OBJNAME']=='FLAT'])
log_start1_names = log_start_names[:int(len(log_start_names)/2)]
log_start2_names = log_start_names[int(len(log_start_names)/2):]
#%%
log_start1_flats = []
for f in log_start1_names:
    new = f.replace('C1.fits.ramp', 'C1.ramp.fits')
    full = '/mnt/d/PRIME_photometry_test_files/flat_testing/twilight_flats/C1/' + new
    log_start1_flats.append(full)

log_start2_flats = []
for f in log_start2_names:
    new = f.replace('C1.fits.ramp', 'C1.ramp.fits')
    full = '/mnt/d/PRIME_photometry_test_files/flat_testing/twilight_flats/C1/' + new
    log_start2_flats.append(full)


start_images1_names_1 = log_start1_flats[:int(len(log_start1_flats)/2)]
start_images1_names_2 = log_start1_flats[int(len(log_start1_flats)/2):]

start_images2_names_1 = log_start2_flats[:int(len(log_start2_flats)/2)]
start_images2_names_2 = log_start2_flats[int(len(log_start2_flats)/2):]
#%%

image_list_1 = []
for f in start_images1_names_1:
    img = fits.getdata(f)
    image_list_1.append(img)
start_images_11 = np.stack(image_list_1)
start_images_11 = start_images_11[:,4:4092,10:]

image_list_2 = []
for f in start_images1_names_2:
    img = fits.getdata(f)
    image_list_2.append(img)
start_images_12 = np.stack(image_list_2)
start_images_12 = start_images_12[:,4:4092,10:]

image_list_3 = []
for f in start_images2_names_1:
    img = fits.getdata(f)
    image_list_3.append(img)
start_images_21 = np.stack(image_list_3)
start_images_21 = start_images_21[:,4:4092,10:]   #remove if cropping issue is fixed

image_list_4 = []
for f in start_images2_names_2:
    img = fits.getdata(f)
    image_list_4.append(img)
start_images_22 = np.stack(image_list_4)
start_images_22 = start_images_22[:,4:4092,10:]   #remove if cropping issue is fixed

#subtraction, averaging, then normalizing
print('subtracting, normalizing, then median combining...')

start_median = []
for i,j in zip(start_images_11,start_images_12):
    sub = i - j
    normsub = sub/np.nanmedian(sub)
    start_median.append(normsub)
start_median_norm1 = np.nanmedian(np.stack(start_median),axis=0)

end_median = []
for k,l in zip(start_images_21,start_images_22):
    sub2 = k - l
    normsub2 = sub2/np.nanmedian(sub2)
    end_median.append(normsub2)
start_median_norm2 = np.nanmedian(np.stack(end_median),axis=0)

#low end cut off
cut_off = 0.45
start_median_norm1[start_median_norm1 < cut_off] = np.nanmedian(start_median_norm1)
start_median_norm2[start_median_norm2 < cut_off] = np.nanmedian(start_median_norm2)

# writing master flat to file
print('creating dir for mflats and writing...')
direct = '/mnt/d/PRIME_photometry_test_files/flat_testing/twilight_flats/halves'
os.chdir(direct)
header_start = fits.getheader(start_images1_names_1[-1])
filter1_start = header_start.get('FILTER1', 'unknown')
filter2_start = header_start.get('FILTER2', 'unknown')
save_name_start = 'mflat.{}-{}.{}-{}.C{}.fits'.format(filter1_start, filter2_start, start_images1_names_1[0][-20:-12],
                                                start_images1_names_2[-1][-20:-12], start_images1_names_1[0][-11])
output_fname_start = os.path.join(direct, save_name_start)
print(output_fname_start + ' created!')

fits.HDUList(fits.PrimaryHDU(header=header_start, data=start_median_norm1)).writeto(output_fname_start, overwrite=True)

header_end = fits.getheader(start_images2_names_1[-1])
filter1_end = header_end.get('FILTER1', 'unknown')
filter2_end = header_end.get('FILTER2', 'unknown')
save_name_end = 'mflat.{}-{}.{}-{}.C{}.fits'.format(filter1_end, filter2_end, start_images2_names_1[0][-20:-12],
                                                start_images2_names_2[-1][-20:-12], start_images2_names_1[0][-11])
output_fname_end = os.path.join(direct, save_name_end)
print(output_fname_end + ' created!')

fits.HDUList(fits.PrimaryHDU(header=header_end, data=start_median_norm2)).writeto(output_fname_end, overwrite=True)
"""

