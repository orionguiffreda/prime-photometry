"""
Creates master flat
"""
from pandas import to_datetime
from collections import Counter
from astropy.io import fits
import os
import numpy as np
import argparse
import sys

from photometrus.getfiles import (get_log_file, get_data_files)
from photometrus.settings import gen_pipeline_file_name
from photometrus.master import getchiplist
from photometrus.utils.defaults import PROCESSING_DEFAULTS as defaults

#%% getting lists of flats for beginning and end of night


def initflatlist_check(date, band=None, flat_col=defaults['flat_col']):
    """
    Checks for any flats (all chips) taken in band (if applicable) during date, returns list of flat files
    """

    if flat_col != defaults['flat_col']:
        flat_val = 'Flat'
    else:
        flat_val = 'FLAT'

    if not band:
        try:
            flatlists, m_lists = get_data_files(date=date, **{flat_col: flat_val})
            if any(lst for lst in m_lists):
                print('Missing files!: ', m_lists)
        except FileNotFoundError:
            print('Error fetching data!')
    else:
        if band == 'Z':
            try:
                flatlists, m_lists = get_data_files(date=date, **{flat_col: flat_val}, filter1=band,
                                                        filter2='Open')
                if any(lst for lst in m_lists):
                    print('Missing files!: ', m_lists)
            except FileNotFoundError:
                print('Error fetching data!')
        else:
            try:
                flatlists, m_lists = get_data_files(date=date, **{flat_col: flat_val}, filter1='Open',
                                                        filter2=band)
                if any(lst for lst in m_lists):
                    print('Missing files!: ', m_lists)
            except FileNotFoundError:
                print('Error fetching data!')

    return flatlists


def flatlistdownload(date, chip, band=None):
    os.chdir(gen_pipeline_file_name())

    flatlists = initflatlist_check(date, band)
    flatlist = getchiplist(flatlists, chip)

    if not flatlist:
        flatlists = initflatlist_check(date, band, flat_col='objtype')
        flatlist = getchiplist(flatlists, chip)

    return flatlist


def flatlists(date, flatlist, chip):
    if flatlist is None:
        raise FileNotFoundError('No flat fields found in storage dir for date? Are you sure you have the right date? '
                                'Or were there missing files in the flat generation?')

    datetime = to_datetime(date)
    date = datetime.strftime('%Y-%m-%d')
    log = get_log_file(date)

    log_start = log.iloc[:int(len(log)/2)]
    log_end = log.iloc[int(len(log)/2):]

    if flatlist[0].endswith('.ramp.fits'):
        orig_ext = 'fits.ramp'
        new_ext = 'ramp.fits'
    else:
        orig_ext = ''
        new_ext = ''

    log_start_orig_names = list(log_start['filename'][log_start['OBJNAME']=='FLAT'])
    log_start_names = []
    for file in log_start_orig_names:
        new = file.replace('C1.'+orig_ext, f'C{chip}.'+new_ext)
        log_start_names.append(new)
    log_end_orig_names = list(log_end['filename'][log_end['OBJNAME']=='FLAT'])
    log_end_names = []
    for file in log_end_orig_names:
        new = file.replace('C1.'+orig_ext, f'C{chip}.'+new_ext)
        log_end_names.append(new)

    log_start_flats = [flatlistfile for flatlistfile in flatlist if
                       any(logfile in flatlistfile for logfile in log_start_names)]
    log_end_flats = [flatlistfile for flatlistfile in flatlist if
                     any(logfile in flatlistfile for logfile in log_end_names)]

    # redundancy check to confirm correct band (protection against some files being named w/ diff band)
    flat_hdulist = fits.open(flatlist[0])
    if len(flat_hdulist) > 1:
        flat_ext = 1
    else:
        flat_ext = 0

    if log_start_flats:
        flatpop = log_start_flats
    else:
        flatpop = log_end_flats

    flatlistbands = []
    for flt in flatpop:
        hdr = fits.getheader(flt, ext=flat_ext)
        bnd = hdr['FILTER2']
        flatlistbands.append(bnd)

    all_same = len(set(flatlistbands)) == 1
    if all_same:
        flat_filter = flatlistbands[0]
        start_idx = 0
    else:
        majority_band = Counter(flatlistbands).most_common(1)[0][0]
        start_idx = next((i for i in range(1, len(flatlistbands)) if flatlistbands[i] != flatlistbands[i - 1]), None)
        print(f'Some starting files taken in diff band!, Majority band: {majority_band}, skipping bad indices...')
        flat_filter = majority_band

    if flat_filter == 'Open':
        flat_filter = 'Z'

    if log_start_flats:
        if log_end_flats:
            if log_start['FILTER1'][start_idx] == 'Z':
                log_filter = flat_filter
            else:
                log_filter = flat_filter
            print('On this night, flats were taken in %s band, both at the start and end of the night' % log_filter)
        if not log_end_flats:
            if log_start['FILTER1'][start_idx] == 'Z':
                log_filter = flat_filter
            else:
                log_filter = flat_filter
            print('On this night, flats were taken in %s band, just at the start of the night' % log_filter)
    elif log_end_flats:
        if log_start['FILTER1'][start_idx] == 'Z':
            log_filter = flat_filter
        else:
            log_filter = flat_filter
        print('On this night, flats were taken in %s band, just at the end of the night' % log_filter)
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

    return start_images_names_1, start_images_names_2, end_images_names_1, end_images_names_2, flat_filter

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

    save_name_start = None
    save_name_end = None

    if start_images_names_1:
        if start_images_names_1[0].endswith('.fz'):
            header_start = fits.getheader(start_images_names_1[-1], ext=1)
            filter1_start = header_start.get('FILTER1', 'unknown')
            filter2_start = header_start.get('FILTER2', 'unknown')
            save_name_start = 'mflat.{}-{}.{}-{}.C{}.fits'.format(filter1_start, filter2_start,
                                                                  start_images_names_1[0][-23:-15],
                                                                  start_images_names_2[-1][-23:-15],
                                                                  start_images_names_1[0][-14])
        else:
            header_start = fits.getheader(start_images_names_1[-1])
            filter1_start = header_start.get('FILTER1', 'unknown')
            filter2_start = header_start.get('FILTER2', 'unknown')
            save_name_start = 'mflat.{}-{}.{}-{}.C{}.fits'.format(filter1_start, filter2_start, start_images_names_1[0][-20:-12],
                                                            start_images_names_2[-1][-20:-12], start_images_names_1[0][-11])
        output_fname_start = os.path.join(direct, save_name_start)
        print(output_fname_start + ' created!')

        fits.HDUList(fits.PrimaryHDU(header=header_start, data=start_median_norm)).writeto(output_fname_start, overwrite=True)

    if end_images_names_1:
        if end_images_names_1[0].endswith('.fz'):
            header_end = fits.getheader(end_images_names_1[-1], ext=1)
            filter1_end = header_end.get('FILTER1', 'unknown')
            filter2_end = header_end.get('FILTER2', 'unknown')
            save_name_end = 'mflat.{}-{}.{}-{}.C{}.fits'.format(filter1_end, filter2_end,
                                                                end_images_names_1[0][-23:-15],
                                                                end_images_names_2[-1][-23:-15],
                                                                end_images_names_1[0][-14])
        else:
            header_end = fits.getheader(end_images_names_1[-1])
            filter1_end = header_end.get('FILTER1', 'unknown')
            filter2_end = header_end.get('FILTER2', 'unknown')
            save_name_end = 'mflat.{}-{}.{}-{}.C{}.fits'.format(filter1_end, filter2_end, end_images_names_1[0][-20:-12],
                                                            end_images_names_2[-1][-20:-12], end_images_names_1[0][-11])
        output_fname_end = os.path.join(direct, save_name_end)
        print(output_fname_end + ' created!')

        fits.HDUList(fits.PrimaryHDU(header=header_end, data=end_median_norm)).writeto(output_fname_end, overwrite=True)

    return save_name_start, save_name_end
#%%


def flatgen(directory, date, chip, band=None):
    flatlist = flatlistdownload(date, chip, band)
    if flatlist:
        start_images_names_1, start_images_names_2, end_images_names_1, end_images_names_2, flat_filter = flatlists(
            date, flatlist, chip)
        save_name_start, save_name_end = flatprocessing(directory, start_images_names_1, start_images_names_2,
                                                        end_images_names_1, end_images_names_2)
    else:
        print('No mflat data taken during this night, halting mflat gen...')
        save_name_start = save_name_end = flat_filter = 0
    return save_name_start, save_name_end, flat_filter


def main():
    parser = argparse.ArgumentParser(description='Generates 2 master flats (1 from start of night & 1 from end) from given twilight flat data')
    parser.add_argument('-dir', type=str, help='[str], directory')
    parser.add_argument('-date', type=str, help='[str], date of observation, in yyyymmdd format')
    parser.add_argument('-chip', type=int, help='[int], which detector number to run on',
                        default=None)
    parser.add_argument('-band', type=str, help='[str], optional, specify filter ex. "J", otherwise it will'
                                                ' try to generate mflats regardless of filter', default=None)
    args, unknown = parser.parse_known_args()

    save_name_start, save_name_end, flat_filter = flatgen(args.dir, args.date, args.chip, args.band)


if __name__ == "__main__":
    main()




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

