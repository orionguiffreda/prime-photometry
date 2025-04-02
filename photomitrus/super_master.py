#%%
import os
import sys
import argparse
import pandas as pd
from pandas import to_datetime
from collections import defaultdict

from photomitrus.getfiles import (get_file_names, get_log_file, log_folder_location)
from photomitrus.multi_combo import combo

#%%


def get_most_current_log_date():
    log_folder = log_folder_location.split('/')[:-1]
    log_folder = '/'.join(log_folder)

    log_list = [log for log in os.listdir(log_folder)]
    df_name = sorted(log_list, reverse=True)[0]

    logdate = df_name[-14:-4]
    logdate = logdate.replace('-','')
    return logdate


def get_fields_from_log(date=None):

    # from given date, pull correct log file
    datetime = to_datetime(date)
    logdate = datetime.strftime('%Y-%m-%d')
    df = get_log_file(logdate)

    # Initialize the storage dictionary
    data_dict = defaultdict(list)

    # Parse file
    midpoint = len(df) // 2
    flat_flag = df.iloc[:midpoint]['OBJNAME'].str.startswith("FLAT").any()

    start_processing = not flat_flag  # If no 'FLAT' in first half, process all rows

    for _, row in df.iterrows():
        entry = row.to_dict()
        objname = entry['OBJNAME']
        if objname.startswith("FLAT"):
            start_processing = True  # Start processing when we hit the first "FLAT"
        if start_processing:
            data_dict[objname].append(entry)

    # Initialize a new dict to store the structured data
    structured_data = {}

    # Process the entries to distribute them into groups based on NINT
    for objname, entries in data_dict.items():
        nint = int(entries[0]['NINT'])
        filter_val = entries[0]['FILTER2']
        obs = entries[0]['OBSERVER']

        key_base = f"{objname}"

        files_list = []
        for i, entry in enumerate(entries):
            files_list.append(entry['filename'])
            if len(files_list) == nint or i == len(entries) - 1:
                int_val = len(files_list)
                key_suffix = (i // nint) + 1
                key = f"{key_base}_{key_suffix}" if key_suffix > 1 else key_base
                structured_data[key] = {
                    'NINT': nint,
                    'INT': int_val,
                    'OBSERVER': obs,
                    'BAND': filter_val,
                    'FILES': files_list.copy()
                }
                files_list.clear()

    # pruning out calibration files, i.e. flats, etc.
    structured_target_observations = {k: v for k, v in structured_data.items() if 'FLAT' not in k and '_test' not in k
                                      and 'sky' not in k}
    structured_target_observations = {k: v for k, v in structured_target_observations.items() if 'CALIB'
                                      not in v['OBSERVER']}

    if structured_target_observations:
        return structured_target_observations
    else:
        sys.exit('No applicable fields were taken during %s!' % logdate)


def full_processing_from_log(observations, date):
    datetime = to_datetime(date)
    logdate = datetime.strftime('%Y%m%d')

    for field in observations:
        print('\n%s' % field)

        # getting band of field
        if observations[field]['BAND'] == 'Open':
            field_band = 'Z'
        else:
            field_band = observations[field]['BAND']

        # getting full filepaths to specific field
        testfieldfiles = [int(file[:8]) for file in observations[field]['FILES']]
        input_ramp_lists = get_file_names(testfieldfiles)

        # running full processing & photometry on field
        combo(target=field, date=logdate, band=field_band, removal=True, auto_mode=True,
              input_ramp_lists=input_ramp_lists)


def supermaster(date=None):
    if not date:
        chosendate = get_most_current_log_date()
    else:
        chosendate = date

    observations = get_fields_from_log(chosendate)
    full_processing_from_log(observations, chosendate)


def main():
    parser = argparse.ArgumentParser(description='Super-master script, designed to check the latest log and '
                                                 'process all targets observed over the night')
    parser.add_argument('-date', type=str, help='[str] optional, date of observation in yyyymmdd or similar '
                                                'format', default=None)
    args, unknown = parser.parse_known_args()

    supermaster(args.date)


if __name__ == "__main__":
    main()