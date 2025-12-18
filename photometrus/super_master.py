"""
Runs full pipeline (proc. & photom.) on whole night's data
"""

import os
import sys
import argparse
import traceback
# import pandas as pd
from pandas import to_datetime
from collections import defaultdict

from photometrus.getfiles import (get_file_names, get_log_file, log_folder_location)
from photometrus.multi_combo import combo

#%%


def get_most_current_log_date():
    log_folder = log_folder_location.split('/')[:-1]
    log_folder = '/'.join(log_folder)

    log_list = [log for log in os.listdir(log_folder)]
    df_name = sorted(log_list, reverse=True)[0]

    logdate = df_name[-14:-4]
    logdate = logdate.replace('-','')
    return logdate

#%%


def get_fields_from_log(date=None, no_bulge=True):

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
        if not entries:
            continue

        # Sort by timestamp (if not already sorted)
        entries.sort(key=lambda x: x['filename'])  # assuming filenames are timestamped

        # Split entries into groups where INT == 1.0 (i.e., start of new observation)
        observation_blocks = []
        current_block = []

        for entry in entries:
            if float(entry['INT']) == 1.0 and current_block:
                observation_blocks.append(current_block)
                current_block = []
            current_block.append(entry)
        if current_block:
            observation_blocks.append(current_block)

        for block_index, block in enumerate(observation_blocks):
            nint = int(block[0]['NINT'])
            filter_val = block[0]['FILTER2']
            obs = block[0]['OBSERVER']
            nframes = [entry['NFRAMES'] for entry in block]
            files_list = [entry['filename'] for entry in block]
            int_val = len(files_list)
            key_suffix = block_index + 1
            key = f"{objname}_{key_suffix}" if key_suffix > 1 else objname

            structured_data[key] = {
                'NINT': nint,
                'INT': int_val,
                'OBSERVER': obs,
                'BAND': filter_val,
                'FILES': files_list,
                'FRAMES': nframes,
            }

    # pruning out calibration files, i.e. flats, etc.
    structured_target_observations = {k: v for k, v in structured_data.items() if 'FLAT' not in k and '_test' not in k
                                      and 'sky' not in k}
    structured_target_observations = {k: v for k, v in structured_target_observations.items() if 'CALIB'
                                      not in v['OBSERVER']}
    if no_bulge:
        structured_target_observations = {k: v for k, v in structured_target_observations.items() if
                                          'GB' not in k}

    if structured_target_observations:
        print('Running automatically on all applicable fields on', logdate)
        return structured_target_observations
    else:
        raise FileNotFoundError('No applicable fields were taken during %s!' % logdate)


#%%


def full_processing_from_log(init_observations, date, start_field=None):
    print('full_processing_from_log date:', date)
    datetime = to_datetime(date)
    logdate = datetime.strftime('%Y%m%d')

    observations = init_observations

    if start_field:
        if start_field in init_observations:
            print(f'\nInput starting field: {start_field} found in log, '
                  f'starting automated processing from that field on!')
            observations = {}
            keep = False

            for k, v in init_observations.items():
                if k == start_field:
                    keep = True
                if keep:
                    observations[k] = v
        else:
            print(
                f'\n**Input starting field: {start_field} NOT found in log!** '
                f'\nContinuing automated processing as normal.')

    for field in observations:
        print('\n%s' % field)

        # getting band of field
        if observations[field]['BAND'] == 'Open':
            field_band = 'Z'
        else:
            field_band = observations[field]['BAND']

        # getting full filepaths to specific field
        testfieldfiles = [int(file[:8]) for file in observations[field]['FILES']]
        nframes = [n for n in observations[field]['FRAMES']]
        input_ramp_lists = get_file_names(testfieldfiles, nframes)

        # running full processing & photometry on field
        combo_dict = dict(
            target=field, date=logdate, band=field_band, removal=True, auto_mode=True,
            input_ramp_lists=input_ramp_lists
        )
        print('combo command dict:')
        print(combo_dict)

        try:
            combo(**combo_dict)
        except Exception:
            tb = traceback.format_exc()
            print(combo_dict)
            print(tb)


def supermaster(date=None, start_field=None, incl_bulge=False):
    if not date:
        chosendate = get_most_current_log_date()
    else:
        chosendate = date
    if incl_bulge:
        observations = get_fields_from_log(chosendate, no_bulge=False)
    else:
        observations = get_fields_from_log(chosendate, no_bulge=True)
    full_processing_from_log(observations, chosendate, start_field=start_field)


def main():
    parser = argparse.ArgumentParser(description='Super-master script, designed to check the latest log and '
                                                 'process all targets observed over the night')
    parser.add_argument('-date', type=str, help='[str] optional, date of observation in yyyymmdd or similar '
                                                'format, default is the most recent date in the list of logs', default=None)
    parser.add_argument('-start_field', type=str, help='[str] optional, put field name (ex. "field1234") of field'
                                                       ' you want the automated processing to start on', default=None)
    parser.add_argument('-bulge', action='store_true', help='Super-master currently prunes out the bulge '
                                                            'fields by default, as we dont yet support stacking of bulge '
                                                            'fields, if you want to include them anyway, use this flag')
    args, unknown = parser.parse_known_args()  # TODO: get default arguments from defaults dict

    supermaster(args.date, args.start_field, args.bulge)


if __name__ == "__main__":
    main()
