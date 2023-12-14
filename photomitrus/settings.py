"""
Settings for pipeline
    - Directory locations
    - Detector specific settings
    - Filter specific settings
    - Catalog settings
"""
import os

def gen_config_file_name(filename):
    base_dir = os.path.dirname(__file__)
    return os.path.join(base_dir, 'configs', filename)

#%% Filter Settings
def filter(filter):
    if filter == 'H':
        # GRB-Open-H settings
        nint = 90
        sky_size = 30
        # nint = 11
        nframes = 4
        start_index = 406203
        # start_index = 406523
        # raw_directory_format = "{}/Users/orion/Desktop/PRIME/HBAND/grb-H"
        raw_directory_format = r"G:\My Drive\PRIME\prime_data\GRB230815A\raw\C2\GRB-Open-H{}"
        cds_dir = '/Users/orion/Desktop/PRIME/GRB/H_band/GRB-Open-H/'
        return nint,nframes,sky_size,start_index #raw_directory_format,cds_dir
    if filter == 'J':
        nint = 9
        nframes = 8
        sky_size = 15
        # start_index = 406955
        # start_index = 407083
        start_index = 407243
        # raw_directory_format = "{}/Users/orion/Desktop/PRIME/HBAND/grb-H"
        raw_directory_format = r"G:\My Drive\PRIME\prime_data\GRB230815A\raw\C2\GRB-Open-J{}"
        # cds_dir = r'G:\My Drive\PRIME\prime_data\GRB230815A\cds\C2\GRB-Open-J'
        cds_dir = '/Users/jdurbak/Documents/prime_data/GRB230815A/cds/C2/GRB-Open-J'
        return nint, nframes, sky_size, start_index
    if filter == 'Y':
        # GRB-Open-Y settings
        # nint = 45
        # nint = 30
        nint = 22
        nframes = 8
        sky_size = 15
        # start_index = 406579
        # start_index = 406699
        # start_index = 406763
        start_index = 406781
        # raw_directory_format = "{}/Users/orion/Desktop/PRIME/HBAND/grb-H"
        raw_directory_format = r"G:\My Drive\PRIME\prime_data\GRB230815A\raw\C2\GRB-Open-Y{}"
        cds_dir = r'G:\My Drive\PRIME\prime_data\GRB230815A\cds\C2\GRB-Open-Y'
        return nint, nframes, sky_size, start_index
    if filter == 'Z':
        # GRB-Z-Open settings
        nint = 45
        nframes = 8
        sky_size = 15
        # start_index = 407331
        start_index = 407451

        # raw_directory_format = "{}/Users/orion/Desktop/PRIME/HBAND/grb-H"
        raw_directory_format = r"G:\My Drive\PRIME\prime_data\GRB230815A\raw\C2\GRB-Z-Open{}"
        cds_dir = r'G:\My Drive\PRIME\prime_data\GRB230815A\cds\C2\GRB-Z-Open'
        return nint, nframes, sky_size, start_index
    else:
        print('Filters must be specified as H, J, Y or Z!')
