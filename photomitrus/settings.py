"""
Settings for pipeline
    - Directory locations
    - Detector specific settings
    - Filter specific settings
    - Catalog settings
"""
#%% Config File Names
import os
#base_dir = os.path.dirname(__file__)
#base_dir = os.path.dirname(os.path.abspath('__file__'))

#export PYTHONPATH=.       if fctn below isnt working from command line, cd to photomitrus and use command
def gen_config_file_name(filename):
    base_dir = os.path.dirname(__file__) #os.path.abspath('__file__')
    return os.path.join(base_dir, 'configs', filename)

def gen_mask_file_name(filename):
    base_dir = os.path.dirname(__file__) #os.path.abspath('__file__')
    return os.path.join(base_dir, 'weightmaps', filename)

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

#%%
import pandas as pd

#Settings for list of directory+filenames
object = 'field4057'
filter = 'H'
chip = 4

def flist(Object=object, Filter=filter, Chip=chip):
    log = pd.read_csv('/mnt/d/PRIME_photometry_test_files/ramp_fit_log_2023-12-10.clean.dat', delimiter=' ') #reads in csv
    if Filter == 'Z':
        fnames = log['filename'][log['CHIP'] == Chip & log['OBJNAME'].str.contains(str(Object)) & ~log['OBJNAME'].str.contains('test') & log['FILTER1'].str.contains('Z')
        & log['Open'].str.contains(str(Filter))]
        fnames = fnames.tolist()
        dir = ('/mnt/d/PRIME_photometry_test_files/C{}/'.format(Chip))
        fullnames = [dir + x for x in sorted(fnames)]
    else:
        #screening log w/ constraints
        fnames = log['filename'][log['OBJNAME'].str.contains(str(Object)) & ~log['OBJNAME'].str.contains('test') & log['FILTER1'].str.contains('Open')
                 & log['FILTER2'].str.contains(str(Filter)) & log['OBSERVER'].str.contains('NASA')]
        fnames = fnames.tolist()
        #adding path
        dir = ('/mnt/d/PRIME_photometry_test_files/C{}/'.format(Chip))
        fullnames = [dir + x for x in sorted(fnames)]
        #fullnames = fnames
    fullnames = [f.replace('fits.ramp','ramp.fits') for f in fullnames]
    fullnames = [f.replace('C1','C%s' % chip) for f in fullnames]
    if not fullnames:
        print('No files found for specified fields!')
    return fullnames

a = flist()
#%% make directories
import os
from pathlib import Path
def makedirs(dir,chip):
    os.chdir(dir)
    a = 'C%i_astrom' % chip
    sub = 'C%i_sub' % chip
    skyexists = os.path.exists('sky')
    if not skyexists:
        os.mkdir('sky')
        print(dir + 'sky')
    if skyexists:
        print(dir + 'sky already exists!')
    stackexists = os.path.exists('stack')
    if not stackexists:
        os.mkdir('stack')
        print(dir + 'stack')
    if stackexists:
        print(dir + 'stack already exists!')
    aexists = os.path.exists(a)
    if not aexists:
        os.mkdir(a)
        print(dir + a)
    if aexists:
        print(dir + a+' already exists!')
    subexists = os.path.exists(sub)
    if not subexists:
        os.mkdir(sub)
        print(dir + sub)
    if subexists:
        print(dir + sub+' already exists!')
    dirnames = os.listdir('.')
    astromlist = [i for i in dirnames if a in i]
    astrom = ' '.join(astromlist)
    sublist = [i for i in dirnames if sub in i]
    sub = ' '.join(sublist)
    skylist = [i for i in dirnames if 'sky' in i]
    sky = ' '.join(skylist)
    stacklist = [i for i in dirnames if 'stack' in i]
    stack = ' '.join(stacklist)
    return dir + astrom +'/', dir + sky +'/', dir + sub +'/', dir + stack +'/'

#%%
#astrom, sky, sub, stack = makedirs('/mnt/d/PRIME_photometry_test_files/GRB240205B/J_Band/flats/',1)

#%% for multiple object fields
"""
object = ['field4037']
filter = 'J'
chip = 1
def flist(Object=object, Filter=filter, Chip=chip):
    log = pd.read_csv('/mnt/d/PRIME_photometry_test_files/ramp_fit_log_2023-12-10.clean.dat', delimiter=' ') #reads in csv
    if Filter == 'Z':
        fnames = []
        for f in Object:
            fobjnames = log['filename'][log['CHIP'] == Chip & log['OBJNAME'].str.contains(str(Object)) & log['FILTER1'].str.contains('Z')
            & log['Open'].str.contains(str(Filter))]
            fobjnames = fobjnames.tolist()
            fnames.extend(fobjnames)
        dir = ('/mnt/d/PRIME_photometry_test_files/C{}/'.format(Chip))
        fullnames = [dir + x for x in sorted(fnames)]
    else:
        #screening log w/ constraints
        fnames = []
        for f in Object:
            fobjnames = log['filename'][log['CHIP'] == Chip & log['OBJNAME'].str.contains(str(f)) & log['FILTER1'].str.contains('Open')
                 & log['FILTER2'].str.contains(str(Filter))]
            fobjnames = fobjnames.tolist()
            fnames.extend(fobjnames)
            #adding path
        dir = ('/mnt/d/PRIME_photometry_test_files/C{}/'.format(Chip))
        fullnames = [dir + x for x in sorted(fnames)]
    if not fullnames:
        print('No files found for specified fields!')
    return fullnames

a = flist()
"""
