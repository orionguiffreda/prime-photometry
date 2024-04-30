"""
Settings for pipeline
    - Directory locations
    - Detector specific settings
    - Filter specific settings
    - Catalog settings
"""
#%% Config File Names
import os
import fnmatch
#base_dir = os.path.dirname(__file__)
#base_dir = os.path.dirname(os.path.abspath('__file__'))

#export PYTHONPATH=.       if fctn below isnt working from command line, cd to photomitrus and use command
def gen_config_file_name(filename):
    base_dir = os.path.dirname(__file__) #os.path.abspath('__file__')
    return os.path.join(base_dir, 'configs', filename)

def gen_mask_file_name(filename):
    base_dir = os.path.dirname(__file__) #os.path.abspath('__file__')
    return os.path.join(base_dir, 'weightmaps', filename)

def gen_pipeline_file_name():
    base_dir = os.path.dirname(os.path.realpath(__file__))
    return base_dir

def gen_mflat_file_name(filter,chip):
    base_dir = os.path.dirname(os.path.realpath(__file__))
    for f in sorted(os.listdir(base_dir+'/mflats/')):
        if fnmatch.fnmatch(f,'*.Open-%s.C%s.fits' % (filter,chip)):
            filename = f
    return os.path.join(base_dir, 'mflats', filename)

#%%
#os.chdir('/mnt/c/PycharmProjects/prime-photometry/photomitrus/')
#test = gen_mflat_file_name('J',1)

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
    astromlist = [i for i in dirnames if i.endswith(a)]
    astrom = ' '.join(astromlist)
    sublist = [i for i in dirnames if i.endswith(sub)]
    sub = ' '.join(sublist)
    skylist = [i for i in dirnames if i.endswith('sky')]
    sky = ' '.join(skylist)
    stacklist = [i for i in dirnames if i.endswith('stack')]
    stack = ' '.join(stacklist)
    return dir + astrom +'/', dir + sky +'/', dir + sub +'/', dir + stack +'/'

#flat field directory creation
def makedirsFF(dir,chip):
    os.chdir(dir)
    FF = 'C%i_FF' % chip
    skyexists = os.path.exists(FF)
    if not skyexists:
        os.mkdir(FF)
        print(dir + FF)
    if skyexists:
        print(dir + FF +' already exists!')
    dirnames = os.listdir('.')
    FFdir = [i for i in dirnames if i.endswith(FF)]
    FFname = ' '.join(FFdir)
    return dir + FFname + '/'


#%%
#astrom, sky, sub, stack = makedirs('/mnt/d/PRIME_photometry_test_files/GRB240205B/J_Band/flats/',1)
#FF = makedirsFF('/mnt/d/PRIME_photometry_test_files/xrf_3_20240318/J_Band/',2)

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
