from astropy.io import ascii
from astropy.table import Table
import os
import numpy as np
import matplotlib.pyplot as plt
#%%

def medianelltable(directory, name):
    chips = [1,2,3,4]
    elltable = Table()
    elliplist = []
    for chip in chips:
        files = [f for f in os.listdir(directory) if f.startswith('coadd.') & f.endswith('.C%s.fits.VHS.ecsv' % chip)]
        if not files:
            files = [f for f in os.listdir(directory) if
                     f.startswith('coadd.') & f.endswith('.C%s.fits.2MASS.ecsv' % chip)]
        file = ''.join(files)
        print(file)
        filepath = os.path.join(directory, file)
        catalog = ascii.read(filepath)

        catell = np.nanmedian(catalog['FLUX_RADIUS'])
        print('C%s median radius = %.3f' % (chip, catell))
        elliplist.append(catell)

    elltable['Radius'] = elliplist
    elltable.write('/mnt/d/PRIME_photometry_test_files/focus_testing/%s_ellip_table.ecsv' % name, overwrite=True)
    print('%s table written!' % name)


medianelltable(directory='/mnt/d/PRIME_photometry_test_files/fornax_survey/field2508_20241224/J/stack/', name='Fornax_field2508_1224')

#%% plotting

raddir = '/mnt/d/PRIME_photometry_test_files/focus_testing/'

C1_arr = []
C2_arr = []
C3_arr = []
C4_arr = []

radlist = [t for t in os.listdir(raddir)]
radlist = sorted(radlist, key=lambda x: x[-6])

for table in radlist:
    elltable = ascii.read(os.path.join(raddir,table))
    C1_arr.append(elltable[0][0])
    C2_arr.append(elltable[1][0])
    C3_arr.append(elltable[2][0])
    C4_arr.append(elltable[3][0])

x_arr = np.linspace(1,6,6)

plt.figure(figsize=(10,8))
plt.plot(x_arr, C1_arr, '-o')
plt.plot(x_arr, C2_arr, '-o', color='r', alpha=0.7)
plt.plot(x_arr, C3_arr, '-o', color='g', alpha=0.7)
plt.plot(x_arr, C4_arr, '-o', color='black', alpha=0.7)
plt.ylim(1,3)
plt.legend(['C1','C2','C3','C4'])
plt.grid()
plt.title('Median Source Radii for Various Observations')
plt.xlabel('Observation Date (mm/dd/yy)')
plt.ylabel('Source Radius (arcsec)')
plt.xticks(x_arr, ['12/01/24','11/19/24','8/19/24','8/14/24','4/25/24','4/25/24'])
plt.show()

#%%
import pandas as pd

fields = ['1022','1023','1024','1024','1025','1136','1137','1138','1139','1140','1141']

maxes = []
mins = []
for f in fields:
    data = pd.read_csv('/mnt/photometry/FebGWEvent/field%s_20250212/J/stack/field%s_20250212_psf_table.csv' % (f,f))
    data = data.select_dtypes(include=['float64'])
    max = data.to_numpy().max()
    min = data.to_numpy().min()
    maxes.append(max)
    mins.append(min)

print(np.median(mins))
print(np.median(maxes))

#%%
from photometrus.super_master import get_fields_from_log
tested = get_fields_from_log(date='20250406')
