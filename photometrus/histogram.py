import os
import numpy as np
from astropy.io import fits
import os
import matplotlib.pyplot as plt
from astropy.io import ascii
import fnmatch
#%%
sci = fits.getdata('/mnt/d/PRIME_photometry_test_files/H_Band/GRBstacks/coadd.Open-H.00747787-00748103.C4.fits')
comp = sci
flatcomp = comp.flatten()
#flatcomp = flatcomp[~np.isinf(flatcomp)]
flatcomp = np.array([x for x in flatcomp if x<=1E4])
flatcomp = np.array([x for x in flatcomp if x>=-1E4])
#%%
start = fits.getdata('/mnt/d/PRIME_photometry_test_files/flat_testing/twilight_flats/mflats/mflat.Open-J.00887641-00887997.C4.fits')
end = fits.getdata('/mnt/d/PRIME_photometry_test_files/flat_testing/twilight_flats/mflats/mflat.Open-J.00890921-00891357.C4.fits')
diff = start - end
#%%
diff = fits.getdata('/mnt/d/PRIME_photometry_test_files/GRB_Followup/GRB_Followup_flat_tests/flatstack/coadd.Open-J.00888557-00888869.C4.fits')
#diff = fits.getdata('/mnt/d/PRIME_photometry_test_files/GRB_Followup/GRBstacks/coadd.Open-J.00888557-00888869.C4.fits')
diff_flat = diff.flatten()
#flatcomp = flatcomp[~np.isinf(flatcomp)]
diff_flat = np.array([x for x in diff_flat if x<=1E4])
diff_flat = np.array([x for x in diff_flat if x>=-1E4])
#%%
plt.hist(diff_flat,bins = 500,density=True,edgecolor='black')
#plt.ylim(1E-6,2)
plt.xlim(-10,10)
plt.yscale('log')
plt.xlabel('Value')
plt.ylabel('Probability')
plt.title('C4 FF & stacked image histogram')
plt.show()

#%%
import pandas as pd
log = pd.read_csv('/mnt/d/PRIME_photometry_test_files/flat_testing/twilight_flats/H_Band/ramp_fit_log_2024-02-04.clean.dat', delimiter=' ', on_bad_lines='warn')

#%%
import shutil

def move(names,i,o):
    for f in names:
        shutil.copyfile(i+f,o+f)
        print(f+' moved!')

#move(names=flist(),i='/mnt/d/PRIME_photometry_test_files/C4/',o='/mnt/d/PRIME_photometry_test_files/C4_GRBramp/')
#%% Read LDAC tables
def get_table_from_ldac(filename, frame=1):
    """
    Load an astropy table from a fits_ldac by frame (Since the ldac format has column
    info for odd tables, giving it twce as many tables as a regular fits BinTableHDU,
    match the frame of a table to its corresponding frame in the ldac file).

    Parameters
    ----------
    filename: str
        Name of the file to open
    frame: int
        Number of the frame in a regular fits file
    """
    from astropy.table import Table
    if frame > 0:
        frame = frame * 2
    tbl = Table.read(filename, hdu=frame)
    return tbl

#%%
from astropy.io import ascii
import matplotlib.pyplot as plt
from astropy.stats import sigma_clipped_stats

psfsourceTable = ascii.read('/mnt/d/PRIME_photometry_test_files/GRB-1-23/J_Band/stack/coadd.Open-J.00924445-00925157.C1.fits.mag.ecsv')
print(psfsourceTable)
#%%
data = fits.getdata('/mnt/d/PRIME_photometry_test_files/GRB-1-23/J_Band/stack/coadd.Open-J.00924445-00925157.C1.fits')

mean, median, sigma = sigma_clipped_stats(data)
fig = plt.figure(figsize=(10, 10))
ax = fig.gca()
plt.imshow(data, vmin=median - 3 * sigma, vmax=median + 3 * sigma)
circles = [plt.Circle((psfsourceTable['X_IMAGE'][i], psfsourceTable['Y_IMAGE'][i]), radius=5, edgecolor='r', facecolor='None') for i in
           range(len(psfsourceTable['X_IMAGE']))]
for c in circles:
    ax.add_artist(c)

plt.show()
#%%
from astropy.io import ascii
from astropy.table import Column
from astropy.table import Table
csv = ascii.read('/mnt/d/PRIME_photometry_test_files/GRB240205B/J_Band/stack/coadd.Open-J.00943667-00943979.C3.fits.2MASS.ecsv')
mag = csv[0]['JMAG_PSF']
data = [1,2,3,4,5,6]
grbdata = Table()
grbdata['ra']= np.array([mag])
#RA = Column(data[0], unit=u.deg, description='GRB RA')

#%% Stacking and averaging for hot pix correction
import os
from astropy.io import fits
direct = '/mnt/d/PRIME_photometry_test_files/GRB_Followup/C3/'
stack = []
for f in os.listdir(direct):
    img = fits.getdata(direct+f)
    img = img[4:4092,4:4092]
    stack.append(img)
stackarray = np.stack(stack)
avgramp = np.nanmedian(stackarray,axis=0)
#%%
directmap = '/mnt/d/PRIME_photometry_test_files/Bad_pix_map/maps/'
fits.writeto(directmap+'avgrampc3.fits',avgramp,header=None,overwrite=True)
#%% rename
path = '/mnt/d/PRIME_photometry_test_files/EP240417a/J_Band/C3_astrom/'
for f in os.listdir(path):
    if f.endswith('ramp.new'):
        fnewname = f.replace('ramp.new','.ramp.new')
        os.rename(path+f,path+fnewname)
        print('renamed!')

#%%
from astropy.coordinates import SkyCoord
import astropy.units as u
def GRB(ra,dec,imageName,survey,filter,thresh):
    from astropy.table import Table
    GRBcoords = SkyCoord(ra=[ra], dec=[dec], frame='icrs', unit='degree')

    mag_ecsvname = '%s.%s.ecsv'  % (imageName,survey)
    mag_ecsvtable = ascii.read(mag_ecsvname)
    mag_ecsvcleanSources = mag_ecsvtable[(mag_ecsvtable['FLAGS'] == 0) & (mag_ecsvtable['FLAGS_MODEL'] == 0)]
    mag_ecsvsourceCatCoords = SkyCoord(ra=mag_ecsvcleanSources['ALPHA_J2000'], dec=mag_ecsvcleanSources['DELTA_J2000'], frame='icrs',
                                  unit='degree')

    photoDistThresh = thresh
    idx_GRB, idx_GRBcleanpsf, d2d, d3d = mag_ecsvsourceCatCoords.search_around_sky(GRBcoords,
                                                                                 photoDistThresh * u.arcsec)
    print('idx size = %d' % len(idx_GRBcleanpsf))
    if len(idx_GRBcleanpsf) == 1:
        print('GRB source at inputted coords %s and %s found!' % (ra,dec))

        grb_mag = mag_ecsvcleanSources[idx_GRBcleanpsf]['%sMAG_PSF' % filter][0]
        grb_magerr = mag_ecsvcleanSources[idx_GRBcleanpsf]['e_%sMAG_PSF' % filter][0]

        grb_ra = mag_ecsvcleanSources[idx_GRBcleanpsf]['ALPHA_J2000'][0]
        grb_dec = mag_ecsvcleanSources[idx_GRBcleanpsf]['DELTA_J2000'][0]
        grb_rad = mag_ecsvcleanSources[idx_GRBcleanpsf]['FLUX_RADIUS'][0]
        grb_snr = mag_ecsvcleanSources[idx_GRBcleanpsf]['SNR_WIN'][0]
        grb_dist = d2d[0]/u.deg

        print('Detected GRB ra = %.6f, dec = %.6f, with 50 percent flux radius = %.3f arcsec and SNR = %.3f' % (
        grb_ra, grb_dec, grb_rad, grb_snr))
        print('%s magnitude of GRB is %.2f +/- %.2f' % (filter, grb_mag, grb_magerr))

        grbdata = Table()
        grbdata['RA'] = np.array([grb_ra])
        grbdata['DEC'] = np.array([grb_dec])
        grbdata['%sMag' % filter] = np.array([grb_mag])
        grbdata['%sMag_Err' % filter] = np.array([grb_magerr])
        grbdata['Radius'] = np.array([grb_rad])
        grbdata['SNR'] = np.array([grb_snr])
        grbdata['Distance'] = np.array([grb_dist])

        grbdata.write('Source_Data.ecsv', overwrite=True)
        print('Generated GRB data table!')
    if len(idx_GRBcleanpsf) > 1:
        print('Multiple sources detected in search radius, refer to .ecsv file for source info!')
        mag_ar = []
        mag_err_ar =[]
        ra_ar =[]
        dec_ar = []
        rad_ar = []
        snr_ar = []
        dist_ar = []
        idx_GRBcleanpsflist = idx_GRBcleanpsf.tolist()
        for i in idx_GRBcleanpsflist:
            grb_mag = mag_ecsvcleanSources[idx_GRBcleanpsf]['%sMAG_PSF' % filter][idx_GRBcleanpsflist.index(i)]
            mag_ar.append(grb_mag)
            grb_magerr = mag_ecsvcleanSources[idx_GRBcleanpsf]['e_%sMAG_PSF' % filter][idx_GRBcleanpsflist.index(i)]
            mag_err_ar.append(grb_magerr)
            grb_ra = mag_ecsvcleanSources[idx_GRBcleanpsf]['ALPHA_J2000'][idx_GRBcleanpsflist.index(i)]
            ra_ar.append(grb_ra)
            grb_dec = mag_ecsvcleanSources[idx_GRBcleanpsf]['DELTA_J2000'][idx_GRBcleanpsflist.index(i)]
            dec_ar.append(grb_dec)
            grb_rad = mag_ecsvcleanSources[idx_GRBcleanpsf]['FLUX_RADIUS'][idx_GRBcleanpsflist.index(i)]
            rad_ar.append(grb_rad)
            grb_snr = mag_ecsvcleanSources[idx_GRBcleanpsf]['SNR_WIN'][idx_GRBcleanpsflist.index(i)]
            snr_ar.append(grb_snr)
            dist = (d2d[idx_GRBcleanpsflist.index(i)]).to(u.arcsec)
            dist = dist/u.arcsec
            dist_ar.append(dist)
        grbdata = Table()
        grbdata['RA (deg)'] = np.array(ra_ar)
        grbdata['DEC (deg)'] = np.array(dec_ar)
        grbdata['%sMag' % filter] = np.array(mag_ar)
        grbdata['%sMag_Err' % filter] = np.array(mag_err_ar)
        grbdata['Radius (arcsec)'] = np.array(rad_ar)
        grbdata['SNR'] = np.array(snr_ar)
        grbdata['Distance (arcsec)'] = np.array(dist_ar)
        grbdata.write('Source_Data.ecsv', overwrite=True)
        print('Generated GRB data table!')
    else:
        print('GRB source at inputted coords %s and %s not found, perhaps increase photoDistThresh?' % (ra,dec))

ra = 122.9204437
dec = -54.6520638
thresh = 0.7
imageName = 'coadd.Open-J.01015179-01015531.C2.fits'
survey = 'VHS'
filter = 'J'

GRB(ra,dec,imageName,survey,filter,thresh)


#%% Crowded field wcs value plots
from astropy.io import ascii

directory = '/mnt/d/PRIME_photometry_test_files/xrf_1/J_Band/C2_sub/'
cd1vals = []
for i in sorted(os.listdir(directory)):
    if i.endswith('.head'):
        opener = open(directory+i,'r')
        file = []
        for f in opener:
            file.append(f)
        phrase = 'CD2_1   ='
        substr = [i for i in file if phrase in i]
        substr = str(substr)
        cd1valsci = substr[13:32]
        cd1val = float(cd1valsci)
        cd1vals.append(cd1val)

cd1vals1 = np.where(np.array(cd1vals) < 1.8e-6)
cd1vals2 = np.where(np.array(cd1vals) > 1.5e-6)
cdvalsnew = np.intersect1d(cd1vals1,cd1vals2)
#%% Pruning out wcs outliers
imglist = []
for i in sorted(os.listdir(directory)):
    if i.endswith('.head'):
        imglist.append(i)

prunedlist = [imglist[i] for i in cdvalsnew]
prunedlistnum = [i[0:-5] for i in prunedlist]
#%% list of matching files and moving
import shutil
fullfilelist = []
for i in sorted(os.listdir(directory)):
    if any(string in i for string in prunedlistnum):
        shutil.copy(directory+i, '/mnt/d/PRIME_photometry_test_files/xrf_1/J_Band/C2_subprune/')

#%% plotting wcs
plt.figure(figsize=(10, 8))
plt.plot(cd1vals, 'ro')
plt.xlabel('Image Number')
plt.ylabel('CD1_2 Value')
#plt.ylim(-0.0001377,-0.0001385)
#plt.xticks(ticks)
plt.title('SCAMP .Head Values - CD1_2')
plt.show()


#%%
img = fits.open('/mnt/d/PRIME_photometry_test_files/BHNS/field9118/osaka/orig/01168799C1.geom_copy.fits')
data = img[0].data
hdr = img[0].header
dec = hdr['DEC-D']
ra = hdr['RA-D']
del hdr[5:-2]
hdr.insert(6, ('RA-D', ra, 'Right ascension'))
hdr.insert(7, ('DEC-D', dec, 'Declination'))
#hdr.insert(8,'END')

fits.writeto('/mnt/d/PRIME_photometry_test_files/BHNS/field9118/osaka/01168799C1.geom_nohdr.fits',data,hdr,overwrite=True)

#%%
opener = open('/mnt/d/PRIME_photometry_test_files/xrf_1/J_Band/flatfielded/C2_flatsub/01015179C2.sky.flat.head', 'r')
file = []
for f in opener:
    file.append(f)
phrasePIX1 = 'CRPIX1  ='
substrPIX1 = [i for i in file if phrasePIX1 in i]
substrPIX1 = str(substrPIX1)
crpix1valsci = substrPIX1[13:33]
crpix1val = float(crpix1valsci)

phrasePIX2 = 'CRPIX2  ='
substrPIX2 = [i for i in file if phrasePIX2 in i]
substrPIX2 = str(substrPIX2)
crpix2valsci = substrPIX2[13:33]
crpix2val = float(crpix2valsci)

phraseVAL1 = 'CRVAL1  ='
substrVAL1 = [i for i in file if phraseVAL1 in i]
substrVAL1 = str(substrVAL1)
crval1valsci = substrVAL1[13:33]
crval1val = float(crval1valsci)

phraseVAL2 = 'CRVAL2  ='
substrVAL2 = [i for i in file if phraseVAL2 in i]
substrVAL2 = str(substrVAL2)
crval2valsci = substrVAL2[13:33]
crval2val = float(crval2valsci)

img = fits.open('/mnt/d/PRIME_photometry_test_files/xrf_1/J_Band/flatfielded/C2_flatsub/01015179C2.sky.flat.fits')
data = img[0].data
hdr = img[0].header
hdr.set('CRPIX1',value = crpix1val)
hdr.set('CRPIX2',value = crpix2val)
hdr.set('CRVAL1',value = crval1val)
hdr.set('CRVAL2',value = crval2val)

fits.writeto('/mnt/d/PRIME_photometry_test_files/xrf_1/J_Band/flatfielded//01015179C2.sky.flat.new',data,hdr,overwrite=True)

#%% sextractor med combine

import os
from astropy.io import fits
direct = '/mnt/d/PRIME_photometry_test_files/sky_testing/sex_backs/'
stack = []
for f in os.listdir(direct):
    img = fits.getdata(direct+f)
    img = img / np.nanmedian(img)
    stack.append(img)
stackarray = np.stack(stack)
sex_med = np.nanmedian(stackarray,axis=0)

fits.writeto('/mnt/d/PRIME_photometry_test_files/sky_testing/sex_med.fits',sex_med,header=None,overwrite=True)

#%% sextractor subtraction
direct = '/mnt/d/PRIME_photometry_test_files/BHNS/field7994_2/C1_FF/'
out = '/mnt/d/PRIME_photometry_test_files/sky_testing/sex_subbed/'
sky = '/mnt/d/PRIME_photometry_test_files/sky_testing/sex_med.fits'
sky_img = fits.getdata(sky)
for f in os.listdir(direct):
    hdu = fits.open(direct+f)
    hdr = hdu[0].header
    img = hdu[0].data
    med_img = np.nanmedian(img)
    scaled = sky_img*med_img
    red_img = img - scaled
    fits.writeto(out+f,red_img,header=hdr,overwrite=True)

#%%
path = '/mnt/d/PRIME_photometry_test_files/BHNS/field9117_2/C4_sub/'
img_list = [f for f in sorted(os.listdir(path)) if f.endswith('.cat')]
img_list = [path+f for f in img_list]
img_list = ' '.join(img_list)

#%%
img = fits.open('/mnt/d/PRIME_photometry_test_files/SWIFT_1246989/C1/01504054C1.ramp.fits')
data = img[0].data
hdr = img[0].header
ra = fits.getval('/mnt/d/PRIME_photometry_test_files/SWIFT_1246989/C1/01504054C1.ramp.fits', 'RA-D')
dec = fits.getval('/mnt/d/PRIME_photometry_test_files/SWIFT_1246989/C1/01504054C1.ramp.fits', 'DEC-D')
fits.writeto('/mnt/d/PRIME_photometry_test_files/SWIFT_1246989/test.fits',data,overwrite=True)

#%%
from astroquery.vizier import Vizier
filter = 'Y'
raImage = 344
decImage = 1
boxsize = 32
catNum = 'II/319/las9'  # changing to 2mass
print('\nQuerying Vizier %s around RA %.4f, Dec %.4f, w/ box size %.2f' % (
    catNum, raImage, decImage, boxsize))
# You can set the filters for the individual columns (magnitude range, number of detections) inside the Vizier query
v = Vizier(columns=['RAJ2000', 'DEJ2000', '%smag' % filter, 'e_%smag' % filter],
           column_filters={"%smag" % filter: ">12"}, row_limit=-1)
print(v)
Q = v.query_region(SkyCoord(ra=raImage, dec=decImage, unit=(u.deg, u.deg)), width= str(boxsize)+ 'm', catalog=catNum, cache=False)
# query vizier around (ra, dec) with a radius of boxsize
print('Queried source total = ', len(Q[0]))
#%%
coordlist_test = ['344.57200,1.02675', '344.5756032,1.0135480']
coordlist_test = tuple(eval(i) for i in coordlist_test)
#%%
coordlist = (245,12),(123,-19),(192,-23)
listlist = {}
keys = np.arange(0,len(coordlist),1)
for i in range(len(keys)):
    print('examine RA =', coordlist[i][0])
    print('examine DEC =', coordlist[i][1])
    listlist[keys[i]] = coordlist[i], 5

for key in listlist:
    values = listlist[key]
    print('ra = ', values[0][0])
    print('location %d' % key)

#%%
fits_file = '/mnt/d/PRIME_photometry_test_files/img_sub_test_imgs/01556556C1.ramp.fits'
f = fits.open(fits_file, 'update')
header = f[0].header
img = f[0].data
del(header[97:100])

fits.writeto('/mnt/d/PRIME_photometry_test_files/img_sub_test_imgs/01556556C1.ramp.astrm.fits',data=img,header=header,overwrite=True)

#%%
fits_file = '/mnt/d/PRIME_photometry_test_files/img_sub_test_imgs/sci.fits'
f = fits.open(fits_file, 'update')
header = f[0].header
data = f[0].data
data = data[4:4092, 4:4092]

fits.writeto('/mnt/d/PRIME_photometry_test_files/img_sub_test_imgs/sci_crop.fits',data=data,header=header,overwrite=True)
#%% vega to ab conversion
mag = 19.5
filter = 'J'

zp_dict = {'zp_J': 1594, 'zp_H': 1024}

pick_zp = 'zp_'+filter
for k,v in zp_dict.items():
    if k == pick_zp:
        zp = v

flx = zp*10**(-mag/2.5)
ab_mag = -2.5*np.log10(flx/3631)
ab_mag = round(ab_mag,3)
print(ab_mag)

#%% new mflat filename gen
band = 'J'
chip = 1

example = ['mflat.J.20241017.C1.fits', 'mflat.J.20241010.C1.fits','mflat.J.20240424.C1.fits']

mflat_list = [f for f in example
              if f.endswith('.fits') if '.%s.' % band in f if 'C%s' % chip in f]

# get latest mflat
mflat_list = sorted(mflat_list, reverse=True)
filename = mflat_list[0]

#%%
import os
chip = 3
skypath = '/mnt/d/PRIME_photometry_test_files/GRB241105a_H/sky/'
filelist = [f for f in os.listdir(skypath) if f.endswith('.C{}.fits'.format(chip))]
if filelist:
    print(filelist[0])
    print('Previous sky found! Skipping sky gen..\n')
else:
    print('no!')

#%%  Check images for bck-bck and bck_sub-bck_sub
dir1 = '/mnt/photometry/supermaster_test/GB63_dither_20250406/J/C1_sub/'
dir2 = '/mnt/photometry/supermaster_test/GB61_dither_20250508/J/C1_sub/'
backdir = '/mnt/photometry/supermaster_test/GB_back_test/'

backs1 = [bck for bck in sorted(os.listdir(dir1)) if bck.endswith('.sky.flat.fits')]
backs2 = [bck for bck in sorted(os.listdir(dir2)) if bck.endswith('.sky.flat.fits')]
backs2 = backs2[:6]

for b1, b2 in zip(backs1, backs2):
    bimg1 = fits.getdata(os.path.join(dir1, b1))
    bimg2 = fits.getdata(os.path.join(dir2, b2))
    bsub = bimg1 - bimg2
    bname = b1[0:8] + '-' + b2[0:8] + '.sky.flat.fits'
    fits.writeto(os.path.join(backdir, bname), bsub, overwrite=True)

#%%
backsdiff = [bck for bck in sorted(os.listdir(backdir)) if bck.endswith('.back.fits')]
backsubdiff = [bck for bck in sorted(os.listdir(backdir)) if bck.endswith('.sky.flat.fits')]

for bd, bsd in zip(backsdiff, backsubdiff):
    bdiff = fits.getdata(os.path.join(backdir, bd))
    bsdiff = fits.getdata(os.path.join(backdir, bsd))

    print('File: ', bd)
    print('   Median Value over Image: %s' % np.median(bdiff))

    print('File: ', bsd)
    print('   Median Value over Image: %s' % np.median(bsdiff))

#%%
fzfile = fits.getheader('/mnt/photometry/S250328ae/field19010_20250424_REDO/Y/C1_astrom/02334960C1.ramp.new')
print(len(fzfile))
# hdr = fzfile[0].header
# if hdr['CHECKSUM']:
#     print('yes!')

