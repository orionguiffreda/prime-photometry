#%%
import os
import numpy as np
from astropy.io import ascii
import matplotlib.pyplot as plt
from astropy.io import fits
import subprocess
from astropy.table import Table
#%%
directory = '/mnt/photometry/S250328ae/field18699_20250419/H/C3_sub/'

# distortparams = ['PV1_0', 'PV1_1', 'PV1_2', 'PV1_4', 'PV1_5', 'PV1_6', 'PV2_0', 'PV2_1', 'PV2_2', 'PV2_4', 'PV2_5', 'PV2_6']
#
# distortvals = {}

astrmrms = []

for hdr in sorted(os.listdir(directory)):
    if hdr.endswith('.head'):
        opener = open(directory+hdr,'r')
        file = []
        for f in opener:
            file.append(f)
        phrase = 'ASTRRMS1='
        substr = [i for i in file if phrase in i]
        substr = str(substr)
        cd1valsci = substr[13:32]
        cd1val = float(cd1valsci)
        nameval = (cd1val, hdr)
        astrmrms.append(nameval)

astrvals = [[f[0],v] for v,f in enumerate(astrmrms)]
med = np.median(astrvals)
stdev = np.std(astrvals)
print('med = %.6f, stdev = %.6f' % (med, stdev))
outliers = [f for f in astrvals if f >= med+stdev or f <= med-stdev]

print('# of outliers = ', len(outliers))
print('outliers = ', outliers)
#%%
plt.figure(figsize=(10,8))
plt.plot(outliers)
#plt.axhline(med)
plt.title('ASTRRMS1 Values = Problematic field (field 18699)')
plt.xlabel('Outliers')
plt.ylabel('Values')
# plt.legend('Median')
plt.show()
#%%

directory = '/mnt/photometry/supermaster_test/GB_sparseastrom/preserve/'

# distortparams = ['PV1_0', 'PV1_1', 'PV1_2', 'PV1_4', 'PV1_5', 'PV1_6', 'PV2_0', 'PV2_1', 'PV2_2', 'PV2_4', 'PV2_5', 'PV2_6']
#
# distortvals = {}

astrmrms = []

for hdr in sorted(os.listdir(directory)):
    if hdr.endswith('flat.head'):
        opener = open(directory+hdr,'r')
        file = []
        for f in opener:
            file.append(f)
        phrase = 'ASTRRMS1='
        substr = [i for i in file if phrase in i]
        substr = str(substr)
        cd1valsci = substr[13:32]
        cd1val = float(cd1valsci)
        astrmrms.append(cd1val)

med = np.median(astrmrms)
std = np.std(astrmrms)
print('ASTRRMS med = %.10f +- %.10f' % (med,std))
#%%
deg2med = 0.0001399214
deg2std = 0.0000040055

deg3med = 0.0000974455
deg3std = 0.0000025863

deg4med = 0.0000959714
deg4std = 0.0000026513

deg5med = 0.0000961605
deg5std = 0.0000027146

deg6med = 0.0000957148
deg6std = 0.0000026503

deg7med = 0.0000957394
deg7std = 0.0000027542

degmeds = [deg3med, deg4med, deg5med, deg6med, deg7med]
degstds = [deg3std, deg4std, deg5std, deg6std, deg7std]

plt.figure(figsize=(10,8))
degx = [3,4,5,6,7]
plt.errorbar(degx, degmeds, yerr=degstds, fmt='ro', capsize=10)
plt.ylim(0.000093, 0.000099)
# plt.xlim(2.5,7.5)
plt.xticks([3,4,5,6,7])
plt.grid()
plt.title('ASTRRMS1 Median Values by SCAMP Polynomial Degrees')
plt.xlabel('DISTORT_DEGREES Value')
plt.ylabel('ASTRRMS1 Median Value')
plt.savefig('/mnt/photometry/TransientEvents/GRB250309B_20250310/J/C1_sub/polydegplot3plus.png', dpi=300)


#%% Crossmatch prime and vvv to find offsets
from photometrus.settings import gen_config_file_name
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astroquery.vizier import Vizier
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable

directory = '/mnt/photometry/supermaster_test/GB_offsets/gb63_again/'
imageName = '02482975C1.sky.flat.fits'
band = 'J'
maglow = 12
maghigh = 15

def sex1(imageName):
    print('Running sextractor for psf...')
    configFile = gen_config_file_name('bulge_new.config')
    paramName = gen_config_file_name('tempsource.param')
    catname = imageName + '.cat'
    try:
        command = 'sex %s -c %s -CATALOG_NAME %s -PARAMETERS_NAME %s -CHECKIMAGE_TYPE NONE' % (imageName, configFile, catname, paramName)
        #print('Executing command: %s' % command)
        rval = subprocess.run(command.split(), check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError as err:
        print('Could not run sextractor with exit error %s'%err)
    return catname

def boxchange(size):
    # size = size of box in arcmin
    boxsize_pix = size * 60 * (1/0.498)
    crop_pix = (4088 - boxsize_pix) / 2
    return size, crop_pix

def imaging(directory, imageName):
    os.chdir(directory)
    f = fits.open(os.path.join(directory, imageName))
    data = f[0].data  # This is the image array
    header = f[0].header

    # strong the image WCS into an object
    w = WCS(header)

    # Get the RA and Dec of the center of the image
    [raImage, decImage] = w.all_pix2world(data.shape[0] / 2, data.shape[1] / 2, 1)

    # Get zero point for image
    chip = header['CHIP']
    if header['FILTER2'] == 'Open':
        band = 'Z'
    else:
        band = header['FILTER2']

    return data, header, w, raImage, decImage

def cat_query(raImage, decImage, band, boxsize=34, maglow=12.5, maghigh=14.5):
    columns = ['RAJ2000', 'DEJ2000', 'RAICRS', 'DEICRS', "%sPSF" % band.lower(), "%smag" % band, "%smag3" % band,
               'srcid']
    coords = SkyCoord(ra=raImage*u.degree, dec=decImage*u.degree, frame='fk5')
    galcoords = coords.galactic
    print('Using bulge field query, converting coords to galactic: l = %s, b = %s' % (galcoords.l.deg, galcoords.b.deg))
    catNum = 'II/348/vvv2'  # changing to VVV
    print('\nQuerying Vizier %s around l %.4f, b %.4f, w/ box size %.2f, mag lim of %s - %s' % (catNum,
        galcoords.l.deg, galcoords.b.deg, boxsize, maglow, maghigh))
    try:
        v = Vizier(columns=columns, column_filters={"%smag3" % band: "%s .. %s" % (maglow, maghigh),
                                                    "%sperrbits" % band: "<=16", "Nd": ">6"}, row_limit=-1)
        Q = v.query_region(SkyCoord(galcoords, unit=(u.deg, u.deg)), width=str(boxsize) + 'm',
                           catalog=catNum, cache=False, frame='galactic')
        print('Queried source total = ', len(Q[0]))
    except:
        print('Issue with bulge-specific query.')
    return Q

def make_tables(directory, catname, Q):
    print('Creating sorted catalog & prime tables...')
    sexcat = Table.read(os.path.join(directory, catname), hdu=2)
    colnames = Q[0].colnames
    RA = colnames[0]
    DEC = colnames[1]

    inner_catsources = Q[0]

    inner_primesources = sexcat[(sexcat['FLAGS'] == 0) & (sexcat['FLUX_RADIUS'] * 0.498 > 1)]
    print('# of prime sources =', len(inner_primesources))

    VVVCatCoords = SkyCoord(ra=inner_catsources[RA], dec=inner_catsources[DEC], frame='icrs', unit='degree')
    PRIMECatCoords = SkyCoord(ra=inner_primesources['ALPHA_J2000'], dec=inner_primesources['DELTA_J2000'], frame='icrs', unit='degree')

    photoDistThresh = 0.6

    idx_PRIME, idx_vvv, d2d, d3d = VVVCatCoords.search_around_sky(PRIMECatCoords,
                                                                          photoDistThresh * u.arcsec)

    print('Crossmatch # = ',len(idx_PRIME))

    return inner_primesources, inner_catsources, idx_PRIME, idx_vvv

def gen_offset_plot(inner_primesources, inner_catsources, idx_PRIME, idx_vvv, directory, imagename):
    # Calculate statistics for the image
    imgpath = os.path.join(directory, imagename)
    fitsimg = fits.getdata(imgpath)
    mean, median, sigma = sigma_clipped_stats(fitsimg)

    # crossmatched catalogs
    vvv_match = inner_catsources[idx_vvv]
    prime_match = inner_primesources[idx_PRIME]

    # Mag offset values
    prime_mags = prime_match['MAG_AUTO']
    vvv_mags = vvv_match['Jmag3']
    diff_mags = abs(prime_mags - vvv_mags)
    # min_diff = np.min(diff_mags)
    # max_diff = np.max(diff_mags)
    min_diff = 21
    max_diff = 25
    """
    # Normalize the ellipticity values
    norm_diff = plt.Normalize(min_diff, max_diff)
    cmap = matplotlib.cm.jet

    # ScalarMappable for the color bar
    sm = plt.cm.ScalarMappable(cmap='jet', norm=norm_diff)

    # Create the plot
    print('Generating plot of crossmatched mag offsets betw. prime and catalogue...')
    fig = plt.figure(figsize=(12, 12))
    ax = fig.gca()

    # Plot the image in grayscale with adjusted contrast
    im = plt.imshow(fitsimg, vmin=-300, vmax=300, cmap='bone', origin='lower')

    # Plot the circles
    for i in range(len(prime_match['X_IMAGE'])):
        color = cmap(norm_diff(diff_mags[i]))  # Normalize the actual ellipticity value
        circle = plt.Circle((prime_match['X_IMAGE'][i], prime_match['Y_IMAGE'][i]), radius=5, color=color, fill=True)
        ax.add_patch(circle)

    # Add the color bar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    plt.colorbar(sm, orientation='vertical', cax=cax)

    ax.set_title('Mag Offset (PRIME v. VVV) Heatmap Over FITS Image: GB61')
    ax.set_xlabel('X Coordinate')
    ax.set_ylabel('Y Coordinate')
    plt.savefig(os.path.join(directory, 'magoffsetmap.png'), dpi=200)
"""
    # histogram

    print('making histogram plot...')
    plt.figure(2, figsize=(10, 8))
    plt.grid()
    plt.hist(diff_mags, bins=100, density=True, edgecolor='black')
    plt.xlim(20,26)
    plt.xticks(np.arange(20,26.5,0.5))
    plt.xlabel('Normalized Difference in Magnitude')
    plt.ylabel('Number of Crossmatched Sources')
    plt.title('Mag Offset (PRIME v. VVV) Histogram: GB63')
    plt.savefig('magoffsethist.png', dpi=200)
    plt.close('all')
    print('done!')



os.chdir(directory)
catname = sex1(imageName)
data, header, w, raImage, decImage = imaging(directory, imageName)
Q = cat_query(raImage, decImage, band, boxsize=34, maglow=maglow, maghigh=maghigh)
inner_primesources, inner_catsources, idx_PRIME, idx_vvv = make_tables(directory, catname, Q)
gen_offset_plot(inner_primesources, inner_catsources, idx_PRIME, idx_vvv, directory, imageName)


#%% bulge check
from photometrus.getfiles import replace_list

gp_name_list = [item for pair in replace_list for item in pair]
plane_list = [name for name in gp_name_list if 'plane' in name]

bulge_list = ['bulge', 'gb', 'gp', 'plane']
# bulge_plane_list = bulge_list + plane_list

case = 'Galactic_Plane'
bulge_check = [name for name in bulge_list if name in case.lower()]
if bulge_check:
    print('true!')
else:
    print('false!')

#%% sextractor sky normalization

sxsky = fits.open('/mnt/photometry/sky_methods_testing/GB94_dither_sexsky/C4_sub/02284954C4.sky.flat.back.fits')
sxskydata = sxsky[0].data

sxskynormdata = sxskydata / np.nanmedian(sxskydata)

fits.writeto('/mnt/photometry/sky_methods_testing/skies/02284954C4.norm.back.fits',
             sxskynormdata, sxsky[0].header)
#%% checkplots
def checkplot(output_directory, save_name):
    print('Generating histogram check plot!\n')
    skypath = os.path.join(output_directory, save_name)
    print('Statistics for %s:' % skypath)
    skyimg = fits.getdata(skypath)
    stackmed = np.nanmedian(skyimg)
    stackstd = np.std(skyimg)
    print(' Med: %.6f +/- stdev: %.4f' % (stackmed, stackstd))

    flat_sky = skyimg.flatten()
    mask = (flat_sky >= -9000) & (flat_sky <= 22000)
    flat_sky = flat_sky[mask]

    plt.figure(figsize=(10, 8))
    plt.hist(flat_sky, bins=25, edgecolor='black')
    plt.xlim(-9000, 22000)
    plt.ylim(1, 10**8)
    plt.yscale('log')
    plt.xlabel('Pixel Value')
    plt.ylabel('Frequency')
    plt.title('%s Image Histogram' % save_name)

    plt.text(10000, 10**7, ' Med: %.3f +/- stdev: %.3f' % (stackmed, stackstd), fontsize=11,
             bbox=dict( alpha=0.5, boxstyle='round' ,pad=0.5))
    plt.savefig('%s.check_plot.png' % skypath, dpi=300)

# checkplot('/mnt/photometry/sky_methods_testing/skies/', 'sky.Open-J.02283146-02283182.C1.fits')
# checkplot('/mnt/photometry/sky_methods_testing/skies/', 'sky.Open-J.02274417-02274437.C1.fits')
# checkplot('/mnt/photometry/sky_methods_testing/skies/', '02284954C1.norm.back.fits')

checkplot('/mnt/photometry/sky_methods_testing/stacks', 'coadd.GCsky.C4.fits')
checkplot('/mnt/photometry/sky_methods_testing/stacks', 'coadd.othersky.C4.fits')
checkplot('/mnt/photometry/sky_methods_testing/stacks', 'coadd.sxsky.C4.fits')

#%%
import astropy.units as u

# Example: AB magnitude array
ABmag = np.array([20.0, 21.5, 23.0]) * u.ABmag  # replace with your data

# Convert AB magnitude to flux density in Jy
# flux_jy = (3631 * 10**(-0.4 * u_ABmag)) * u.Jy

# Convert to mJy
flux_mjy = ABmag.to(u.microjansky)

print(flux_mjy)

#%%
tab = Table.read('/mnt/photometry/sky_methods_testing/GB110/GB110_dither_allsky/H/C4_astrom/02714994C4.sky.flat.psf.cat', hdu=2)
mask = [str(y).startswith('1875.') for y in tab['Y_IMAGE']]
filtered = tab[mask]
xims = filtered['X_IMAGE']
xs = [f for f in xims]
#%%
chosen = filtered['MAG_AUTO'][21]
cal_mag = 24.111 + chosen
print('Cal mag = ',cal_mag)
