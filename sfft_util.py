import os
import numpy as np
import os.path as pa
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from mpl_toolkits.axes_grid1 import ImageGrid
from photutils.background import Background2D, MedianBackground
from scipy.ndimage import gaussian_filter


from astropy.io import fits, ascii
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.visualization import astropy_mpl_style, ZScaleInterval
from astropy.nddata import Cutout2D
import astropy.units as u
from astropy.wcs import WCS
from regions import Regions


from sfft.EasySparsePacket import Easy_SparsePacket
from sfft.EasyCrowdedPacket import Easy_CrowdedPacket
from sfft.CustomizedPacket import Customized_Packet
from sfft.utils.pyAstroMatic.PYSWarp import PY_SWarp
import time
from pathlib import Path
import glob

import photometrus.photometry.photometry as photometry
import photometrus_utils as util


def get_sci_ref(sci_stack, ref_stack):
	
	f_ref_p = Path(ref_stack)	   # reference 
	for j in f_ref_p.iterdir():
		# print(j)
		if 'coadd' in str(j) and '.fits' in str(j) and not '.fits.' in str(j) and not '.aligned.' in str(j) and not ".fits_" in str(j): 
			FITS_REF = str(j)
			print('Ref. Image: '+FITS_REF)
			
	
	f_sci_p = Path(sci_stack)	   # science 
	# print(f_sci_p)
	for j in f_sci_p.iterdir():
		if 'coadd' in str(j) and '.fits' in str(j) and not '.fits.' in str(j) and not '.aligned.' in str(j) and not ".fits_" in str(j):
			FITS_SCI = str(j)
			print('Sci. Image: '+FITS_SCI)


	return FITS_SCI, FITS_REF

def get_ra_dec(FITS_SCI):
	reg_path = f"{pa.dirname(FITS_SCI)}/GRB_query_thresh.reg"
	coord = Regions.read(reg_path)[0].center
	ra, dec = coord.ra.deg, coord.dec.deg
	return ra, dec


def get_force_conv(sci_stack,ref_stack):

	def get_med_flux_rad(path):
		ecsvs = glob.glob(f"{path}/*.ecsv")
		coadd_ecsvs = [f for f in ecsvs if os.path.basename(f).startswith("coadd")]
			
		ecsv_file = coadd_ecsvs[0]
		# print(ecsv_file)
		ecsv = pd.read_csv(ecsv_file, sep='\s+', comment='#')
		# print(ecsv.keys())
		flux_rad = np.median(ecsv['FLUX_RADIUS'])
		return flux_rad
		
	sci_flux_rad = get_med_flux_rad(sci_stack)
	ref_flux_rad = get_med_flux_rad(ref_stack)
	
	# def blur_image(path):
	#	 data = fits.getdata(path)
	#	 blurred = gaussian_filter(data, sigma=0.2)
		
	#	 new_header = fits.getheader(path).copy()
	#	 new_header.update(cutout.wcs.to_header())
	#	 fits.writeto(f"blur_{new_path}", blurred, new_header, overwrite=True)
	#	 return new_path
	
	# swap_ref_sci=False
	# print("Conv(Sci(flux rad ", sci_flux_rad, ")) - Ref(flux rad ", ref_flux_rad, ")")
	
	if ref_flux_rad < sci_flux_rad:
		ForceConv = 'REF'		  # FIXME {'AUTO', 'REF', 'SCI'}
		# FITS_REF = blur_image(FITS_REF)
		
	else:
		ForceConv = 'SCI'
		# FITS_SCI = blur_image(FITS_SCI)
	return ForceConv
	

# def get_bkg(data):
# 	bkg_estimator = MedianBackground()
# 	bkg = Background2D(data, box_size=64, filter_size=3,
# 						   bkg_estimator=bkg_estimator)
# 	return bkg
	

def get_bkg(image):
	data = fits.getdata(image)
	dir_data = os.path.dirname(image)
	
	bkg_estimator = MedianBackground()
	bkg = Background2D(data, box_size=64, filter_size=3,
					   bkg_estimator=bkg_estimator)

	
	plt.imshow(bkg.background, cmap='gray')
	bkg_path = f"{dir_data}/background.png"
	plt.savefig(bkg_path)
	return bkg, bkg_path

def remove_bkg_crop(img, filename, header, ra, dec, crop = 1000):
	data = fits.getdata(img)
	dir_data = os.path.dirname(img)

	bkg, bkg_path = get_bkg(img)
	data_sub = data #- bkg.background


	new_path = f"{dir_data}/{filename}"

	wcs = WCS(header)
	xpx, ypx = wcs.world_to_pixel(SkyCoord(ra=ra*u.deg, dec=dec*u.deg))
	pxthresh = 50
	xpxmax, ypxmax = data_sub.shape

	print("transient at", xpx, ypx)
	print("image dims", xpxmax, ypxmax)



	def crop_img(data, crop):
		# crop image to 1000, or largest image that i can 
		if crop == 0:
			fits.writeto(new_path, data, header, overwrite=True)
			return data
		
			# 3500 < 50 + 400	   (4500 - 3500 ) < 50 + 400		 700 < 50 + 400		   4500 - 700 < 50 + 400
		if xpx < pxthresh + crop or xpxmax - xpx < pxthresh + crop or ypx < pxthresh + crop or ypxmax - ypx < pxthresh + crop:
			return crop_img(data, crop-100)
			
		else:
			size = (xpxmax-crop, ypxmax-crop)
			center = (xpxmax//2, ypxmax//2)
			print(size, center)
			cutout = Cutout2D(data, center, size, wcs=wcs, copy=True)

			new_header = header.copy()
			new_header.update(cutout.wcs.to_header())

			fits.writeto(new_path, cutout.data, new_header, overwrite=True)
			
			print(xpxmax,ypxmax, "->", cutout.data.shape)	
			print("wrote to", new_path)


			return cutout.data
		
	data_crop = crop_img(data_sub, 0) # dont crop for now

	return new_path, bkg_path


def get_all_coords(sci_stack, ref_stack):
	
	cat_files = glob.glob(f"{sci_stack}/*.cat")
	cat_files.extend(glob.glob( f"{ref_stack}/*.cat"))
	
	all_coords = []
	for f in cat_files:
		cat = Table.read(f, hdu=2)
		px_per_arcsec = 1 / 0.49766 # 49766 arcsec / px
	
		# FIX which flux rad col? first seems to be smaller so ill go with that
	
		flux_rad_col = cat['FLUX_RADIUS'] if len(cat['FLUX_RADIUS'].shape)==1 else cat['FLUX_RADIUS'][:,0]
		mask = (cat['FLAGS'] == 0)  & (abs(flux_rad_col) > px_per_arcsec) & (cat['ELONGATION'] < 2) # no issues and flux radius > 1 arcsec
		# crop_mask = 
		
		x_good = cat['X_IMAGE'][mask]
		y_good = cat['Y_IMAGE'][mask]
		coords = np.column_stack((x_good, y_good)).tolist()
		all_coords.extend(coords) 
	
	print(np.shape(all_coords))
	all_coords = np.array(all_coords)
	return all_coords


def rms(sky, name):
	sky = sky[np.isfinite(sky)] # remove nans
	median = np.median(sky)
	mad = np.median(np.abs(sky - median)) # mean absolute deviation
	rms = 1.4826 * mad
	print(name, "RMS:",rms)
	return rms

def rms_around_sources():
	# calculate stdev around sources only, mask space between
	for f in cat_files:
		cat = Table.read(f, hdu=2)
		for source in cat.rows():
			ra, dec, flux_rad = source['RA'], source['dec']
			# mask = mask & (ra>
	

def get_diff_stats(SFFTPrepDict, PixA_DIFF):
	
	print(SFFTPrepDict.keys())
	
	# mag offset: flux scaling
	for key in ['SATLEVEL_REF', 'SATLEVEL_SCI', 'MAG_OFFSET', 'FWHM_REF', 'FWHM_SCI']:
		print(key, ":",SFFTPrepDict[key])
	
	# use masked images?
	print()
	rms_ref = rms(SFFTPrepDict['PixA_REF'], "REF")
	rms_sci = rms(SFFTPrepDict['PixA_SCI'], "SCI")
	rms_diff_predicted = (rms_ref**2 + rms_sci**2)**(1/2)
	print("dif should be:", rms_diff_predicted)
	
	mask_diff = (~SFFTPrepDict['Active-Mask'] & ~SFFTPrepDict['Union-NaN-Mask'])
	rms_diff = rms(PixA_DIFF[mask_diff], "DIFF")

	return rms_ref, rms_sci, rms_diff_predicted, rms_diff



def ztf_cuts(cat, PixA_DIFF):

	def count_neg(source):
		x,y = source['X_IMAGE'], source['Y_IMAGE']
		x, y = round(x), round(y)
		# print(x, y)
		# stamp = PixA_DIFF[y-10:y+10, x-10:x+10]
		
		stamp = PixA_DIFF[y-2:y+3, x-2:x+3]
		n_neg = np.sum(stamp < 0)
		
		# hot pixels, dead pixels, ...
		#n_bad = np.sum(badmask[y-2:y+3, x-2:x+3])
		
		return n_neg
	

	# these are for the difference, but they should certainly be true of the science image then
	
	# remove multi dimensional cols
	
	names = [name for name in cat.colnames if len(cat[name].shape) <= 1]
	cat = cat[names]
	df = cat.to_pandas()
	
	# # missing:
	#	   # SNR for 8 px diameter aperture >5
	#	   # 0 < Flux (8 px) / Flux (18 px) < 1.5
	#	   # Number of bad pixels in a 5 × 5 pixel area <= 7; 

	# 300 px from edge

	# uncomment for cuts
	crop = 300
	x_px_max, y_px_max = np.shape(PixA_DIFF)
	
	crop_mask = (df['X_IMAGE'] > crop) & (x_px_max - df['X_IMAGE'] > crop) & (df['Y_IMAGE'] > crop) & (y_px_max - df['Y_IMAGE'] > crop)
	df = df[crop_mask]

	df = df[df["SNR_WIN"]>5] 
   
	df = df[df["ELONGATION"] <= 2]
	
	df["N_NEG"] = df.apply(count_neg,axis=1) #wrong axis?
	df = df[df["N_NEG"] <= 13] #13]  
	
	return df
	
	
def source_extract(img):
	grb_radius = "4.0"
	crop = 150 
	comp_lvl = 0.3
	directory = os.path.dirname(img)
	name = os.path.basename(img)
	grb_thresh = photometry.grb_rad_convert(grb_radius)
	
	data, header, w, raImage, decImage, bulge, det_thresh, chip = photometry.img(directory, name, crop) # basic image info
	# Q, chosen_survey, mag_low_cutoff = photometry.query(raImage, decImage, band, w, data, crop) # queries catalog (vizier)
 
	catalogName = util.sex1(name, det_cut=det_thresh)
	print("wrote to ",  catalogName)
	cat = Table.read(catalogName, hdu=2)

	# psf_file = glob.glob(FITS_SCI_OLD+'.psf')[0]
	# # if psf_files == []:
	# #	 psf_files = glob.glob(FITS_REF_OLD+'.psf')
	# # psf_file = psf_files[0]
	# print("PSF file:", psf_file)

	# psfcatalogName = util.sex2(name, det_cut=det_thresh, catalogName=catalogName, psfName = psf_file) # raw source extractor catalog
	# psf_cat = Table.read(psfcatalogName, hdu=2) # use this! 

	# good_cat_stars, cleanPSFSources, PSFSources, idx_psfmass, idx_psfimage, massCatCoords = photometry.tables(Q, data, w, psfcatalogName,
	#																				 crop, None) # this runs a crossmatch, returns good psf sources (not in middle, point sources, ...)
	# idx_psfmass, idx_psfimage, are idx with crossmatches in vizier
	# PSFSources

	# magnitudes in sextractor arent correct without photometry, run grb photometry on source 

	return cat

def closest_source(good_ztf_sources, FITS_DIFF, ra, dec):
	threshold = 0.005
	good_ztf_sources["distance"] = ((good_ztf_sources['ALPHA_J2000']-ra)**2 + (good_ztf_sources['DELTA_J2000']-dec)**2)**(1/2)
	min_source_idx = good_ztf_sources["distance"].argmin()
	
	source = good_ztf_sources.iloc[min_source_idx]
	source_ra = source['ALPHA_J2000']
	source_dec = source['DELTA_J2000']
	savename, threshname = util.make_cutout(FITS_DIFF, source_ra, source_dec, display_file=False, name_ext="extracted")
	
	# for i in range(len(good_ztf_sources)):
	   
	# 	source = good_ztf_sources.iloc[i]
	# 	source_ra = source['ALPHA_J2000']
	# 	source_dec = source['DELTA_J2000']
	# 	if abs(source_ra-ra) < threshold and abs(source_dec-dec) < threshold:
	# 		print("Found source!")
	# 		print(source_ra, source_dec)
	# 		print(source)
		  

	return savename, threshname, source_ra, source_dec


# def SFFT(dir_data, f_sci, f_ref, band):
def SFFT(FITS_SCI, FITS_REF):

	dir_data = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(FITS_SCI))))
	print(dir_data)
	
	t = time.time()
	FILE_BASENAME = f"sub_3_26_{t}"
	
	# * computing backend and resourse 
	BACKEND_4SUBTRACT = 'Numpy'	 # FIXME {'Cupy', 'Numpy'}, Use 'Numpy' if you only have CPUs
	CUDA_DEVICE_4SUBTRACT = '0'	 # FIXME ONLY work for backend Cupy
	NUM_CPU_THREADS_4SUBTRACT = 10   # FIXME ONLY work for backend Numpy
	
	# * required info in FITS header
	GAIN_KEY = 'GAIN'  # not GAIN_CAL			 # NOTE Keyword of Gain in FITS header
	SATUR_KEY = 'SATURATE'		  # NOTE Keyword of Saturation in FITS header
	
	# * how to subtract
	GKerHW = None				   # FIXME given matching kernel half width
	KerHWRatio = 2.0				# FIXME Ratio of kernel half width to FWHM (typically, 1.5-2.5).
	KerPolyOrder = 2			  # FIXME {0, 1, 2, 3}, Polynomial degree of kernel spatial variation
	BGPolyOrder = 0 # trivial for sparse (already sky subtracted)			# As above but for CROWDED field
	ConstPhotRatio =  False	  #False	# FIXME Constant photometric ratio between images? dont scale them
	PriorBanMask = None			 # FIXME None or a boolean array with same shape of science/reference.
	
	COARSE_VAR_REJECTION = False #True	 # FIXME Coarse Variable Rejection? {True, False}
	CVREJ_MAGD_THRESH = 0.12		# FIXME magnitude threshold for Coarse Variable Rejection
	ELABO_VAR_REJECTION = False #True	  # FIXME Elaborate Variable Rejection? {True, False}

	sci_stack = os.path.dirname(FITS_SCI)
	ref_stack = os.path.dirname(FITS_REF)
	
	# FITS_SCI, FITS_REF, ra, dec = get_sci_ref(sci_stack, ref_stack)
	ra, dec = get_ra_dec(FITS_SCI)
	ForceConv =  get_force_conv(sci_stack,ref_stack)

	sci_header = fits.getheader(FITS_SCI)
	ref_header = fits.getheader(FITS_REF)
	
	FITS_REF_OLD, FITS_SCI_OLD = FITS_REF, FITS_SCI
	FITS_SCI, sci_bkg_path = remove_bkg_crop(FITS_SCI, f"{FILE_BASENAME}_science_bkgsub.fits", sci_header, ra, dec)
	FITS_REF, ref_bkg_path = remove_bkg_crop(FITS_REF, f"{FILE_BASENAME}_reference_bkgsub.fits", ref_header, ra, dec)
	
	# align with Orion's method
	FITS_REF_OLD = FITS_REF
	FITS_SCI, FITS_REF = util.multi_epoch_astrom(FITS_SCI, FITS_REF)

	
	FITS_DIFF = dir_data+'/'+FILE_BASENAME+'.sfftdiff.fits'			# difference
	FITS_REF_al = FITS_REF[:-5] + '.aligned.fits'   # refernce aligned
	
	
	# align with swarp
	PY_SWarp.PS(FITS_obj=FITS_REF, FITS_ref=FITS_SCI, FITS_resamp=FITS_REF_al, \
		GAIN_KEY=GAIN_KEY, SATUR_KEY=SATUR_KEY, OVERSAMPLING=1, RESAMPLING_TYPE='LANCZOS3', \
		SUBTRACT_BACK='N', FILL_VALUE=np.nan, VERBOSE_TYPE='NORMAL', VERBOSE_LEVEL=2)
	
	print('\nMeLOn CheckPoint: IMAGE ALIGNMENT WITH SWARP DONE!\n')
	
	
	print('Ref. Image aligned: '+FITS_REF_al)
	print('Diff. Image:'+FITS_DIFF)

	all_coords = get_all_coords(sci_stack, ref_stack)
	

	PixA_DIFF, SFFTPrepDict = Easy_SparsePacket.ESP(FITS_REF=FITS_REF_al, FITS_SCI=FITS_SCI,
								FITS_DIFF=FITS_DIFF, FITS_Solution=None, ForceConv=ForceConv, GKerHW=GKerHW,
								KerHWRatio=KerHWRatio, KerHWLimit=(5, 20), KerPolyOrder=KerPolyOrder, 
								BGPolyOrder=BGPolyOrder, ConstPhotRatio=ConstPhotRatio, MaskSatContam=False, 
								GAIN_KEY=GAIN_KEY, SATUR_KEY=SATUR_KEY, BACK_TYPE='MANUAL', BACK_VALUE=0.0, 
								BACK_SIZE=64, BACK_FILTERSIZE=2, DETECT_THRESH=2, DETECT_MINAREA=5, 
								DETECT_MAXAREA=0, DEBLEND_MINCONT=1e-4, BACKPHOTO_TYPE='LOCAL', 
								ONLY_FLAGS=[0], BoundarySIZE=30, XY_PriorSelect=all_coords, PointSource_MINELLIP=0.3, MatchTol=None, 
								MatchTolFactor=3.0, StarExt_iter=4, XY_PriorBan=None,
								PostAnomalyCheck=False, PAC_RATIO_THRESH=5.0, BACKEND_4SUBTRACT=BACKEND_4SUBTRACT, 
								CUDA_DEVICE_4SUBTRACT=CUDA_DEVICE_4SUBTRACT,
								NUM_CPU_THREADS_4SUBTRACT=NUM_CPU_THREADS_4SUBTRACT)[:2]
	
	# (2, 20)
	
	# MAG_OFFSET = -2.5 * np.log10(ConstPhotRatio) # 0.5
	# SFFTPrepDict["MAG_OFFSET"] = MAG_OFFSET



	  
	# PixA_DIFF, SFFTPrepDict = Easy_SparsePacket.ESP(FITS_REF=FITS_REF_al, FITS_SCI=FITS_SCI, \
	# 								 FITS_DIFF=FITS_DIFF, FITS_Solution=None, ForceConv=ForceConv, GKerHW=None, \
	# 								 KerHWRatio=KerHWRatio, KerHWLimit=(2, 20), KerPolyOrder=KerPolyOrder, \
	# 								 BGPolyOrder=BGPolyOrder, ConstPhotRatio=ConstPhotRatio, MaskSatContam=False, \
	# 								 GAIN_KEY=GAIN_KEY, SATUR_KEY=SATUR_KEY, BACK_TYPE='MANUAL', BACK_VALUE=0.0, \
	# 								 BACK_SIZE=64, BACK_FILTERSIZE=2, DETECT_THRESH=2, DETECT_MINAREA=5, \
	# 								 DETECT_MAXAREA=0, DEBLEND_MINCONT=1e-4, BACKPHOTO_TYPE='LOCAL', \
	# 								 ONLY_FLAGS=[0], BoundarySIZE=30, XY_PriorSelect=None, Hough_MINFR=0.1, \
	# 								 Hough_PeakClip=0.7, BeltHW=0.2, PointSource_MINELLIP=0.3, MatchTol=None, \
	# 								 MatchTolFactor=3.0, COARSE_VAR_REJECTION=COARSE_VAR_REJECTION, \
	# 								 CVREJ_MAGD_THRESH=CVREJ_MAGD_THRESH, ELABO_VAR_REJECTION=ELABO_VAR_REJECTION, \
	# 								 EVREJ_RATIO_THREH=5.0, EVREJ_SAFE_MAGDEV=0.04, StarExt_iter=4, XY_PriorBan=None, \
	# 								 PostAnomalyCheck=False, PAC_RATIO_THRESH=5.0, BACKEND_4SUBTRACT=BACKEND_4SUBTRACT, \
	# 								 CUDA_DEVICE_4SUBTRACT=CUDA_DEVICE_4SUBTRACT, \
	# 								 NUM_CPU_THREADS_4SUBTRACT=NUM_CPU_THREADS_4SUBTRACT)[:2]

	# util.display(FITS_SCI)
	# util.display(FITS_REF)
	# util.display(FITS_DIFF)
	


	# get_diff_stats(SFFTPrepDict, PixA_DIFF)

	# sources = source_extract(FITS_DIFF)
	
	
	# #TODO add a check here to see if source ra, dec in sources 
	# good_ztf_sources = ztf_cuts(sources)
	
	# print(np.shape(PixA_DIFF))
	# print(len(good_ztf_sources), "sources")
	# print(ra, dec)
	
	# threshold = 0.005
	# for i in range(len(good_ztf_sources)):
	   
	# 	source = good_ztf_sources.iloc[i]
	# 	source_ra = source['ALPHA_J2000']
	# 	source_dec = source['DELTA_J2000']
	# 	if abs(source_ra-ra) < threshold and abs(source_dec-dec) < threshold:
	# 		print("Found source!")
	# 		print(source_ra, source_dec)
	# 		print(source)
		  
	# 		savename, threshname = util.make_cutout(FITS_DIFF, source_ra, source_dec)

	# sci_bkg, ref_bkg, sci_cutout, ref_cutout
	return FITS_SCI, FITS_REF, FITS_DIFF, PixA_DIFF, SFFTPrepDict, sci_bkg_path, ref_bkg_path



	

	
	# sci_cat = sfft_util.source_extract(FITS_SCI)
	# # print(sci_cat)
	# ref_cat = sfft_util.source_extract(FITS_REF)
	
	# print("Sources in sci:",len(sci_cat))
	# print("Sources in ref:",len(ref_cat))
	# print("Sources in dif:",len(sources))











