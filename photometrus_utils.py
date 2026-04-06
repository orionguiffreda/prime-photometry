import photometrus.photometry.photometry as photometry
from photometrus.stack.stack import swarp_sx, swarp_missfits
from photometrus.astrom.astrom_img_sub import multi_epoch_scamp

import subprocess
import os
import numpy as np
import matplotlib.pyplot as plt
import glob
import shutil


import astropy
from astropy.io import fits
from astropy.nddata import Cutout2D

from regions import CircleSkyRegion
from astropy.stats import sigma_clip, sigma_clipped_stats
import astropy.units as u
from astropy.coordinates import Angle, SkyCoord
from astropy.wcs import WCS
from astropy.visualization import astropy_mpl_style, ZScaleInterval


import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from photutils.aperture import aperture_photometry, SkyCircularAperture
import re



def psfex(catalogName):
	print('Running PSFEx…')

	psfConfigFile = photometry.gen_config_file_name('default.psfex')
	base = os.path.splitext(os.path.basename(catalogName))[0]
	psfImageName = f'PSF_{base}'

	command = [
		'psfex',
		catalogName,
		'-c', psfConfigFile,
		'-CHECKIMAGE_TYPE', 'SNAPSHOTS',
		'-CHECKIMAGE_NAME', psfImageName
	]

	subprocess.run(command, check=True)

	
	psf_files = glob.glob(f'PSF_{base}*.fits')
	if not psf_files:
		raise FileNotFoundError("No PSFEx output found")

	print(psf_files)
	psfImageName = psf_files[0]

	return psfImageName


def sex1(imageName, det_cut, magtype='PSF'):
	print('Running sextractor on img to initially find sources...')
	if magtype == 'PSF':
		configFile = photometry.gen_config_file_name('sex2.config')
		paramName = photometry.gen_config_file_name('tempsource.param')
		catalogName = imageName + '.cat'
	else:
		configFile = photometry.gen_config_file_name('sex2.config')
		paramName = photometry.gen_config_file_name('photomAUTO.param')
		catalogName = imageName + '.photom.cat'

	weightName = 'weight'+imageName[5:]
	if os.path.isfile(weightName):
		# imghdr = fits.getheader(imageName)
		# if 'BUNIT' in imghdr:
		#	 scale_fac = imghdr['CONV_FAC']
		# else:
		scale_fac = 1
		weightdata = fits.getdata(weightName)
		weightdata = weightdata / scale_fac**2
		weight_med = np.nanmedian(weightdata)
		weight_std = np.nanstd(weightdata)
		detect_cutoff = weight_med - (weight_std * det_cut)
		with fits.open(weightName, mode='update') as hdu:
			whdr = hdu[0].header
			whdr.set('MEDIAN', weight_med, 'Median of weight image', after='EQUINOX')
			whdr.set('DET_CUT', detect_cutoff, 'Pix value cutoff for source detection', after='MEDIAN')
			hdu.close()
		try:
			print('Including weight map!')
			command = ('sex %s -c %s -CATALOG_NAME %s -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_THRESH %s -WEIGHT_IMAGE %s -PARAMETERS_NAME %s' %
					   (imageName, configFile, catalogName, detect_cutoff, weightName, paramName))
			# print('Executing command: %s' % command)
			subprocess.run(command.split(), check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
		except subprocess.CalledProcessError as err:
			print('Could not run sextractor with exit error %s' % err)
	else:
		try:
			command = ('sex %s -c %s -CATALOG_NAME %s -PARAMETERS_NAME %s' %
					   (imageName, configFile, catalogName, paramName))
			print('Executing command: %s' % command)
			subprocess.run(command.split(), check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
		except subprocess.CalledProcessError as err:
			print('Could not run sextractor with exit error %s' % err)
	return catalogName

def sex2(imageName, det_cut, catalogName, psfName):
	# dynamic aperture photometry adjustment
	init_cat = photometry.get_table_from_ldac(catalogName)
	if isinstance(init_cat['FLUX_RADIUS'][0], np.ndarray):
		acc_sources = init_cat[(init_cat['FLUX_RADIUS'][:, 0] > 1 / 0.498)]
		hwhm = np.nanmedian(acc_sources['FLUX_RADIUS'][:, 0])
		fwhm = 2 * hwhm
	else:
		acc_sources = init_cat[(init_cat['FLUX_RADIUS'][:, 0] > 1 / 0.498)]
		fwhm = np.nanmedian(acc_sources['FLUX_RADIUS'])

	aper_arr = [round(1 * fwhm,2), round(1.25 * fwhm,2), round(1.5 * fwhm,2), round(1.75 * fwhm,2), round(2 * fwhm,2)]
	aper_str = f'-PHOT_APERTURES {aper_arr[0]},{aper_arr[1]},{aper_arr[2]},{aper_arr[3]},{aper_arr[4]}'
	# print(f'Photometric apertures (pix) adjusted dynamically to seeing: {aper_arr}')
	with fits.open(imageName, mode='update') as hdul:
		hdr = hdul[0].header
		hdr.set('APERS', f'{aper_arr[0]},{aper_arr[1]},{aper_arr[2]},{aper_arr[3]},{aper_arr[4]}',
				'adjusted aperture sizes (pix)')
		hdul.close()

	print('Feeding psf model back into sextractor for fitting and flux calculation...')


	psfcatalogName = imageName.replace('.fits', '.photom.cat')
	configFile = photometry.gen_config_file_name('sex2.config')
	psfparamName = photometry.gen_config_file_name('photomPSF.param')
	weightName = 'weight' + imageName[5:]

	if os.path.isfile(weightName):
		# imghdr = fits.getheader(imageName)
		# if 'BUNIT' in imghdr:
		#	 scale_fac = imghdr['CONV_FAC']
		# else:
		scale_fac = 1
		weightdata = fits.getdata(weightName)
		weightdata = weightdata / scale_fac**2
		weight_med = np.nanmedian(weightdata)
		weight_std = np.nanstd(weightdata)
		detect_cutoff = weight_med - (weight_std * det_cut)
		try:
			# We are supplying SExtactor with the PSF model with the PSF_NAME option
			command = (f'sex {imageName} -c {configFile} -CATALOG_NAME {psfcatalogName} -WEIGHT_TYPE MAP_WEIGHT '
					   f'-WEIGHT_THRESH {detect_cutoff} -WEIGHT_IMAGE {weightName} -PSF_NAME {psfName} '
					   f'-PARAMETERS_NAME {psfparamName} {aper_str}') # -DET_THRESH X can change configs here (tune sextractor)
			print("Executing command: %s" % command)
			subprocess.run(command.split(), check=True) #, stdout=None, stderr=None)
		except subprocess.CalledProcessError as err:
			print('Could not run sextractor with exit error %s' % err)
			print('Is there a problem with the sextractor configs or PSF model? Recommend temporarily removing '
				  '"stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL" from this subprocess command to investigate.')
	else:
		try:
			# We are supplying SExtactor with the PSF model with the PSF_NAME option
			command = (f'sex {imageName} -c {configFile} -CATALOG_NAME {psfcatalogName} -PSF_NAME {psfName} '
					   f'-PARAMETERS_NAME {psfparamName} {aper_str}')
			# print("Executing command: %s" % command)
			subprocess.run(command.split(), check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
		except subprocess.CalledProcessError as err:
			print('Could not run sextractor with exit error %s' % err)
			print('Is there a problem with the sextractor configs or PSF model? Recommend temporarily removing '
				  '"stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL" from this subprocess command to investigate.')
	return psfcatalogName

def get_reg(imageName):
	dir_name = os.path.dirname(imageName)
	reg_files = glob.glob(f"{dir_name}/GRB_query_thresh.reg")
	plt_primeregs_all = []

	for reg_file in reg_files:
		primeregs = open(reg_file, 'r')
		plt_primeregs = []
		primeallregs = [reg for reg in primeregs if reg != 'fk5\n']
		for reg in primeallregs:
			nums = re.findall(r'[-+]?\d*\.?\d+', reg)
			srcra = float(nums[0])
			srcdec = float(nums[1])
			srcrad = float(nums[2])
			srccoords = SkyCoord(ra=[srcra], dec=[srcdec], frame='icrs', unit='degree')
			srcreg = CircleSkyRegion(center=srccoords[0], radius=Angle(srcrad, unit='arcsec'))
			plt_primeregs.append(srcreg)
			plt_primeregs_dict = {reg_file:plt_primeregs}
		plt_primeregs_all.append(plt_primeregs_dict)
	return plt_primeregs_all

def grb_cutout(imageName, GRBcoords, photoDistThresh, grb_thresh, name_ext, regprimename=None, regsurvname=None):
		imgdata = fits.getdata(imageName)
		img = fits.open(imageName)
		head = img[0].header
		w = WCS(head)
		grb_name = f"GRB_{name_ext}"
		# savename = '%s_%s_Cutout_%s_%s' % (grb_name, band, survey, num)
		savename = '%s_Cutout' % grb_name
		threshname = '%s_query_thresh.reg' % grb_name

		size = 4 * photoDistThresh * u.arcsec  # smaller radius?
		try:
			cutout = Cutout2D(imgdata, GRBcoords, size, wcs=w, copy=True)
			region = CircleSkyRegion(center=GRBcoords[0], radius=Angle(grb_thresh, unit='arcsec'))
			pix_region = region.to_pixel(cutout.wcs)
		except astropy.nddata.utils.NoOverlapError:
			print(' Area of GRB threshold not found within image, cannot generate cutout!')
			return savename, threshname

		mean, median, sigma_cut = sigma_clipped_stats(cutout.data)
		plt.figure(figsize=(8, 8))
		plt.imshow(cutout.data, vmin=median - 3 * sigma_cut, vmax=median + 3 * sigma_cut, origin='lower', cmap='viridis')

		regs_all = get_reg(imageName)
		colors = ['red', 'cyan', 'green']
		for i, reg_dict in enumerate(regs_all):
			for file, reg_list in reg_dict.items():
				for reg in reg_list:
					pix_reg = reg.to_pixel(cutout.wcs)
					pix_reg.plot(color=colors[i], ls='-', label=os.path.basename(file))

		handles, labels = plt.gca().get_legend_handles_labels()
		seen = set()
		filtered_handles = []
		filtered_labels = []
		for h, l in zip(handles, labels):
			if l not in seen:
				filtered_handles.append(h)
				filtered_labels.append(l)
				seen.add(l)
		plt.legend(filtered_handles, filtered_labels, loc='best')
		
		plt.savefig(os.path.dirname(imageName) + '/' + savename + '.png', dpi=300)
		plt.clf()

		fits.writeto(savename + '.fits', cutout.data, cutout.wcs.to_header(), overwrite=True)
		return savename, threshname


def multi_epoch_astrom(base_epoch_path, matching_epoch_path):
	"""
	Runs sextractor on base epoch stacked image, uses the resulting .cat file as a supplier to
	scamp, which is run on the other epoch.  This is in order to match astrometry between both.
	The matching epoch image is copied to the base epoch directory, and all processes are run there.

	Parameters
	----------
	base_epoch_path: str
		Full filepath to epoch you want to base the astrometry on
	matching_epoch_path: str
		Full filepath to epoch you want to match the base epoch astrometry to
	"""

	base_epoch_dir, base_epoch_name = os.path.split(base_epoch_path)
	os.chdir(base_epoch_dir)
	base_epoch_hdr = fits.getheader(base_epoch_path)
	# if len(base_epoch_hdr['FILTER2']) > 1:
	#	 band = 'Z'
	# else:
	#	 band = base_epoch_hdr['FILTER2']

	match_epoch_dir, match_epoch_name = os.path.split(matching_epoch_path)
	# match_epoch_new_path = os.path.join(base_epoch_dir, match_epoch_name)
	# shutil.copyfile(matching_epoch_path, match_epoch_new_path)
	match_epoch_hdr = fits.getheader(matching_epoch_path)

	base_chip = base_epoch_hdr['CHIP']
	match_chip = match_epoch_hdr['CHIP']
	print(f'Sextracting base epoch: {base_epoch_name}...')
	base_cat_path = swarp_sx(imgpath=base_epoch_path, chip=base_chip)
	print(f'Sextracting matching epoch: {match_epoch_name}...')
	match_cat_path = swarp_sx(imgpath=matching_epoch_path, chip=match_chip)
	multi_epoch_scamp(input_epoch_cat_path=match_cat_path, base_epoch_cat_path=base_cat_path)
	swarp_missfits(imgpath=matching_epoch_path, chip=match_chip)
	return base_epoch_path, matching_epoch_path




def display(file, png=False):
	if png:
		# Read the image data into a NumPy array
		plt.style.use('default')
		cutout = mpimg.imread(file)
		plt.imshow(cutout)

	else:
		plt.style.use(astropy_mpl_style)
	
		image_data = fits.getdata(file , ext=0)
		# print("sum of pixels:",sum(image_data))
		zscale = ZScaleInterval()
		scaled_data = zscale(image_data)
		plt.figure()
		plt.grid(False)
		plt.imshow(scaled_data, cmap='gray')
		plt.colorbar()
	plt.show()

def make_cutout(img, source_ra, source_dec, png=True, display_file=True, name_ext=""):
	if source_ra and source_dec: 
		deci_sky_coords = SkyCoord(ra=[source_ra], dec=[source_dec], frame='icrs', unit='degree')
	
		savename, threshname =  grb_cutout(img, deci_sky_coords, photoDistThresh=4.0, grb_thresh=1.0, name_ext=name_ext) 
		filename = f"{os.path.dirname(img)}/{savename}.png"
		print("made cutout: ",filename)

		if display_file:
			display(filename, png)
		return savename, threshname



def forced_photometry(ra, dec, imageName):

	data = fits.getdata(imageName)
	img = fits.open(imageName)
	header = img[0].header
	w = WCS(header)
	# print(header.keys)

	GRBcoords = SkyCoord(ra=[ra], dec=[dec], frame='icrs', unit='degree')
	aper = SkyCircularAperture(GRBcoords, 1.0*u.arcsec)

	table = aperture_photometry(data, aper, wcs=w)
	flux_counts = table['aperture_sum'][0]
	
	zp = header.get('ZP_APER')

	mag_ab = zp - 2.5 * np.log10(flux_counts)

	return mag_ab


def distance(ra1, dec1, ra2, dec2):
	# print("ra1, dec1, ra2, dec2:", ra1, dec1, ra2, dec2)
	c1 = SkyCoord(ra1, dec1, frame='icrs', unit='deg')
	c2 = SkyCoord(ra2, dec2, frame='icrs', unit='deg')
	ang_sep = c1.separation(c2)
	# print(ang_sep, ang_sep.arcsec)
	
	ang_sep = round(ang_sep.arcsec, 3)
	return ang_sep



	

