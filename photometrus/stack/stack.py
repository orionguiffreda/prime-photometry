"""
Stacks astrometrically calibrated files using swarp
"""
import datetime
#%%
import os
import shutil
import argparse
import subprocess
import sys
from astropy.io import fits
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
import numpy as np
import fnmatch
import matplotlib.pyplot as plt
import pandas as pd

# sys.path.insert(0,'C:\PycharmProjects\prime-photometry\photometrus')
from photometrus.astrom.astrometry import sextract, scamp
from photometrus.photometry.photometry import photometry
from photometrus.settings import gen_config_file_name, auto_bulge_detect, gen_mask_file_name
from photometrus.utils.utils import combine_header_and_fits, remove_wcs_headers

from astropy.stats import SigmaClip
from photutils.background import Background2D, SExtractorBackground

#%%


def get_astrom_files(imgdir):
	image_fnames = [os.path.join(imgdir, f) for f in sorted(os.listdir(imgdir)) if
					f.endswith('.flat.fits') or f.endswith('.flat.new')]
	return image_fnames


def datetime_to_julian_date(datetime):
	ts = Time(datetime)
	jd = ts.jd1 + ts.jd2
	return jd


def julian_date_to_modified_julian_date(julian_date):
	return julian_date - 2400000.5


def datetime_to_modified_julian_date(datetime):
	jd = datetime_to_julian_date(datetime)
	return julian_date_to_modified_julian_date(jd)


def get_time_stack_header(header_df):
	dt_fmt = '%Y-%m-%dT%H:%M:%S.%f'
	return_header = fits.Header()
	start_frame = header_df.iloc[0]
	end_frame = header_df.iloc[-1]
	start_time = datetime.datetime.strptime(start_frame['DATE-BEG'], dt_fmt)
	end_time = datetime.datetime.strptime(end_frame['DATE-END'], dt_fmt) + datetime.timedelta(seconds=end_frame['EXPTIME'])
	return_header.set('DATE-BEG', start_time.strftime(dt_fmt), 'UT DT at start of 1st exp in stack')
	return_header.set('DATE-END', end_time.strftime(dt_fmt), 'UT DT at end of last exp in stack')
	return_header.set('MJD-BEG', datetime_to_modified_julian_date(start_time), 'modified jd at start of 1st exp')
	return_header.set('MJD-END', datetime_to_modified_julian_date(end_time), 'modified jd at end of last exp')
	return_header.set('JD-BEG', datetime_to_julian_date(start_time), 'julian date at start of 1st exp')
	return_header.set('JD-END', datetime_to_julian_date(end_time), 'julian date at end of last exp')
	return_header.set('EXPSTART', datetime_to_modified_julian_date(start_time), 'modified jd at start of 1st exp')
	return_header.set('EXPEND', datetime_to_modified_julian_date(end_time), 'modified jd at end of last exp')
	return_header.set('EXPMID', (return_header['EXPSTART'] + return_header['EXPEND']) / 2, 'modified jd at middle of stack exp')
	return_header.set('TELAPSE', (return_header['EXPEND'] - return_header['EXPEND']) * 24 * 3600, '[s] time elapsed during stack')
	return_header.set('EXPERR', return_header['TELAPSE'] / 2, '[s] time error around EXPMID for stack')
	return return_header


def gen_stack_header(image_fnames):
	# image_fnames = get_astrom_files(imgdir)
	image_fnames.sort()
	all_headers = [remove_wcs_headers(fits.getheader(f)) for f in image_fnames]
	header_template = all_headers[0]
	all_headers_df = pd.DataFrame(all_headers)
	delete_headers = [
		'NAXIS1', 'NAXIS2', 'FRAME', 'COMMENT', 'ISRESET', 'SIZAXIS1', 'SIZAXIS2', 'NCOLS', 'NROWS', 'UTDATE',
		'UTSTART', 'DATE', 'ASDFNAME', 'RSTSAVE', 'SCIWRD1', 'SCIWRD2', 'SCIWRD3', 'SCIWRD4', 'SCIWRD5', 'SCIWRD6',
		'ASCIXPID', 'NSCIHEAD', 'ABLKLEN', 'ABLKCNT', 'AASICID', 'AHDRLEN', 'ASDPCNT', 'AFF', 'AEXPVIDL', 'ASCIFRM',
		'SATURATE', 'GAIN', 'SLICE', 'BLOCKID', 'DITHPH', 'DITH_IDX', 'INT', 'USEDNUMS', 'USEDNUML', 'CHECKSUM',
		'DATASUM',
	]
	mean_headers = [
		'NRESETS', 'ASICADDR', 'ASICINDX', 'DETECTOR', 'ASICSN', 'FPAPOS', 'MACIESN', 'VRESET', 'CELLDRAI', 'VBIASGAT',
		'VBIASPOW', 'SUB', 'DSUB', 'DRAIN', 'VDDA', 'VDD', 'GND', 'GNDA', 'VREF', 'CHIP', 'RAOFF', 'DECOFF', 'ROTOFF',
		'DITHRAD', 'DITH_REP', 'TEMPASIC', 'TEMPPLAT', 'TEMPSTRP', 'TEMPDET', 'TEMPMOTI', 'TEMPMOTO', 'PRESSURE',
		'ALT', 'AZI', 'RA-D', 'DEC-D', 'TSDOME', 'FOCUS', 'ROT', 'TIMEOFF', 'TSSECZ', 'ETMPCCP', 'ETMPCC1', 'ETMPCCC1',
		'ETMPCC2', 'ETMPCCC2', 'ETMPDECK', 'ETMPHK', 'ETMPMAC', 'ETMPMANI', 'ETMPMANO', 'ETMPPDU', 'ETMPPRES',
		'TEMPHED1', 'PWRHED1C', 'PWRHED1M', 'TEMPREJ1', 'CRYO1STA', 'TMPHED1T', 'TEMPHED2', 'PWRHED2C', 'PWRHED2M',
		'TEMPREJ2', 'CRYO2STA', 'TMPHD2T', 'SKY_FAC', 'AIRMASS', 'X_SHIFT', 'Y_SHIFT',
	]
	sum_headers = ['NFRAMES', 'EXPTIME', 'EXPTIMEE', 'EXPTIMEC', 'NINT', 'DITH_TOT']
	first_headers = [
		'FRTIME', 'TFRAME', 'REDXMODE', 'REFOUT', 'NOUTPUTS', 'AMPMODE', 'MACIEINT', 'TIMEUNIT', 'DEINTERL', 'RSTTYPE',
		'LODFIL0', 'LODFIL1', 'LODFIL2', 'LODFIL3', 'LODFIL4', 'LODFIL5', 'LODFIL6', 'LODFIL7', 'LODFIL8', 'LODFIL9',
		'LODFIL10', 'RMVSCIHD', 'PIXSCALE', 'INSTRUME', 'LATITUDE', 'LONGITUD', 'ALTITUDE', 'BINNING', 'BINX', 'BINY',
		'WAVELENG', 'PIXSIZE', 'OBSERVER', 'OBJNAME', 'OBJTYPE', 'DITHTYP', 'COMMENT1', 'COMMENT2', 'DEC', 'RA',
		'FILTER1', 'FILTER2', 'SKY_FILE', 'COER', 'COED', 'DARKLIM', 'SPBIAS', 'SATULIM',
	]
	last_headers = []
	time_headers = get_time_stack_header(all_headers_df)
	for k in delete_headers:
		try:
			del header_template[k]
		except KeyError:
			print('Could not find header key "%s" skipping delete' % k)
	for k in mean_headers:
		try:
			header_template.set(k, np.mean(all_headers_df[k]), 'Stack mean of '+ header_template.comments[k])
		except KeyError:
			print('Could not find header key "%s" skipping mean' % k)
		except TypeError:
			first_headers.append(k)
	for k in sum_headers:
		try:
			header_template.set(k, np.sum(all_headers_df[k]), 'Stack sum of ' + header_template.comments[k])
		except KeyError:
			print('Could not find header key "%s" skipping sum' % k)
		except TypeError:
			first_headers.append(k)
	for k in first_headers:
		try:
			header_template[k] = all_headers_df[k].tolist()[0]
		except KeyError:
			print('Could not find header key "%s" skipping first' % k)
	for k in last_headers:
		try:
			header_template[k] = all_headers_df[k].tolist()[-1]
		except KeyError:
			print('Could not find header key "%s" skipping last' % k)
	header_template.update(time_headers)
	return header_template



def badpixmask(subpath, chip):
	old_storage_dir = os.path.join(subpath, 'no_mask')
	if not os.path.exists(old_storage_dir):
		os.mkdir(old_storage_dir)
	mask = gen_mask_file_name('badpixmask_c%i.fits' % chip)
	badmask = fits.getdata(mask)
	badmask = badmask.astype(bool)
	for i in sorted(os.listdir(subpath)):
		if i.endswith('.sky.flat.fits') or i.endswith('.sky.flat.new'):
			shutil.move(os.path.join(subpath, i), os.path.join(old_storage_dir, i))
	print('Applying bad pixel mask...')
	for i in sorted(os.listdir(old_storage_dir)):
		img = fits.open(os.path.join(old_storage_dir, i))
		hdr = img[0].header
		data = img[0].data
		data[~badmask] = 65000
		fits.writeto(os.path.join(subpath, i), data, hdr)
	return old_storage_dir


def astromfin(directory, chip):
	if not chip:
		raise ValueError('Specify a chip when using this functionality!')
	print('Re-running astrometry on swarped image! Running sextractor...')
	catpath = swarp_sx(directory, chip)
	print('Applying 4th order scamp fit to stacked image...')
	scamp(directory, swarpcat=catpath)
	print('Combining scamp .head and stacked image...')
	swarp_missfits(directory, chip)
	try:
		os.remove(catpath)
	except Exception as e:
		print(f"Error removing file: {catpath} - {e}")
	print('Absolute astrometry complete!')


def astrom_check(imgdir):
	"""
	Check to determine if all images are w/in the same area in the sky (2x dither rad) before attempting stacking.

	Parameters
	----------
	imgdir: str
		Directory where input images are stored
	"""

	image_fnames = get_astrom_files(imgdir)

	image_hdrs = [fits.getheader(img) for img in image_fnames]

	image_ras = [(float(hdr['CRVAL1'])) for hdr in image_hdrs]
	image_decs = [(float(hdr['CRVAL2'])) for hdr in image_hdrs]
	image_coords = SkyCoord(ra=image_ras, dec=image_decs, frame='icrs', unit='degree')

	med_ra = np.nanmedian(image_ras)
	med_dec = np.nanmedian(image_decs)
	med_coords = SkyCoord(ra=med_ra, dec=med_dec, frame='icrs', unit='degree')

	med_dith_rad = np.nanmedian([(float(hdr['DITHRAD'])) for hdr in image_hdrs])
	acc_radius = (2 * med_dith_rad) * u.arcsec

	seps = image_coords.separation(med_coords)

	acc_mask = seps < acc_radius

	if not acc_mask.all():
		print(' Pre-stacking astrom check shows 1 or more images are *NOT* w/in acceptable area!')
		bad_idxs = np.where(~acc_mask)[0]
		bad_imgs = [item for index, item in enumerate(image_fnames) if index in bad_idxs]
		bad_img_names = [os.path.split(img)[1] for img in bad_imgs]
		print(f' Recommend checking quality / astrometry on offending images: {bad_img_names}')

		if len(bad_imgs) > len(image_fnames) // 2:
			raise Exception(f'*WARNING* {len(bad_imgs)}/{len(image_fnames)} images (> 1/2 total # of images) are NOT '
							f'in acceptable area, examine images!  Is there an issue with image acquisition tracking '
							f'or astrometry?')

		print(' Renaming offending images to avoid stacking issues...')
		for img_name in bad_imgs:
			os.rename(img_name, img_name.replace('.flat.','.flat.EXCL.'))


def make_weight_map(image_adu, chip, use_badmask=True, coadd=True, bkg_percent=1.3):
	"""
	construct weight map in ADU
	"""
	
	data_adu = fits.getdata(image_adu)
	# fill_value = 10000
	# data_adu = np.nan_to_num(data_adu, nan=fill_value, posinf=fill_value, neginf = -1*fill_value) # fill NaNs to high noise default


	hdr = fits.getheader(image_adu)
	new_path = f"{image_adu}.weightmap.fits"

	# gain = 1.8 # e- / ADU
	# # shot noise from electron statistics
	# sky_file = os.path.dirname(os.path.dirname(image_adu)) +'/sky/'+ hdr['SKY_FILE']
	# raw_adu = data_adu + fits.getdata(sky_file) * hdr["SKY_FAC"]
	# poisson_noise = np.nan_to_num(np.sqrt((raw_adu) / gain), 0) # sigma_ADU = sigma_e- / gain = sqrt(ADU * g) / g = sqrt(ADU / g)

	# weight_map = fits.getdata(os.path.dirname(image_adu)+'/weight'+os.path.basename(image_adu)[5:])
	bad_mask = data_adu == 0

	# background noise (read noise, flat fielding, sky subtraction, nonlinearity correction)
	bkg_rms = Background2D(data_adu, 
						   box_size=32,
						   filter_size=5, 
						   sigma_clip=SigmaClip(sigma=3), 
						   bkg_estimator=SExtractorBackground(),
						   mask=bad_mask).background_rms # uses source extractor background estimator


	fits.writeto(os.path.dirname(image_adu)+'/rms_coadd.fits', bkg_rms, hdr, overwrite=True)
	print("wrote to", os.path.dirname(image_adu)+'/rms_coadd.fits')

	bkg_rms[bkg_rms==0] = np.percentile(bkg_rms, 99)


	# print(f"Zero pixels: {np.sum(bkg_rms == 0)}")
	# print(f"bkg_rms range: {np.min(bkg_rms):.4f} - {np.max(bkg_rms):.4f}")
	# print(f"bkg_rms percentiles: {np.percentile(bkg_rms, 99)}")

	sigma_2 = bkg_rms ** 2 #+ poisson_noise**2
	weight = 1/sigma_2

	print("MEANS data:", np.mean(data_adu), ", bkg:", np.mean(bkg_rms), ", noise:",np.mean(sigma_2), ", weight:", np.mean(weight))


	if not coadd:
		mask = gen_mask_file_name('/home/alex/PycharmProjects/prime-photometry/photometrus/weightmaps/badpixmask_c%i.fits' % chip)
		badmask = fits.getdata(mask)
		badmask = np.fliplr(badmask) # noticing the maps were flipped 
		badmask = badmask.astype(bool) # 1 is good, use
		weight[~badmask] = 0
	
		

	 # set any weights to zero which are below a certain value
	
	print(np.median(weight), np.std(weight))

	bkg_thresh = bkg_percent * np.median(bkg_rms)
	noisy_areas = bkg_rms > bkg_percent * bkg_thresh # mask noisy areas 
	print("setting rms areas above", bkg_thresh, "to zero")
	weight[noisy_areas] = 0

	# print("weight after mask:", np.mean(weight))

	fits.writeto(new_path, weight, hdr, overwrite=True)
	return new_path
	
def make_weight_maps(images, chip):
	"""
	Make weight maps for input list of dither images, based on poisson noise, background (read) noise, and bad pixel maps

	Parameters
	----------
	images: list
		list of filepaths for input images
	chip: int
		chip for observations
	"""

	print("Making weight map from poisson noise, bkg noise, bad px mask, ...")
	# maybe im applying this rotated 90 deg because the corners look super noisy 
	weight_maps = ""
	for image_adu in images:
		weight_maps+= make_weight_map(image_adu, chip)+","
	return weight_maps

def remove_small_weights(weight_file):
	weight = fits.getdata(weight_file)
	print(np.median(weight)/1.3)
	weight_cutoff = np.median(weight)/1.3
	# print("setting weights below", weight_cutoff, "to zero")
	# weight[weight < weight_cutoff] = 0
	fits.writeto(weight_file, weight, fits.getheader(weight_file), overwrite=True)
	

	
def swarp(imgdir, finout, chip):
	print('SWARP Stacking!')
	image_fnames = get_astrom_files(imgdir)
	image_fnames.sort()

	# header = fits.getheader(image_fnames[-1])
	header = gen_stack_header(image_fnames)
	filter1 = header.get('FILTER1', 'unknown')
	filter2 = header.get('FILTER2', 'unknown')
	ext = os.path.splitext(image_fnames[-1])[1]
	if ext == '.new':
		save_name = 'coadd.{}-{}.{}-{}.C{}.fits'.format(filter1, filter2, image_fnames[0][-23:-15],
														image_fnames[-1][-23:-15], image_fnames[0][-14])
		print(save_name)
		weight_name = 'weight.{}-{}.{}-{}.C{}.fits'.format(filter1, filter2, image_fnames[0][-23:-15],
														image_fnames[-1][-23:-15], image_fnames[0][-14])
	else:
		save_name = 'coadd.{}-{}.{}-{}.C{}.fits'.format(filter1, filter2, image_fnames[0][-24:-16],
														image_fnames[-1][-24:-16], image_fnames[0][-15])
		print(save_name)
		weight_name = 'weight.{}-{}.{}-{}.C{}.fits'.format(filter1, filter2, image_fnames[0][-24:-16],
														image_fnames[-1][-24:-16], image_fnames[0][-15])

	os.chdir(str(finout))
	#save_name = 'coaddastr.fits'
	#weight_name = 'coaddastrweight.fits'



	sw = gen_config_file_name('default.swarp')
	
	use_weight_map = True
	if use_weight_map:
		weight_maps = make_weight_maps(image_fnames, chip)
		print('Generating weight map:', weight_name)

		com = f"swarp {os.path.join(imgdir, '*.flat'+ext)} -c {sw} -IMAGEOUT_NAME {save_name} -WEIGHTOUT_NAME {weight_name} -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE {weight_maps}"

	else:
		com = f"swarp {os.path.join(imgdir, '*.flat'+ext)} -c {sw} -IMAGEOUT_NAME {save_name} -WEIGHTOUT_NAME {weight_name}"


	try:

	# WEIGHT_SUFFIX .weight.fits (instead of weight image)

		
		subprocess.run(com, shell=True, check=True, capture_output=True, text=True)
	except subprocess.CalledProcessError as err:
		print(f"SWARP failed with exit code {err.returncode}")
		print(f"STDERR:\n{err.stderr}")
		print(f"STDOUT:\n{err.stdout}")
	print('Co-added image created, all done!')
	remove_small_weights(weight_name) 
	with fits.open(save_name, mode='update') as fin:
		fin[0].header.update(header)



def swarp_increm(imgdir, finout, im_num=5):
	batch_size = im_num
	files = get_astrom_files(imgdir)
	total_files = len(files)

	stack_files = [os.path.join(finout, f) for f in sorted(os.listdir(finout)) if f.startswith('coadd') and f.endswith('.fits')]
	if len(stack_files) < 1 - round(total_files / batch_size):

		for i in range(0, total_files, batch_size):
			batch = files[:i + batch_size]
			print(batch)
			print('images taken = ', len(batch))

			image_fnames = batch
			header = fits.getheader(image_fnames[-1])
			filter1 = header.get('FILTER1', 'unknown')
			filter2 = header.get('FILTER2', 'unknown')
			ext = os.path.splitext(image_fnames[-1])[1]
			if ext == '.new':
				save_name = 'coadd.{}-{}.{}-{}.C{}.fits'.format(filter1, filter2, image_fnames[0][-23:-15],
																image_fnames[-1][-23:-15], image_fnames[0][-14])
				print(save_name)
				weight_name = 'weight.{}-{}.{}-{}.C{}.fits'.format(filter1, filter2, image_fnames[0][-23:-15],
																   image_fnames[-1][-23:-15], image_fnames[0][-14])
			elif ext == '.fits':
				save_name = 'coadd.{}-{}.{}-{}.C{}.fits'.format(filter1, filter2, image_fnames[0][-24:-16],
																image_fnames[-1][-24:-16], image_fnames[0][-15])
				print(save_name)
				weight_name = 'weight.{}-{}.{}-{}.C{}.fits'.format(filter1, filter2, image_fnames[0][-24:-16],
																   image_fnames[-1][-24:-16], image_fnames[0][-15])
			os.chdir(str(finout))
			sw = gen_config_file_name('default.swarp')
			#save_name = 'coaddastr.fits'
			#weight_name = 'coaddastrweight.fits'
			image_list = (',').join(image_fnames)
			com = ["swarp ", image_list, ' -c '+sw
				   , ' -IMAGEOUT_NAME '+save_name, ' -WEIGHTOUT_NAME '+weight_name]
			s0 = ''
			com = s0.join(com)
			out = subprocess.Popen([com], shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
			out.wait()
			
			print('Co-added image created, all done!')

	stack_files = [os.path.join(finout, f) for f in sorted(os.listdir(finout)) if f.startswith('coadd') and f.endswith('.fits')]
	total_stacks = len(stack_files)
	band = fits.getheader(stack_files[0])['FILTER2']
	if len(band) > 1:
		band = 'Z'

	ecsv_files = [f for f in sorted(os.listdir(finout)) if f.startswith('coadd') and f.endswith('.ecsv')]

	lim_mags = []
	for stack in stack_files:
		if len(ecsv_files) < 1 - round(total_files / batch_size):
			print(f'\nRunning photometry on {stack}\n')
			photometry(full_filename=stack, band=band)
		hdr = fits.getheader(stack)
		lim_mag = hdr['LM_AUTO']
		lim_mags.append(lim_mag)

	exptime = fits.getheader(files[0])['EXPTIMEC']
	single_stack_exp = exptime * im_num
	max_exp = (total_stacks * single_stack_exp)
	exptimes_axis = np.arange(single_stack_exp, max_exp + single_stack_exp, single_stack_exp)

	print('\nGenerating incremental lim mag plot!')
	plt.figure(1, figsize=(10,8))
	plt.plot(exptimes_axis, lim_mags, 'ro')
	plt.grid()
	# plt.yscale('log')
	# plt.xscale('log')
	plt.xticks(exptimes_axis, rotation=45, ha='right')
	plt.xlabel('Exposure Time (s)')
	plt.ylabel(f'{band}MAG_AUTO Limiting Magnitude (AB)')
	plt.title('PRIME Limiting Magnitude Evolution')
	plt.savefig(os.path.join(finout, 'lim_mag_increm.png'), dpi=200)
	plt.clf()
	print('Generated!')

	plt.close('all')


def swarp_alt(imgdir, imout):
	image_fnames = get_astrom_files(imgdir)
	image_fnames.sort()
	image_fnames_1 = image_fnames[::2]
	image_fnames_2 = image_fnames[1::2]
	print('1st stack = ', len(image_fnames_1))
	print('2nd stack = ', len(image_fnames_2))

	header = fits.getheader(image_fnames[-1])
	filter1 = header.get('FILTER1', 'unknown')
	filter2 = header.get('FILTER2', 'unknown')
	ext = os.path.splitext(image_fnames_1[-1])[1]
	if ext == '.new':
		save_name1 = 'coadd.{}-{}.{}-{}.C{}.fits'.format(filter1, filter2, image_fnames_1[0][-23:-15],
														image_fnames_1[-1][-23:-15], image_fnames_1[0][-14])
		print(save_name1)
		weight_name1 = 'weight.{}-{}.{}-{}.C{}.fits'.format(filter1, filter2, image_fnames_1[0][-23:-15],
														   image_fnames_1[-1][-23:-15], image_fnames_1[0][-14])

		save_name2 = 'coadd.{}-{}.{}-{}.C{}.fits'.format(filter1, filter2, image_fnames_2[0][-23:-15],
														image_fnames_2[-1][-23:-15], image_fnames_2[0][-14])
		print(save_name2)
		weight_name2 = 'weight.{}-{}.{}-{}.C{}.fits'.format(filter1, filter2, image_fnames_2[0][-23:-15],
														   image_fnames_2[-1][-23:-15], image_fnames_2[0][-14])
	elif ext == '.fits':
		save_name1 = 'coadd.{}-{}.{}-{}.C{}.fits'.format(filter1, filter2, image_fnames_1[0][-24:-16],
														image_fnames_1[-1][-24:-16], image_fnames_1[0][-15])
		print(save_name1)
		weight_name1 = 'weight.{}-{}.{}-{}.C{}.fits'.format(filter1, filter2, image_fnames_1[0][-24:-16],
														   image_fnames_1[-1][-24:-16], image_fnames_1[0][-15])

		save_name2 = 'coadd.{}-{}.{}-{}.C{}.fits'.format(filter1, filter2, image_fnames_2[0][-24:-16],
														image_fnames_2[-1][-24:-16], image_fnames_2[0][-15])
		print(save_name2)
		weight_name2 = 'weight.{}-{}.{}-{}.C{}.fits'.format(filter1, filter2, image_fnames_2[0][-24:-16],
														   image_fnames_2[-1][-24:-16], image_fnames_2[0][-15])
	os.chdir(str(imout))
	sw = gen_config_file_name('default.swarp')
	# save_name = 'coaddastr.fits'
	# weight_name = 'coaddastrweight.fits'
	image_list1 = (',').join(image_fnames_1)
	image_list2 = (',').join(image_fnames_2)
	com = ["swarp ", image_list1, ' -c ' + sw
		, ' -IMAGEOUT_NAME ' + save_name1, ' -WEIGHTOUT_NAME ' + weight_name1]
	com = ''.join(com)
	out = subprocess.Popen([com], shell=True)
	out.wait()
	print('First set co-added image created, all done!')

	com2 = ["swarp ", image_list2, ' -c ' + sw
		, ' -IMAGEOUT_NAME ' + save_name2, ' -WEIGHTOUT_NAME ' + weight_name2]
	com2 = ''.join(com2)
	out2 = subprocess.Popen([com2], shell=True)
	out2.wait()
	print('2nd set co-added image created, all done!')


def swarp_sx(imgpath, chip):
	# imgpath can be either full file path to stack or directory where stack is in

	if os.path.isdir(imgpath):
		os.chdir(str(imgpath))
		bulge = auto_bulge_detect(imgpath)
		if bulge:
			sx = gen_config_file_name('bulge_new.config')
			ap = gen_config_file_name('tempsource.param')
		else:
			sx = gen_config_file_name('sex.config')
			ap = gen_config_file_name('astrom_coadd.param')

		stackimg = [f for f in sorted(os.listdir(imgpath)) if fnmatch.fnmatch(f, 'coadd.*.C%i.fits' % chip)]
		coaddimg = stackimg[0]
		catname = coaddimg.replace('.fits', '.cat')
		catpath = os.path.join(imgpath, catname)
		weightname = 'weight' + coaddimg[5:]
	elif os.path.isfile(imgpath):
		coaddimgdir, coaddimg = os.path.split(imgpath)
		os.chdir(str(coaddimgdir))
		bulge = auto_bulge_detect(coaddimgdir)
		if bulge:
			sx = gen_config_file_name('bulge_new.config')
			ap = gen_config_file_name('tempsource.param')
		else:
			sx = gen_config_file_name('sex.config')
			ap = gen_config_file_name('astrom_coadd.param')

		catname = coaddimg.replace('.fits', '.cat')
		catpath = os.path.join(coaddimgdir, catname)
		weightname = 'weight' + coaddimg[5:]
	else:
		print('Must specify a directory to stacked image or full file path to stacked image! '
			  'Cannot continue with stack image sextraction!')
		catname = catpath = weightname = None
	if catname:
		if os.path.isfile(weightname):
			print('Including weight map!')
			command = ('sex %s -c %s -CATALOG_NAME %s -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE %s -PARAMETERS_NAME %s' %
					   (coaddimg, sx, catname, weightname, ap))
			# print('Executing command: %s' % command)
			subprocess.run(command.split(), check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
		else:
			command = ('sex %s -c %s -CATALOG_NAME %s -PARAMETERS_NAME %s' %
					   (coaddimg, sx, catname, ap))
			# print('Executing command: %s' % command)
			subprocess.run(command.split(), check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
		return catpath


def swarp_missfits(imgpath, chip):
	# imgpath can be either full file path to stack or directory where stack is in
	if os.path.isdir(imgpath):
		stackimg = [f for f in sorted(os.listdir(imgpath)) if fnmatch.fnmatch(f, 'coadd.*.C%i.fits' % chip)]
		stackhdr = [f for f in sorted(os.listdir(imgpath)) if fnmatch.fnmatch(f, 'coadd.*.C%i.head' % chip)]
		combine_header_and_fits(stackhdr[0], stackimg[0], remove_header_file=True)
	else:
		stackimg = imgpath
		stackhdr = imgpath.replace('.fits', '.head')
		combine_header_and_fits(stackhdr, stackimg, remove_header_file=True)


#%%


def stack(subpath, stackpath, chip, num=5, no_astrom=False, astrom_only=False, increm=False, alt=False, mosaic=False):
	# if args.mask:
	#	 otherdir = badpixmask(args.parent,args.sub,args.chip)
		# print('removing temp dir...')
		# shutil.rmtree(otherdir)
	if not mosaic:
		astrom_check(subpath)

	if no_astrom:
		swarp(subpath, stackpath)
	elif astrom_only:
		astromfin(stackpath, chip)
	elif increm:
		swarp_increm(subpath, stackpath, num)
	elif alt:
		swarp_alt(subpath, stackpath)
		astromfin(stackpath, chip)
	else:
		if not chip:
			raise ValueError('Remember to specify chip number using default stacking behavior!  It is required for '
							 'absolute astrometry check!')
		# badpixmask(subpath=subpath, chip=chip)
		swarp(subpath, stackpath, chip)
		astromfin(stackpath, chip)


def main():
	parser = argparse.ArgumentParser(description='Runs swarp to stack imgs, then reruns astrometry for improved wcs *NOTE* if using -mask flag, do not keyboard interrupt')
	# parser.add_argument('-mask', action='store_true', help='optional flag, use if you want to utilize a bad pixel mask')
	parser.add_argument('-chip', type=int, help='[int] Detector chip number', default=None)
	parser.add_argument('-no_astrom', action='store_true', help='optional flag, use if you just want the swarped image, not the image with improved astrometry')
	parser.add_argument('-astrom_only', action='store_true',
						help='optional flag, use if you already have the swarped image, but want improved astrometry')
	parser.add_argument('-increm', action='store_true',
						help='optional flag, use if you want to generate a stacked img from increments of images, ex. 5 stack, then 10 stack, etc.')
	parser.add_argument('-alt', action='store_true',
						help='create 2 stacked images from 1 set of data, alternating images used')
	parser.add_argument('-mosaic', action='store_true',
						help='Use this flag if you are attempting to make a large mosaic, will disable the default '
							 'image location screening')
	parser.add_argument('-sub', type=str, help='[str] Processed images path')
	parser.add_argument('-stack', type=str, help='[str] Output stacked image path')
	parser.add_argument('-num', type=int, help='*USE ONLY W/ -INCREM* [int] # of individual imgs to increment by, default=5', default=5)
	#parser.add_argument('-parent', type=str, help='*USE ONLY W/ -MASK FLAG* [str] Parent directory where all img folders are stored', default=None)
	args, unknown = parser.parse_known_args()

	stack(args.sub, args.stack, args.chip, args.num, args.no_astrom, args.astrom_only, args.increm, args.alt, args.mosaic)


if __name__ == "__main__":
	main()
