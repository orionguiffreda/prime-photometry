import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
from pathlib import Path
import ast
import math
import random 
import glob

from astropy.table import Table
from astropy.io import fits

import photometrus_utils as util
import sfft_util
import generate_image_list as gen_img
import assess_processed_data_util as gen_util
import subprocess

def run_bog_photometry(file, band, ra, dec):
	"""
	Make a bogus directory if it doesnt exist and run photometry methods on bogus source
	"""
	stack_path = os.path.dirname(file)
	coadd_file = os.path.basename(file)
	header = fits.getheader(file)
	
	bogus_dir = f"{stack_path}/bogus"
	try:
		os.listdir(bogus_dir)
	except IOError:
		os.mkdir(bogus_dir)
		os.system(f'cp {file} {bogus_dir}/{coadd_file}')
		files = os.listdir(bogus_dir)
		print(files)


	ecsv = gen_img.get_ecsv(bogus_dir)
	if ecsv== "":
		try:

		# bigger grb radius to contain overlapping sources, get mags of those
			print("photometrus", "photometry",
						"-stackpath",f"{bogus_dir}/",
						"-band", band,
						"-chip", f"{header['CHIP']}",
						"-grb_ra", ra,
						"-grb_dec", dec,
						"-grb_radius", "4",)
			
			subprocess.run([
						"photometrus", "photometry",
						"-stackpath",f"{bogus_dir}/",
						"-band", f"{band}",
						"-chip", f"{header['CHIP']}",
						"-grb_ra", f"{ra}",
						"-grb_dec", f"{dec}",
						"-grb_radius", "4",
						])

		except Exception as e:
			print(e)
			
	ecsv = gen_img.get_ecsv(bogus_dir)
	row = {}
	if ecsv!= "":
		df = pd.read_csv(f"{bogus_dir}/{ecsv}", sep='\\s+', comment='#')
		ecsv_data = df.to_dict(orient='records')
		for ecsv_dict in ecsv_data:
			row.update(ecsv_dict)
	
	return row


def limiting_mag(header):
	"""
	Retrieve limiting mag from file header
	"""
	try:
		psflim, autolim = float(header['HIERARCH lim_mag_psf']), float(header['HIERARCH lim_mag_auto'])
	except KeyError as e:
		try:
			psflim, autolim = float(header['LM_PSF']), float(header['LM_AUTO'])
		except KeyError as e2:
			psflim, autolim = 0, 0
	return np.average([psflim, autolim])
	

def make_bogus_df(triplets):
	"""
	Make bogus DataFrame with a randomly selected source in difference image. Enforce boguses to have detection in sci image
	"""
	
	# triplets = triplets[triplets['full_name'].isin(['SN2024ugv'])]
	data = []
	for i, triplet in triplets.iterrows():

		image_keys = ['full_name', 'field', 'truncated_field', 'band', 'chip', 'coadd_path', 'discoverydate', 'EXPTIMEE', 'DITHTYP', 'DITHRAD', 'DITH_REP', 'NINT', 'ROTOFF', 'RA', 'DEC', 'field_date', 'days_since_discovery', 'role', 'ref_coadd_path', 'ref_field', 'ref_chip',  'ref_days_since_discovery', 'objid', 'name_prefix', 'name', 'edge', 'px_from_edge', 'b', 'sfft_diff', 'diff_cat', 'sci_rms', 'ref_rms', 'diff_rms', 'sci_fwhm', 'ref_fwhm', 'sci_cat_len', 'ref_cat_len', 'diff_cat_len']
			
		row = {key: val for key, val in triplet[image_keys].items()}
		
		ref = ast.literal_eval(triplet['ref_coadd_path'])[0]
		sci = triplet['coadd_path']

		try: # missing diff info?
			print(triplet, triplet["sfft_diff"])
			diff_file = triplet["sfft_diff"]
			diff_data = fits.getdata(diff_file)
			diff_cat_name = triplet["diff_cat"]
			diff_cat = Table.read(diff_cat_name, hdu=2)
			cat_good = diff_cat[sfft_util.get_good_mask(diff_cat)]

		except OSError as e:
			print("read file error, Malformed diff? :", e)
			print(diff_file)
			print()
			try:
				glob_str = f"/mnt/photometry/{triplet['full_name']}/*.sfftdiff.fits"
				print(glob_str)
				diff_files = glob.glob(glob_str)
				if len(diff_files)==0:
					print("no diff found?")
					
				diff_file = glob.glob(glob_str)[-1]
				diff_data = fits.getdata(diff_file)
				cat_good = sfft_util.get_catalog(diff_file)
				triplet["sfft_diff"] = diff_file
			except IndexError as e:
				print("missing diff or cat!")
				continue

		cat_good_df = sfft_util.make_df(cat_good)
		band = triplet['band']

		# crop
		xdim, ydim = np.shape(diff_data)
		cat_good_df = sfft_util.crop_mask(cat_good_df, xdim, ydim, crop=1000)
		cat_good_df = cat_good_df[cat_good_df["SNR_WIN"]>10] 

		detected_sci = False
		print("# good sources:",len(cat_good_df))

		ite = 0 
		random_idx_list = random.sample(range(0, len(cat_good_df)), len(cat_good_df))

		# iterate through sources in image randomly 
		while(not detected_sci and ite < len(cat_good_df)):
			source_idx = random_idx_list[ite]
			source = cat_good_df.iloc[source_idx]

			# print()
			# print("NEW SOURCE:", triplet['full_name'])
			# print("idxs", ite, source_idx)

			ra, dec = source['ALPHA_J2000'], source['DELTA_J2000']
			if ra-triplet['ra'] <0.0001 and dec-triplet['declination'] <0.0001:
				print(f"selected same source at ra={ra}, dec={dec}!")
				ite += 1
				continue

			forced_photometry_mag = util.forced_photometry(ra, dec, sci)
			lim_mag = limiting_mag(fits.getheader(sci))

			print("detected mag:",forced_photometry_mag, " lim mag:", round(lim_mag, 3))

			# if nothing is detected via forced photometry, dont try to run photometry
			if math.isnan(forced_photometry_mag) or forced_photometry_mag > lim_mag:
				print("Bogus not detectable in sci")
				print("detected mag:",forced_photometry_mag, " lim mag:", round(lim_mag, 3), "ra:", ra, "dec:", dec)
				ite +=1
				continue
				
			# savename, threshname = util.make_cutout(sci, ra, dec)
			# savename, threshname = util.make_cutout(ref, ra, dec)
			# savename, threshname = util.make_cutout(diff_file, ra, dec)
					
			sci_photometry_row = run_bog_photometry(sci, band, ra, dec)

			if sci_photometry_row == {}: # no sci detection, don't include this example! 
				ite+=1
				continue
			else:
				# print("detected sci:", sci_photometry_row)
				detected_sci = True

		ref_photometry_row = run_bog_photometry(ref, band, ra, dec)
		if ref_photometry_row == {}:
			ref_photometry_row = {f"{key}_ref":0 for key, val in sci_photometry_row.items()}
		else:
			ref_photometry_row = {f"{key}_ref":val for key, val in ref_photometry_row.items()}

		# print(sci_photometry_row)

		try:
			row.update({'ra': sci_photometry_row['RA'], 'declination':sci_photometry_row['DEC']})
		except KeyError as e:
			print("Error with non detection row:", e)
			print(sci_photometry_row)
			print()
			continue

		
		row.update(sci_photometry_row)
		row.update(ref_photometry_row)

		data.append(row)

	df = pd.DataFrame(data)
	df = gen_util.combine_mag(pd.DataFrame(df))

	
	return df
	

def make_df():

	"""
	Assemble a dataframe of real and bogus metadata for machine learning
	"""
	# triplets = pd.read_csv("data_results/triplets_small_radius.csv")
	triplets = pd.read_csv("subtracted_triplets.csv")
	bogus = make_bogus_df(triplets.copy())
	triplets['bogus']=0
	bogus['bogus']=1
	
	print(len(triplets), "triplets found!")
	print(len(bogus), "boguses found!")

	triplet_df = pd.concat([triplets, bogus])
	triplet_df.to_csv("data_results/full_data.csv")

	return triplet_df


if __name__ == "__main__":

	make_df()


	

