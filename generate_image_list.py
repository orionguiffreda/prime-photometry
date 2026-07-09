import os
import re
from pathlib import Path
import pandas as pd
import time
import matplotlib.pyplot as plt
import math
import base64
import numpy as np
from datetime import datetime, timedelta
import subprocess

import glob
from ast import literal_eval
from collections import Counter

from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u


import photometrus_utils as util
import sfft_util
import assess_processed_data_util as gen_util
from assess_processed_data_util import list_dirs, list_files, contains, combine_mag, constrain_dt


def chip_search(file_list):
	"""
	identify files with matching chip label
	"""
	chips = []
	files = []
	for file in file_list:
		file = file_list[0]
		chip_search = re.search(r'C\d', file)  
		if chip_search: # if there IS a cutout labelled with a chip
			chip = chip_search.group()
			if chip not in chips:
				chips.append(chip)
				files.append(file)

	if len(chips)==1:
		return files[0], chips[0]
	if len(chips)>1:
		print("something went wrong, multiple chips in cutouts")
	return None, None
		


def fix_dates(df):
	"""
	Convert date field to datetime, create days since discovery field
	"""

	corrupted_dates = False
	if corrupted_dates: # merge with TNS csv that has correct dates
		date_csv = "~/PycharmProjects/prime-photometry-fiona/relevant_objs_2.csv"
		date_csv = pd.read_csv(date_csv)
		date_csv['full_name'] = date_csv['name_prefix']+date_csv['name']
		date_csv = date_csv[['full_name', 'discoverydate']]
		print(df.keys())
	
		df = df.drop(columns='discoverydate')
		df = df.merge(date_csv, 'left', on=['full_name'])
	
	df['discoverydate'] = pd.to_datetime(df['discoverydate'], format='%m/%d/%y')
	df['field_date'] = pd.to_datetime(df['field'].str[-10:], format='%Y-%m-%d')
	df['days_since_discovery'] = (df['field_date'] - df['discoverydate']).dt.days

	return df


def fetch_prime_data(transient_file = "~/PycharmProjects/prime-photometry-fiona/tns_objects_4_13.csv"):
	"""
	Search through processed data (starting with transient ids in TNS) and make a dataframe of all matching observations
	"""
	
	BASE = "/mnt/photometry"

	folders = 0
	
	id_tags= list_dirs(BASE)
	if not id_tags:
		print("No directories found under /mnt/photometry.")
		return

	candidates = pd.read_csv(transient_file)

	# add in discovery dates
	dates = pd.read_csv(transient_file)
	del candidates['discoverydate']
	candidates = candidates.merge(dates[['objid', 'discoverydate']], on='objid', how='left')


	print(candidates['discoverydate'])
	candidate_ids = list(candidates['full_name'])

	candidates['status']=''
	# candidates['GRB']=''
	triplets = pd.read_csv("data_results/triplets_recent.csv")
	print("num triplets:", len(triplets))
	print("num triplets small rad:", len(pd.read_csv("data_results/triplets_small_radius.csv")))
	
	triplet_names = list(triplets["full_name"])
	
	print(triplet_names)
	
	obs_rows = []
	
	for i,id_tag in enumerate(candidate_ids):
		# candidates.loc[i, 'GRB'] = gen_util.GRB_SN_Match(candidates, i) # clunky
					
		if id_tag not in id_tags: 
			continue
			
		id_path = f"{BASE}/{id_tag}"
		# print(f"\n ID Tag: {id_tag}")

		fields = list_dirs(id_path)
		for field in fields:
			field_path = f"{id_path}/{field}"
			bands = list_dirs(field_path)

			for band in bands:
				stack_path = f"{field_path}/{band}/stack" #/bogus"
				try:
					files = os.listdir(stack_path)
				except IOError:
					continue

				


				folders += 1
				# get files

				# photometry has been run such that coadds are copied into bogus directory and photometry is run 

				coadd_files = [f for f in files if f.startswith("coadd") and f.endswith('.fits')]
				coadd_files = [f for f in coadd_files if not f.endswith('aligned.fits') and not f.endswith('masked.fits') and not f.endswith('weightmap.fits')]

				coadd_chip_file, coadd_chip = chip_search(coadd_files)
				

				# enforce chip match if GRB cutout is labelled with a chip
				if coadd_chip:
					coadd_files = [f for f in coadd_files if contains(coadd_chip, f)]

								
				if len(coadd_files) == 0:
					# print('no coadd found, data not fully processed!')
					continue

				coadd_file = coadd_files[0] # select first one
				print(f"found {stack_path}/{coadd_file}")

				# get header from sky file
				fits_files = glob.glob(f"{field_path}/{band}/sky/*.fits")
				if len(fits_files) == 0:
					print("no sky fits files found!")
					continue
					
				fits_file = fits_files[0] # select the first one
				header = fits.getheader(fits_file)

				ecsv_files = [f for f in files if f.startswith('GRB') and f.endswith('.ecsv')]
				
				if False: # rerun photometry methods to determine if detection

					print("running photometrus photometry")	
					print("photometrus", "photometry",
						"-stackpath",f"{stack_path}/",
						"-band", band,
						"-chip", f"{header['CHIP']}",
						"-grb_ra", f"{candidates.loc[i, 'ra']}",
						"-grb_dec", f"{candidates.loc[i, 'declination']}",
						"-grb_radius", "1",
						"-sx_cfg", "sex2_aggr.config")
					try:
						
						subprocess.run([
							"photometrus", "photometry",
							"-stackpath",f"{stack_path}/",
							"-band", band,
							"-chip", f"{header['CHIP']}",
							"-grb_ra", f"{candidates.loc[i, 'ra']}",
							"-grb_dec", f"{candidates.loc[i, 'declination']}",
							"-grb_radius", "1",
							"-sx_cfg", "sex2_aggr.config",
						])
					except Exception as e:
						print(e)

				files = os.listdir(stack_path) # new files


				ecsv_files = [f for f in files if f.startswith('GRB') and f.endswith('.ecsv')]
				ecsv_chip_file, ecsv_chip = chip_search(ecsv_files)

				
				has_ecsv = len(ecsv_files) > 0
				ecsv = None
				if ecsv_chip:
					ecsv = ecsv_chip_file
				elif has_ecsv:
					ecsv = ecsv_files[0] # select first one
				else:
					print("non detection (no ecsv):")
					

				# bogus_dir = f"{stack_path}/bogus"
				# try:
				# 	os.listdir(bogus_dir)
				# except IOError:
				# 	s.mkdir(bogus_dir)
				# 	os.system(f'cp {stack_path}/{coadd_file} {bogus_dir}/{coadd_file}')
				# 	# files = os.listdir(bogus_dir)
				# 	print(files)

		
					
				row={
					"full_name": id_tag,
					"field": field,
					"truncated_field": field[:-11],
					"band": band,
					"chip": coadd_chip,
					"coadd_path": f"{stack_path}/{coadd_file}",
					"ecsv_path": f"{stack_path}/{ecsv}",
					"discoverydate":candidates.loc[i, 'discoverydate'],
				}
				
				header_keys = ["EXPTIMEE", "DITHTYP", "DITHRAD", "DITH_REP", "NINT", "ROTOFF"]
				
				row.update({key: header[key] for key in header_keys})

				# add in data from ecsv
				if len(ecsv_files)>0:
					df = pd.read_csv(row["ecsv_path"], sep='\\s+', comment='#')
					ecsv_data = df.to_dict(orient='records')  # list of dicts, one per row
					for ecsv_dict in ecsv_data:
						row.update(ecsv_dict)
						# print(ecsv_dict.get('autoS_CM'), type(ecsv_dict.get('autoS_CM')))

				# bogus_coords = {"AT2024acdq":[52.43301900536419, -32.098768098394586], 
				# 				"AT2025syx":[352.97641781649025, -28.65682938159262],
				# 				"AT2025wdq":[36.21625305425166, -24.817318526608844],
				# 				"AT2025foi":[193.42981039864407, -15.619197073527094],
				# 			   }

				# print(row)
				obs_rows.append(row)
				print("added row!")
				
				# year = id_tag[4:6]
	
				if False:# or year=="24" or year=="25" or year=="26": 

					files = glob.glob(f"{os.path.dirname(row['coadd_path'])}/*")
					files = [os.path.basename(f) for f in files]
					
					coadd_files = [f for f in files if f.startswith('coadd') and f.endswith('.ecsv')] 
					
					auto_files = [f for f in coadd_files if f.endswith('AUTO.ecsv')] 
					psf_files = [f for f in coadd_files if f.endswith('PSF.ecsv')] 
					aper_files = [f for f in coadd_files if f.endswith('APER.ecsv')] 
				
					updated_photometry = len(coadd_files)>0 and len(auto_files)==0 and len(psf_files)==0 and len(aper_files)==0
					
					b = gen_util.get_latitude(candidates.loc[i, 'ra'], candidates.loc[i, 'declination'], field)
					sparse = abs(b) > 16
					print("Sparse field:", sparse, "| b =",b)

					# Commands for rerunning stacking and bogus photometry
					
					# if not updated_photometry and sparse:
				
					
					# print( "running: photometrus", "process", "stack"
					# 					" -astrom_only"
					# 					" -stack", stack_path,
					# 				"-chip", f"{header['CHIP']}")
					# try:
					# 	subprocess.run([
					# 				"photometrus", "process", "stack",
					# 				"-astrom_only",
					# 				"-stack", stack_path,
					# 				"-chip", f"{header['CHIP']}"
					# 			], check=True)
					# except Exception as e:
					# 		print(f"Error: {e}")
						
						
					# subprocess.run([
					# 	"photometrus", "photometry",
					# 	"-stackpath",f"{bogus_dir}/",
					# 	"-band", band,
					# 	"-chip", f"{header['CHIP']}",
					# 	"-grb_ra", f"{bogus_coords[id_tag][0]}",
					# 	"-grb_dec", f"{bogus_coords[id_tag][1]}",
					# 	"-grb_radius", "1",
					# ])
					# except Exception as e:
					# 	print(e)


	
	obs_df = pd.DataFrame(obs_rows)
	obs_df = combine_mag(obs_df) # 	combine magnitude columns like JMag, HMag into one column
	
	print("folders:",folders)
	obs_df.to_csv("data_results/prime_observations.csv", index=False)
	candidates.to_csv("data_results/transient_candidates.csv", index=False)

	return candidates, obs_df, folders

def label_edge_transients(df):
	"""
	calculate distance (px) to edge for each observation
	"""
	df["edge"] = False
	df["px_from_edge"]=0.0
	
	for idx, row in df.iterrows():
		header = fits.getheader(str(row["coadd_path"]))
		wcs = WCS(header)
		
		coord = SkyCoord(ra=row["ra"]*u.deg, dec=row["declination"]*u.deg)
		x_pix, y_pix = wcs.world_to_pixel(coord)
		
		nx = header['NAXIS1']
		ny = header['NAXIS2']

		cd11 = header['CD1_1']
		cd12 = header['CD1_2']
		cd21 = header['CD2_1']
		cd22 = header['CD2_2']
		
		pixel_scale_x = np.sqrt(cd11**2 + cd12**2)
		pixel_scale_y = np.sqrt(cd21**2 + cd22**2)
		# should be 0.00013824 deg / px (0.49766 arcsec / px)

		edge_thresh_x = 0.05 / pixel_scale_x # 3 arcmin = 0.05 deg / deg per px
		edge_thresh_y = 0.05 / pixel_scale_y

		near_edge = (
		(x_pix < edge_thresh_x) |
		(x_pix > nx - edge_thresh_x) |
		(y_pix < edge_thresh_y) |
		(y_pix > ny - edge_thresh_y)
		)

		df.loc[idx, "edge"] = near_edge
		df.loc[idx, "px_from_edge"] = min(x_pix, y_pix, nx-x_pix, ny-y_pix)
		
	return df
	
def find_single_detections(scis):
	"""
	Find detections which do not appear redundantly in the combo commands (each unique id, band, and field combination appears only once)

	scis: pandas DataFrame of PRIME detections
	"""
	print()
	print("finding single detections...")
	
	seen = Counter()
	parent_re = re.compile(r"'parentdir':\s*'([^']+)'")
	band_re = re.compile(r"'band':\s*'([^']+)'")
	
	with open("combo_commands/combo_commands_2026-01-28T15:38:20.474500") as f:
		for line in f:
			parent = parent_re.search(line).group(1)
			field = os.path.basename(os.path.dirname(parent))[:-11] # truncate to remove date
			id_tag = os.path.basename(os.path.dirname(os.path.dirname(parent)))
			band = band_re.search(line).group(1)
			seen[(id_tag, field, band)] += 1

	scis["key"] = list(zip(scis["full_name"], scis["truncated_field"], scis["band"]))
	scis["count"] = scis["key"].map(seen).fillna(0).astype(int)
	df = scis[scis["count"] == 1].copy()
	print(len(df), "detections with (id, field, and band) which appear once combo commands from 1/28")

	df = constrain_dt(df) 

	return df

def assign_roles(group):
	"""
	assign sci and ref roles to observations based on crossmatch flags
	"""
	
	group = group.copy()
	group["role"] = "ref"

	# if any mag has identified a change
	sci_idx = -1
	for idx,item in group.iterrows():
		# if any magnitude method registered a detection
		def detection(key):
			S_CM = float(item.get(key))
			return S_CM!=1 and not math.isnan(S_CM)
			
		# many of these are all NaN	
		if detection('autoS_CM') or detection('aperS_CM') or detection('psfS_CM'):
			# if this is a new detection, or closer detection to TNS
			if abs(item['days_since_discovery']) < 20:
				if sci_idx == -1 or item['autoMag'] < group.loc[sci_idx, 'autoMag']: # brighter in this image
					sci_idx = idx

	if sci_idx != -1:
		group.loc[sci_idx, "role"] = "sci"
	

	return group
	

def make_triplet_df(candidates, obs_df, folders):
	"""
	turn df of observations into a list of triplets 
	"""
	
	obs_df = fix_dates(obs_df)

	candidates = candidates.drop(columns=['discoverydate'])
	candidate_groups = ["full_name", "band", "truncated_field"] # not chip

	# create sci and ref labels, deprecated behavior
	obs_df = (obs_df.groupby(candidate_groups, group_keys=False).apply(assign_roles))
	
	scis = obs_df.query("role=='sci'")
	refs = obs_df.query("role=='ref'")


	keys_for_ref = ["coadd_path", "field", "chip", "autoMag", "autoMag_Err", "psfMag", "psfMag_Err", "aperMag", "aperMag_Err", "days_since_discovery", "autoS_CM", "psfS_CM", "aperS_CM"]

	agg_keys = {key:list for key in keys_for_ref}
	
	refs_agg = (
		refs.groupby(candidate_groups).agg(agg_keys)
		.rename(columns=lambda col_name: f"ref_{col_name}").reset_index()
	)

	# identify triplets
	triplets = scis.merge(
		refs_agg,
		on=candidate_groups,
		how="inner",
	)
	
	# label triplets with candidate info 
	triplets = triplets.merge(
		candidates,
		how='left',
		on=["full_name"],
	)

	scis = scis.merge(
		candidates,
		how='left',
		on=["full_name"],
	)


	sci_only = find_single_detections(scis)
	sci_only = label_edge_transients(sci_only)
	obs_df = label_edge_transients(triplets)


	# save files 
	timestamp = time.time()

	sci_only.to_csv(f'data_results/to_observe_7_7.csv', index=False)
	print("Number of sci only fields that exist:", len(sci_only))


	triplets_subset = constrain_dt(triplets) 
	triplets = add_fields(triplets)
	triplets.to_csv("data_results/triplets.csv", index=False)
	
	triplet = len(triplets)
	sci = len(scis)
	ref = len(refs)
	print("distinct transients (excluding repeats with different bands, fields):",len(set(list(triplets["full_name"]))))
	print()
	print(f"processed images: triplets: {triplet}, only sci: {sci-triplet}, only ref: {ref-triplet}, partially processed: {folders-sci-ref+triplet}")

	return candidates, triplets, sci_only


def add_fields(df, basic_cuts = False):
	"""
	Add analysis fields assessing data quality: difference between magnitude methods, difference between epochs, galactic coordinates 
	"""

	# Other potential constraints:
		  # SNR for 8 px diameter aperture >5
		  # 0 < Flux (8 px) / Flux (18 px) < 1.5
		  # Number of negative pixels in a 5 × 5 pixel area <= 13;
		  # Number of bad pixels in a 5 × 5 pixel area <= 7; 

	df["mag_dif"]= None
	df["sci_ref_dif"] = None
	df["average_sci_mag"] = None
	df["b"]= None
	
	for i, row in df.iterrows():
		df.loc[i,"mag_dif"] = abs(float(row['psfMag']) - float(row['autoMag']))
		mag_keys = ["autoMag", "aperMag"] # not psf?
		sci_ref_difs = [float(row[key]) - np.nanmean(row[f'ref_{key}']) for key in mag_keys]
		df.loc[i,"sci_ref_dif"] = np.nanmean(sci_ref_difs)
		df.loc[i,"average_sci_mag"] = np.nanmean([float(row[key]) for key in mag_keys])
		df.loc[i,"b"] = gen_util.get_latitude(row['ra'], row['declination'], row['field'])

		diff_path = glob.glob(f"/mnt/photometry/{row['full_name']}/*.fits")	

		diff_path = [f for f in diff_path if os.path.basename(f).startswith("sub_4_16")] # or f.endswith(".sfftdiff.fits")]
		cat = ""
		for path in diff_path:
			cat_path = glob.glob(f"{path[:-4]}*.cat")
			cat_path.extend(glob.glob(f"{path}*.cat"))
			if len(cat_path)>0:
				FITS_DIFF = path
				cat = cat_path[0]


	if basic_cuts: 
		df = df[df["px_from_edge"]>300]
		df = df[df["SNR"]>5] 
		df = df[df["Elongation"] <= 2]
		df = df[abs(df["b"]) > 16]
		df = df[abs(df['days_since_discovery'])<20]

	return df



	
if __name__ == "__main__":
	# candidates, obs_df, folders = fetch_prime_data()

	candidates = pd.read_csv("data_results/transient_candidates.csv")
	obs_df = pd.read_csv("data_results/prime_observations.csv")
	folders=5610

	
	# good_sci = ["SN2025wzh", "AT2025umy", "AT2025wdt", "SN2025abgi", "AT2025abcc", "AT2025aahv", "AT2025zwq", "SN2025yrt", "AT2025zgx","AT2025yzj", "SN2025vzq", "AT2025vma", "AT2025vdr", "AT2025uwp", "AT2025ukj",  "AT2025szq", "AT2025swv", "AT2025rkt", "AT2025smc", "AT2025rxk", "AT2025rnf", "AT2025qsr", "AT2025nqj", "AT2025gee", "SN2025hkm", "AT2025ika", "AT2025ihi", "AT2025dil", "SN2024qje", "AT2024xwo", "AT2024vuk", "AT2024vpy", "AT2024vlx", "AT2024uwv", "AT2024uwy", "AT2024uwk", "SN2024sfg", "AT2024sjl", "AT2024scc", "AT2024rwy", "AT2024qir", "AT2024dtu", "AT2024eel", "AT2024bcj", "AT2024axt", "AT2024ase"]
	# to_obs = obs_df[obs_df['full_name'].isin(good_sci)]
	# to_obs.to_csv("to_observe_4_28.csv")	
	
	# sci_df = pd.read_csv("data_results/to_observe_4_13.csv")
	
	# timestamp = time.time()


	candidates, triplet_df, sci_df = make_triplet_df(candidates, obs_df, folders)
	# print(sci_df)
	
	# gen_util.save_summary_statistics(candidates, triplet_df)
	# generate_plots(triplet_df)


