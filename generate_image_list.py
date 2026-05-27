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


import photometrus.photometry.photometry as photo
import photometrus_utils as util
import sfft_util
import assess_processed_data_util as gen_util
# import SFFT

# from photometrus.archive_reducer import get_all_date_logs


# ---------------- CONFIG ----------------
BASE = "/mnt/photometry"
# ----------------------------------------

def list_dirs(path):
	"""Return list of directory names under a path."""
	try:
		return [f for f in os.listdir(path) if os.path.isdir(os.path.join(path,f))]
	except FileNotFoundError:
		return []

def list_files(path):
	"""Return list of file names under a path."""
	try:
		return [f for f in os.listdir(path) if not os.path.isdir(os.path.join(path,f))]
	except FileNotFoundError:
		return []

def contains(mystr, file):
	return re.search(mystr, file) is not None


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
	# 'full_tns_objects_11_15.csv?
	# date_csv = "~/PycharmProjects/prime-photometry-fiona/relevant_objs_2.csv"
	# date_csv = pd.read_csv(date_csv)
	# date_csv['full_name'] = date_csv['name_prefix']+date_csv['name']
	# date_csv = date_csv[['full_name', 'discoverydate']]
	# print(df.keys())

	# df = df.drop(columns='discoverydate') # in werid format i dont understand
	# df = df.merge(date_csv, 'left', on=['full_name'])
	
	# df['discoverydate'] = pd.to_datetime(df['discoverydate'],format='%Y-%m-%d %H:%M:%S.%f')
	df['discoverydate'] = pd.to_datetime(df['discoverydate'], format='%m/%d/%y')
	df['field_date'] = pd.to_datetime(df['field'].str[-10:], format='%Y-%m-%d')
	df['days_since_discovery'] = (df['field_date'] - df['discoverydate']).dt.days

	return df


def fetch_prime_data(transient_file = "~/PycharmProjects/prime-photometry-fiona/tns_objects_4_13.csv"):
	"""
	search through processed data (starting with transient ids in TNS) and make a dataframe of all matching observations
	"""
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
					
		if id_tag not in ['AT2024acdq', 'AT2025syx', 'AT2025wdq', 'AT2025foi']:#id_tags: 
			continue
			
		id_path = f"{BASE}/{id_tag}"
		print(f"\n ID Tag: {id_tag}")

		fields = list_dirs(id_path)
		for field in fields:
			field_path = f"{id_path}/{field}"
			bands = list_dirs(field_path)

			for band in bands:
				stack_path = f"{field_path}/{band}/stack/bogus"
				try:
					files = os.listdir(stack_path)
				except IOError:
					continue

				try:
						bogus_coords = {"AT2024acdq":[52.43301900536419, -32.098768098394586], 
										"AT2025syx":[352.97641781649025, -28.65682938159262],
										"AT2025wdq":[36.21625305425166, -24.817318526608844],
										"AT2025foi":[93.42981039864407, -15.619197073527094],
									   }
						
						print("running photometrus photometry")
						# bogus_directory = stack_path+'/bogus'
						# if not os.path.exists(bogus_directory):
						# 	os.mkdir(bogus_directory)
						# os.system(f'cp {stack_path}/{coadd_file} {bogus_directory}/{coadd_file}')

						subprocess.run([
							"photometrus", "photometry",
							"-stackpath",f"{stack_path}/",
							"-band", band,
							"-chip", f"{header['CHIP']}",
							"-grb_ra", f"{bogus_coords[id_tag][0]}",
							"-grb_dec", f"{bogus_coords[id_tag][1]}",
							"-grb_radius", "2"
						])
				except Exception as e:
					print(e)

				
				folders += 1
				# get files
				# cutout_files = [f for f in files if contains('Cutout', f) and f.endswith(".fits")]

				ecsv_files = [f for f in files if f.startswith('GRB') and f.endswith('.ecsv')]
				ecsv_chip_file, ecsv_chip = chip_search(ecsv_files)

				coadd_files = [f for f in files if f.startswith("coadd") and f.endswith('.fits')]
				coadd_chip_file, coadd_chip = chip_search(coadd_files)

				
				if len(coadd_files) ==0:
					print('no coadd found, data not fully processed!')
					continue
				
				# if len(cutout_files) == 0:
				# 	print('no cutout found, data not fully processed!')
				# 	continue

				# enforce chip match if GRB cutout is labelled with a chip
				if coadd_chip:
					coadd_files = [f for f in coadd_files if contains(coadd_chip, f)]
				
				coadd_file = coadd_files[0] # select first one
				print(f"found {stack_path}/{coadd_file}")

				has_ecsv = len(ecsv_files) > 0
				ecsv = None
				if ecsv_chip:
					ecsv = ecsv_chip_file
				elif has_ecsv:
					ecsv = ecsv_files[0] # select first one
					
				row={
					"full_name": id_tag,
					"field": field,
					"truncated_field": field[:-11],
					"band": band,
					"chip": coadd_chip,
					"coadd_path": f"{stack_path}/{coadd_file}",
					# "has_ecsv": len(ecsv_files) > 0,
					"ecsv_path": f"{stack_path}/{ecsv}",
					"discoverydate":candidates.loc[i, 'discoverydate'],
					# "n_cat": count_cat_lines(f"{field_path}/{band}", candidates, i),
				}
				

				if not has_ecsv:
					print("non detection (no ecsv):",row)

				fits_files = glob.glob(f"{field_path}/{band}/sky/*.fits")
				if len(fits_files) == 0:
					print("no sky fits files found!")
					continue
					
				fits_file = fits_files[0] # select the first one
				header = fits.getheader(fits_file)
				header_keys = ["EXPTIMEE", "DITHTYP", "DITHRAD", "DITH_REP", "NINT", "ROTOFF"]
				
				row.update({key: header[key] for key in header_keys})

				# add in data from ecsv
				if len(ecsv_files)>0:
					df = pd.read_csv(row["ecsv_path"], sep='\s+', comment='#')
					ecsv_data = df.to_dict(orient='records')  # list of dicts, one per row
					for ecsv_dict in ecsv_data:
						row.update(ecsv_dict)
						# print(ecsv_dict.get('autoS_CM'), type(ecsv_dict.get('autoS_CM')))

										
				obs_rows.append(row)
				print("added row!")
				year = id_tag[4:6]
				# print(year)

				
				if year=="24" or year=="25" or year=="26": 

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
					
					# reg_files = glob.glob(
					

					# print( "running: photometrus", "process", "stack"
					# 				" -astrom_only"
					# 				" -stack", stack_path,
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
					
					# if not updated_photometry and sparse:
						
					try:
						bogus_coords = {"AT2024acdq":[52.43301900536419, -32.098768098394586], 
										"AT2025syx":[352.97641781649025, -28.65682938159262],
										"AT2025wdq":[36.21625305425166, -24.817318526608844],
										"AT2025foi":[93.42981039864407, -15.619197073527094],
									   }
						
						# print("running photometrus photometry")
						# bogus_directory = stack_path+'/bogus'
						# if not os.path.exists(bogus_directory):
						# 	os.mkdir(bogus_directory)
						# os.system(f'cp {stack_path}/{coadd_file} {bogus_directory}/{coadd_file}')

						# subprocess.run([
						# 	"photometrus", "photometry",
						# 	"-stackpath",f"{bogus_directory}/",
						# 	"-band", band,
						# 	"-chip", f"{header['CHIP']}",
						# 	"-grb_ra", f"{bogus_coords[id_tag][0]}",
						# 	"-grb_dec", f"{bogus_coords[id_tag][1]}",
						# 	"-grb_radius", "2"
						# ])

						# subprocess.run([
						# 	"photometrus", "photometry",
						# 	"-stackpath",{stack_path}/{coadd_file},
						# 	"-band", band,
						# 	"-chip", f"{header['CHIP']}",
						# 	"-grb_ra", f"{candidates.loc[i, 'ra']}",
						# 	"-grb_dec", f"{candidates.loc[i, 'declination']}",
						# 	"-grb_radius", "2"
						# ])

					except Exception as e:
		
						print(f"Error: {e}")
						# continue

	# sftp.close()
	# ssh.close()

	# print(obs_rows)
	
	obs_df = pd.DataFrame(obs_rows)
	obs_df = combine_mag(obs_df)
	#obs_df = fix_dates(obs_df)

	
	print("folders:",folders)
	obs_df.to_csv("data_results/prime_observations_bogus.csv", index=False)
	candidates.to_csv("data_results/transient_candidates_bogus.csv", index=False)

	return candidates, obs_df, folders

def label_edge_transients(df):
	"""
	calculate distance (px) to edge for each observation
	"""
	df["edge"] = False
	df["px_from_edge"]=0
	
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
			# if id_tag == 'AT2025xfk':
			# 	print('saw AT2025xfk:', (id_tag, field, band))

	scis["key"] = list(zip(scis["full_name"], scis["truncated_field"], scis["band"]))
	scis["count"] = scis["key"].map(seen).fillna(0).astype(int)
	print(scis[scis['full_name'] == 'AT2025xfk' ]["count"])
	df = scis[scis["count"] == 1]
	print(len(df), "detections with (id, field, and band) which appear once combo commands from 1/28")
	# print("AT2025xfk in scis:", "AT2025xfk" in list(df['full_name']))

	# enforce date constraints
	start_date = datetime(2024, 1, 1, 0, 0, 0)
	df = df[pd.to_datetime(df['discoverydate']) > start_date]
	print(len(df), "detections which after 2023")

	df = df[df['days_since_discovery'] < 20]
	print(len(df), "detections within 20 day window of ZTF measurement")

	return df

def assign_roles(group):
	group = group.copy()
	group["role"] = "ref"

	# if any mag has identified a change
	sci_idx = -1
	# print(group.keys())
	for idx,item in group.iterrows():
		# if any magnitude method registered a detection
		# print(item)
		def detection(key):
			S_CM = float(item.get(key))
			return S_CM!=1 and not math.isnan(S_CM)
			
		# many of these are all NaN	
		if detection('autoS_CM') or detection('aperS_CM') or detection('psfS_CM'):
			# if this is a new detection, or closer detection to TNS
			if abs(item['days_since_discovery']) < 20:
				if sci_idx == -1 or item['autoMag'] < group.loc[sci_idx, 'autoMag']: # brighter in this image
					# print("real detection")
					sci_idx = idx

	if sci_idx != -1:
		group.loc[sci_idx, "role"] = "sci"
	# print(group[["full_name", "aperMag","days_since_discovery","chip","role"]])
	# print()

	return group



def make_triplet_df(candidates, obs_df, folders):
	obs_df = fix_dates(obs_df)
	print(obs_df.head(10))

	# print("das:",obs_df[obs_df["full_name"]=="AT2025das"][['field', 'autoMag', 'autoMag_Err', 'psfMag', 'psfMag_Err', 'ra', 'declination']])

	candidates = candidates.drop(columns=['discoverydate'])
	candidate_groups = ["full_name", "band", "truncated_field"] # not chip
	# print(obs_df.head())
	# create sci and ref labels
	obs_df = (
		obs_df
		.groupby(candidate_groups, group_keys=False)
		.apply(assign_roles)
	)
	print("keys!", obs_df.keys())

	scis = obs_df.query("role=='sci'")
	refs = obs_df.query("role=='ref'")

	# group references that go with the same triplet

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
	print("initial merge:",triplets.head())
	# label triplets with candidate info 
	triplets = triplets.merge(
		candidates,
		how='left',
		on=["full_name"],
	)
	print("second merge:",triplets.head())

	scis = scis.merge(
		candidates,
		how='left',
		on=["full_name"],
	)

	print(scis.keys())

	
	sci_only = find_single_detections(scis)
	sci_only = label_edge_transients(sci_only)
	# sci_only_minus_edge = sci_only_minus_edge[sci_only_minus_edge["px_from_edge"] > 300]
	obs_df = label_edge_transients(triplets)


	# save files 
	timestamp = time.time()

	print("saving obs df")


	# obs_df.to_csv("data_results/obs_with_groups.csv")
	sci_only.to_csv(f'data_results/to_observe_4_13.csv', index=False)
	# print("sci only:", len(sci_only), list(sci_only["full_name"]))


	def constrain_dt(df):
		"""
		enforce: discovery date is after 2023 and within 20 days of sci
		"""
		df['discoverydate'] = pd.to_datetime(df['discoverydate'])
		print(len(df), "total triplets")
		start_date = datetime(2024, 1, 1, 0, 0, 0)
		df = df[df['discoverydate'] > start_date]
		print(len(df), "detections after 2023")
		df['days_since_discovery_abs'] = abs(df['days_since_discovery'])
		df = df[df['days_since_discovery_abs'] < 8]
		print(len(df), "detections within 7 day window of ZTF measurement")
		return df

	# triplets_subset = constrain_dt(triplets) #UNCOMMENT
	triplets.to_csv("data_results/triplets_bogus.csv", index=False)
	
	triplet = len(triplets)
	sci = len(scis)
	ref = len(refs)
	print("distinct transients:",len(set(list(triplets["full_name"]))))
	print(len(set(list(triplets["full_name"]))))
	print(f"triplets: {triplet}, only sci: {sci-triplet}, only ref: {ref-triplet}, partially processed: {folders-sci-ref+triplet}")
	print("\n Complete.")

	return candidates, triplets, sci_only



def get_relevant_info(id_tag, df):


	row = df[df["full_name"]==id_tag]
	
	def get_item(key):
		return row[key].iloc[0]
		
	ra, dec = get_item('ra'), get_item('declination')
	
	
	FITS_SCI = str(row["coadd_path"].iloc[0])
	FITS_REF = str(row["ref_coadd_path"].iloc[0][0])

	sci_stack = os.path.dirname(FITS_SCI)
	ref_stack = os.path.dirname(FITS_REF)

	sci_cat, sci_cat_png, fwhm_sci = sfft_util.get_cat(FITS_SCI)
	ref_cat, ref_cat_png, fwhm_ref = sfft_util.get_cat(FITS_REF)


	
	do_subtractions = False
		
	diff_path = glob.glob(f"/mnt/photometry/{id_tag}/*.fits")	


	diff_path = [f for f in diff_path if os.path.basename(f).startswith("sub_4_16")] # or f.endswith(".sfftdiff.fits")]
	cat = ""
	for path in diff_path:
		cat_path = glob.glob(f"{path[:-4]}*.cat")
		cat_path.extend(glob.glob(f"{path}*.cat"))
		if len(cat_path)>0:
			FITS_DIFF = path
			cat = cat_path[0]

	
	if id_tag=="AT2025xfk":# getting stuck
		return ""
		
	if cat=="" or do_subtractions:

		print('SKIPPING UNSUBTRACTED FIELD')
		return ""
		
		print("doing subtractions, extracting sources...")
		FITS_SCI, FITS_REF, FITS_DIFF, PixA_DIFF, SFFTPrepDict, sci_bkg_path, ref_bkg_path =  sfft_util.SFFT(FITS_SCI, FITS_REF)
		rms_ref, rms_sci, rms_diff_predicted, rms_diff = sfft_util.get_diff_stats(SFFTPrepDict, PixA_DIFF)
		udx, sdx, udy, sdy = sfft_util.sfft_source_reg_gen(SFFTPrepDict)
		sources = sfft_util.source_extract(FITS_DIFF)

	else: 
		print()

		print("reading in diff image and diff source catalog")
		print("catalog:", cat)
		
		# sci_bkg, sci_bkg_path = sfft_util.get_bkg(FITS_SCI)
		# ref_bkg, ref_bkg_path = sfft_util.get_bkg(FITS_REF)

		# sci_bkg_path, ref_bkg_path = glob.glob(f"{sci_stack}/background.png")[0], glob.glob(f"{ref_stack}/background.png")[0]

		# sci_bkg, sci_bkg_path = sfft_util.get_bkg(FITS_SCI)
		# ref_bkg, ref_bkg_path = sfft_util.get_bkg(FITS_REF)
		PixA_DIFF = fits.getdata(FITS_DIFF)

		rms_sci, rms_ref, rms_diff = sfft_util.rms(fits.getdata(FITS_SCI), "sci"), sfft_util.rms(fits.getdata(FITS_REF), "ref"), sfft_util.rms(PixA_DIFF, "diff")
		sources = Table.read(cat, hdu=2)


	print("sci, ref, dif:",FITS_SCI, FITS_REF, FITS_DIFF)
	diff_dir = os.path.dirname(FITS_DIFF)


	html_content = ""

	def prnt(string_content):
		return "<p>"+string_content+"</p>"
		
	
	html_content += prnt(f"\nSTATS for {id_tag} (discovery mag: {get_item('discoverymag')}, sci ref difference: {round(get_item('sci_ref_dif'),3)}, {get_item('band')} band, ra={round(ra,4)}, dec={round(dec,4)}, chip={get_item('chip')})")



	def print_psf(header, img_str):

		html_content = ""

		try:
			html_content += prnt(f"{img_str} (psf, auto) lim mag: {float(header['HIERARCH lim_mag_psf'])}, {float(header['HIERARCH lim_mag_auto'])}")
		except KeyError as e:
			try:
				html_content += prnt(f"{img_str} (psf, auto) lim mag: {float(header['LM_PSF'])}, {float(header['LM_AUTO'])}")
			except KeyError as e2:
				html_content += prnt(f"Error:{e}, {e2}")
		return html_content


	def print_mags(img_str, image):
		header = fits.getheader(image)
		html_content = ""
		def get_item(key):
			if img_str == "sci":
				return row[key].iloc[0] # get instead of idx?
			else:
				return row[f"ref_{key}"].iloc[0]
			
		try:

			CM_list = [get_item('autoS_CM'), get_item('aperS_CM'), get_item('psfS_CM')]
			CM_arr = np.array(CM_list).flatten()
			# print(CM_arr)
			color = "yellow" if 2 in CM_arr else "green"
			color = "cyan" if 3 in CM_arr else color
			
			 #ff0000
			open_span, close_span = f"<span style=\"background-color:{color};\">", "</span>"
			html_content += prnt(f"{open_span}Survey Crossmatch (auto, aper, psf): {get_item('autoS_CM')}, {get_item('aperS_CM')}, {get_item('psfS_CM')}{close_span}")
			
			mag_arr = np.array([get_item('autoMag'),get_item('aperMag'),get_item('psfMag')]).flatten()
			avg_mag = round(np.average(mag_arr),3)
			stdev = round(np.std(mag_arr),3)
			html_content += prnt(f"<b>{img_str} mag: &mu;={avg_mag}, &sigma;={stdev}</b> \n(auto: {get_item('autoMag')}+-{get_item('autoMag_Err')}, aper:  {get_item('aperMag')} +- {get_item('aperMag_Err')}, psf: {get_item('psfMag')} +- {get_item('psfMag_Err')})")

			html_content += print_psf(header, img_str)
			html_content += prnt(f"{img_str}: {get_item('field')} | {get_item('days_since_discovery')} days since discovery")
			

		except KeyError as e:
			html_content += prnt(e)

		html_content += prnt("")
		return html_content


	
	def diff(sci_stack, ref_stack, FITS_DIFF, PixA_DIFF, ra, dec):
		html_content = ""

		flipped_name = os.path.dirname(FITS_DIFF) + '/' + 'diff_flipped.fits'
		neg_residuals = -1 * PixA_DIFF
		fits.writeto(flipped_name, neg_residuals, fits.getheader(FITS_DIFF), overwrite=True)

		mag = util.forced_photometry(ra, dec, FITS_DIFF)
		flipped_mag = util.forced_photometry(ra, dec, flipped_name)
		name_ext = "diff"
		if math.isnan(mag) or not math.isnan(flipped_mag) and flipped_mag < mag:
			print("FLIPPING MAG")
			PixA_DIFF *= -1
			FITS_DIFF = flipped_name
			name_ext += "flipped"
		
		good_ztf_sources = sfft_util.ztf_cuts(sources, PixA_DIFF, crop=500)
		source_savename, source_threshname, closest_ra, closest_dec = sfft_util.closest_source(good_ztf_sources, FITS_DIFF, ra, dec)
		distance = util.distance(ra, dec, closest_ra, closest_dec)
		diff_savename, diff_threshname = util.make_cutout(FITS_DIFF, ra, dec, png=True, display_file=False, name_ext=name_ext)

		
			
		html_content += f"<p>AB mag: {mag}</p><p>AB mag flipped: {flipped_mag}</p>"

		html_content += f"<p>distance to closest source: {distance} arcsec (ra={round(closest_ra,4)}, dec={round(closest_dec,4)})</p>"
		html_content += f"<p>RMS: sci={round(rms_sci,3)}, ref={round(rms_ref,3)}, diff={round(rms_diff,3)}</p>" # pred={round(rms_diff_predicted,3)}
		# html_content += f"<p>&mu;x={udx}, &sigma;x={sdx}, &mu;x={udy}, &sigma;x={sdy}</p>"
		html_content += f"{len(good_ztf_sources)} good / {len(sources)} raw sources in diff image"
		html_content += f"<p>pixels from edge: {round(get_item('px_from_edge'), 3)}</p>"



		return html_content, diff_savename #(diff_dir, source_savename)

	def source_info():
		html='distance to source: '
		for i,item in enumerate(row['RA']):
			prime_ra, prime_dec = row['RA'].iloc[i], row['DEC'].iloc[i]
			html+= f"{util.distance(ra, dec, prime_ra, prime_dec)},"
		return html

	diff_html_content, diff_savename = diff(sci_stack, ref_stack, FITS_DIFF, PixA_DIFF, ra, dec)

	
	sci_savename, sci_threshname = util.make_cutout(FITS_SCI, ra, dec, png=True, display_file=False, name_ext="sci")
	ref_savename, ref_threshname = util.make_cutout(FITS_REF, ra, dec, png=True, display_file=False, name_ext="ref")
	# diff_savename, diff_threshname = util.make_cutout(FITS_DIFF, ra, dec, png=True, display_file=False, name_ext="diff")
	diff_bigger_savename, diff_bigger_threshname = util.make_cutout(FITS_DIFF, ra, dec, png=True, display_file=False, name_ext="diff_bigger", photoDistThresh=30.0)
	# diff_png = util.display(FITS_DIFF, png=True, show=False, savename="diff_whole")




	def append_table_row(contents):
		html_content = "<tr>"
		for item in contents:
			html_content +=f"<td>{item}</td>"
		html_content += "<tr>"
		return html_content

	def seeing(stack):
		return round(sfft_util.seeing(stack), 3)


	style = """
	<style>
		table, th, td {
		  border: 1px solid black; /* Sets border width, style, and color for all elements */
		  border-collapse: collapse; /* Merges adjacent borders into a single line */
		}
		th, td {
		  padding: 10px; /* Adds space between content and borders */
		  text-align: left;
		}
		</style>
	"""

	# could ADD source extractor, will make slow

	html_content+=f"""
	<table>
	{style}
	  <thead>
		{append_table_row(["SCI", "REF", "DIFF"])}
	  </thead>
	  <tbody>
		{append_table_row([gen_util.GRB_png_html(sci_stack, sci_savename)+source_info(), gen_util.GRB_png_html(ref_stack, ref_savename), gen_util.GRB_png_html(diff_dir,diff_savename)])}
		
		{append_table_row([gen_util.GRB_png_html_full(sci_cat_png)+f"fwhm:{fwhm_sci}", gen_util.GRB_png_html_full(ref_cat_png)+f"fwhm:{fwhm_ref}", gen_util.GRB_png_html(diff_dir, diff_bigger_savename)])}
		
		{append_table_row([print_mags("sci",FITS_SCI) , print_mags("ref",FITS_REF), diff_html_content])}
	  </tbody>
	</table>
	"""	

	return html_content
		

	
def get_relevant_info_sci_only(id_tag, df):

	row = df[df["full_name"]==id_tag]
	
	def get_item(key):
		return row[key].iloc[0]
		
	ra, dec = get_item('ra'), get_item('declination')
	
	FITS_SCI = str(row["coadd_path"].iloc[0])
	sci_stack = os.path.dirname(FITS_SCI)
	sci_bkg, sci_bkg_path = sfft_util.get_bkg(FITS_SCI)

	rms_sci = sfft_util.rms(fits.getdata(FITS_SCI), "sci")


	html_content = ""

	def prnt(string_content):
		return "<p>"+string_content+"</p>"
		
	
	html_content += prnt(f"\nSTATS for {id_tag} (discovery mag: {get_item('discoverymag')}, {get_item('band')} band, ra={round(ra,4)}, dec={round(dec,4)}, chip={get_item('chip')})")

	def print_psf(header, img_str):
		html_content = ""
		try:
			html_content += prnt(f"{img_str} psf lim mag: {float(header['HIERARCH lim_mag_psf'])} | {img_str} auto lim mag: {float(header['HIERARCH lim_mag_auto'])}")
		except KeyError as e:
			try:
				html_content += prnt(f"{img_str} psf lim mag: {float(header['LM_PSF'])} | {img_str} auto lim mag: {float(header['LM_AUTO'])}")
			except KeyError as e2:
				html_content += prnt(f"Error:{e}, {e2}")
		return html_content


	def print_mags(img_str, image):
		header = fits.getheader(image)
		html_content = ""
		def get_item(key):
			if img_str == "sci":
				return row[key].iloc[0] # get instead of idx?
			else:
				return row[f"ref_{key}"].iloc[0]
			
		try:

			CM_list = [get_item('autoS_CM'), get_item('aperS_CM'), get_item('psfS_CM')]
			CM_arr = np.array(CM_list).flatten()
			# print(CM_arr)
			color = "yellow" if 2 in CM_arr else "green"
			color = "cyan" if 3 in CM_arr else color
			
			 #ff0000
			open_span, close_span = f"<span style=\"background-color:{color};\">", "</span>"
			html_content += prnt(f"{open_span}Survey Crossmatch (auto, aper, psf): {get_item('autoS_CM')}, {get_item('aperS_CM')}, {get_item('psfS_CM')}{close_span}")
			
			mag_arr = np.array([get_item('autoMag'),get_item('aperMag'),get_item('psfMag')]).flatten()
			avg_mag = round(np.average(mag_arr),3)
			stdev = round(np.std(mag_arr),3)
			html_content += prnt(f"<b>{img_str} mag: &mu;={avg_mag}, &sigma;={stdev}</b> \n(auto: {get_item('autoMag')}+-{get_item('autoMag_Err')}, aper:  {get_item('aperMag')} +- {get_item('aperMag_Err')}, psf: {get_item('psfMag')} +- {get_item('psfMag_Err')})")

			html_content += print_psf(header, img_str)
			html_content += prnt(f"{img_str}: {get_item('field')} | {get_item('days_since_discovery')} days since discovery")
			

		except KeyError as e:
			html_content += prnt(e)

		html_content += prnt("")
		return html_content


	def source_info():
		html=''
		html+= f"distance to source:"
		for i,item in enumerate(row['RA']):
			prime_ra, prime_dec = row['RA'].iloc[i], row['DEC'].iloc[i]
			html+= f"{util.distance(ra, dec, prime_ra, prime_dec)}|"
		return html

	def get_num_sources(image):
		sources = sfft_util.source_extract(image)
		
		return f"<p>{len(sources)} sources total</p>"


	sci_savename, sci_threshname = util.make_cutout(FITS_SCI, ra, dec, png=True, display_file=False, name_ext="sci")
	# sci_bigger_savename, sci_bigger_threshname = util.make_cutout(FITS_SCI, ra, dec, png=True, display_file=False, name_ext="sci_bigger", photoDistThresh=30.0)

	def append_table_row(contents):
		html_content = "<tr>"
		for item in contents:
			html_content +=f"<td>{item}</td>"
		html_content += "<tr>"
		return html_content


	style = """
	<style>
		table, th, td {
		  border: 1px solid black; /* Sets border width, style, and color for all elements */
		  border-collapse: collapse; /* Merges adjacent borders into a single line */
		}
		th, td {
		  padding: 10px; /* Adds space between content and borders */
		  text-align: left;
		}
		</style>
	"""

	html_content+=f"""
	<table>
	{style}
	  <thead>
		{append_table_row(["Cutout", "Bkg", "Stats"])}
	  </thead>
	  <tbody>
		{append_table_row([gen_util.GRB_png_html(sci_stack, sci_savename)+source_info()+ get_num_sources(FITS_SCI), gen_util.GRB_png_html_full(sci_bkg_path), print_mags("sci",FITS_SCI)])}
	  </tbody>
	</table>
	"""	

	return html_content
		
		

def combine_mag(triplet_df):
	"""
	combine magnitude columns like JMag, HMag into one column
	"""
	mag_cols = [c for c in triplet_df.columns if c[0] in ['J','H','Y','Z']]

	for base in {c[1:] for c in mag_cols}: # "Mag" instead of "JMag"
		cols = [c for c in mag_cols if c[1:] == base] # ["JMag", "HMag", ...]
		if len(cols) == 1:
			triplet_df[base] = triplet_df[cols[0]]
		else:
			if (triplet_df[cols].notnull().sum(axis=1) > 1).any(): # If multiple cols are not null
				raise ValueError(f"Multiple variants for '{base}' in same row!\n{triplet_df[cols]}")
			triplet_df[base] = triplet_df[cols].bfill(axis=1).iloc[:,0] # Fill columns across the board and then just copy the first one into the new column
		triplet_df = triplet_df.drop(columns=cols)

	return triplet_df

def apply_ZTF_cuts(df):

	# these are for the difference, but they should certainly be true of the science image then

	# missing:
		  # SNR for 8 px diameter aperture >5
		  # 0 < Flux (8 px) / Flux (18 px) < 1.5
		  # Number of negative pixels in a 5 × 5 pixel area <= 13;
		  # Number of bad pixels in a 5 × 5 pixel area <= 7; 
	# # df = df[df["px_from_edge"]>10]

	# df = df[df["SNR"]>5] #UNCOMMENT
	# df = df[df["Elongation"] <= 2] #UNCOMMENT

	
	

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


	# constraints
	df = df[abs(df["b"]) > 16]
	df = df[abs(df['days_since_discovery'])<20]

	return df

def get_html_header():
	html_content = """
		<!DOCTYPE html>
		<html>
		<head>
			<link rel="stylesheet" href="styles.css">
		</head>
		<style>
		body {background-color: powderblue;}
		h1   {color: black;}
		p	{color: black;}
		
		#wrapper {width: 1100px; border: 1px solid black, margin;15px}
		#first {width: 300px; float:left; padding:15px;}
		#second {width: 300px; float:left; padding:15px;}
		#third {padding:15px;}
		</style>
		"""
	
	return html_content
	
	
def make_html(triplet_df):

	# triplet_df = combine_mag(triplet_df)
	# combine mag into one col?
	print(triplet_df.keys())
	print("num triplets with small dt",len(triplet_df))

	best_candidates = apply_ZTF_cuts(triplet_df)
	
	best_candidates = best_candidates.sort_values('sci_ref_dif', ascending=True)

	detections=["SN2025tpu", "SN2024aeqs", "SN2024xyu", "AT2025foi", "AT2024ahg", "AT2024acdq", "AT2025gex", "AT2025syx", "AT2025zee", "SN2025ihl","AT2025sox", "AT2025ably", "AT2025rpv"]
	
	# neg_residuals=["AT2024aevn", "AT2025bnu", "AT2024dfx", "AT2025pwa", "AT2025bnv", "AT2025fdc", "AT2024svc", "AT2024zoz", "AT2025vjl", "AT2024fje", "AT2024xwh", "AT2025wog", "AT2024yfr", "AT2024yfg"]
	# really_bad = []
	# processed = ["AT2024aevn"]#["AT2025admo", "AT2025bob", "AT2024ahec", "AT2024aevn", "AT2025aayv", "SN2024aeqs", "SN2025tpu", "AT2025btg", "AT2025bnx"]

	# best_candidates = best_candidates[best_candidates["full_name"].isin(processed)]
	# detections = []

	
	# noisy = ["AT2025admo", "AT2025adoj", "AT2024dfx", "AT2025ssh", "AT2025zx", "AT2024ahtd", "AT2024ahxv", "AT2025abgu"] # edge, poor image, ..

	# best_candidates = best_candidates[~best_candidates["full_name"].isin(dense)]
	# best_candidates = best_candidates[~best_candidates["full_name"].isin(noisy)]

	# best_candidates = best_candidates[best_candidates["full_name"].isin(["AT2025gex", "AT2025ssh"])]
				  
	print(best_candidates[['full_name', 'days_since_discovery', 'sci_ref_dif']])
	# best_candidates = pd.concat([best_candidates,best_candidates_picked])
	# print(best_candidates[['full_name', 'SNR', 'Elongation', 'Mag', 'apMag2', 'autoMag', 'psfMag', 'coadd_path', 'ref_coadd_path']].head(20))

	html_content = get_html_header() 
	html_content += """
	<body>
		<h1>Sci-Ref pairs currently processed with updated photometry methods</h1>
	"""
#500px width before
	
	for name in best_candidates['full_name']:
		try:
		    html_content += get_relevant_info(name, triplet_df)
		except Exception as e:
		    print(f"Error: {e}")
			

	html_content+="""
		</body>
		</html>
		"""

	filename= "5_15_test"
	
	with open(f'/home/alex/PycharmProjects/prime-photometry-fiona/triplet_{filename}.html', 'w') as f:
		f.write(html_content)
		print(f"HTML file '{filename}' generated successfully.")


# def make_html_only_sci(sci_df):

# 	sci_df['days_since_discovery_abs'] = abs(sci_df['days_since_discovery'])
# 	best_candidates = sci_df[sci_df['days_since_discovery_abs']<20]
# 	print(best_candidates[['full_name', 'days_since_discovery']])
	
# 	html_content = get_html_header() 
# 	html_content += """
# 	<body>
# 		<h1>Science only candidates currently processed with updated photometry methods</h1>
# 	"""
# #500px width before
	
# 	for name in best_candidates['full_name']:
# 		# try:
# 		html_content += get_relevant_info_sci_only(name, sci_df)
# 		# except Exception as e:
# 		# 	print(f"Error: {e}")
			

# 	html_content+="""
# 		</body>
# 		</html>
# 		"""

# 	filename= "4_27_sci"
	
# 	with open(f'/home/alex/PycharmProjects/prime-photometry-fiona/triplet_{filename}.html', 'w') as f:
# 		f.write(html_content)
# 		print(f"HTML file '{filename}' generated successfully.")


	

	
if __name__ == "__main__":
	candidates, obs_df, folders = fetch_prime_data()

	# discovery mag not deeper than 21
	# detection 

	
	# candidates = pd.read_csv("data_results/transient_candidates_results_default_recent.csv")
	# obs_df = pd.read_csv("data_results/prime_observations_recent.csv")
	
	# good_sci = ["SN2025wzh", "AT2025umy", "AT2025wdt", "SN2025abgi", "AT2025abcc", "AT2025aahv", "AT2025zwq", "SN2025yrt", "AT2025zgx","AT2025yzj", "SN2025vzq", "AT2025vma", "AT2025vdr", "AT2025uwp", "AT2025ukj",  "AT2025szq", "AT2025swv", "AT2025rkt", "AT2025smc", "AT2025rxk", "AT2025rnf", "AT2025qsr", "AT2025nqj", "AT2025gee", "SN2025hkm", "AT2025ika", "AT2025ihi", "AT2025dil", "SN2024qje", "AT2024xwo", "AT2024vuk", "AT2024vpy", "AT2024vlx", "AT2024uwv", "AT2024uwy", "AT2024uwk", "SN2024sfg", "AT2024sjl", "AT2024scc", "AT2024rwy", "AT2024qir", "AT2024dtu", "AT2024eel", "AT2024bcj", "AT2024axt", "AT2024ase"]
	# to_obs = obs_df[obs_df['full_name'].isin(good_sci)]
	# to_obs.to_csv("to_observe_4_28.csv")	
	
	# sci_df = pd.read_csv("data_results/to_observe_4_13.csv")
	folders=20 #5567
	
	# timestamp = time.time()


	candidates, triplet_df, sci_df = make_triplet_df(candidates, obs_df, folders)
	# print(sci_df)
	# # make_html_only_sci(sci_df)
	# make_html(triplet_df)
	
	# gen_util.save_summary_statistics(candidates, triplet_df)
	# generate_plots(triplet_df)


