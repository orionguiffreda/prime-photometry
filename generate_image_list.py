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

def count_cat_lines(band_path, candidates, i):
	line_count = 0 # some have multiple chips
	for chip in ["C1", "C2", "C3", "C4"]:
		astrom_path = f"{band_path}/{chip}_astrom"
		astrom_files = list_files(astrom_path)
		astrom_cat_files = [f for f in astrom_files if f.endswith(".cat")]
		line_count = 0
		for file in astrom_cat_files:
			# print("found", file, "cat file")
			with open(f"{astrom_path}/{file}", 'r', errors='ignore') as f:
				lines = f.readlines()
				line_count += len(lines)
				# if len(lines) >0:
					# print("added", len(lines), "lines from", file)
	print("cat lines:", line_count)	
	candidates.loc[i,'cat_length'] = max(line_count, candidates.loc[i,'cat_length'])


def chip_search(file_list):
	"""
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
		
def GRB_SN_Match(candidates, i):
	GRB_SN_Data = pd.read_csv('~/PycharmProjects/prime-photometry-fiona/GRBSNData/GRBSNdbdata.txt')
	# print(GRB_SN_Data.iloc[-30:])
	if candidates.loc[i, 'name'] in GRB_SN_Data['SNe']: # GRB data might start with AT/SN?
			idx = GRB_SN_Data['SNe']==candidates.loc[i, 'name']
			GRB = GRB_SN_Data[idx]['GRB']
			print(GRB)
			is_NaN = math.isnan(float(GRB)) #issue
			# candidates.loc[i, 'GRB'] = GRB if not is_NaN else ''
			print("FOUND GRB")
			return GRB if not is_NaN else ''

def fix_dates(df):
	# 'full_tns_objects_11_15.csv?
	# date_csv = "~/PycharmProjects/prime-photometry-fiona/relevant_objs_2.csv"
	# date_csv = pd.read_csv(date_csv)
	# date_csv['full_name'] = date_csv['name_prefix']+date_csv['name']
	# date_csv = date_csv[['full_name', 'discoverydate']]
	# print(df.keys())

	# df = df.drop(columns='discoverydate') # in werid format i dont understand
	# df = df.merge(date_csv, 'left', on=['full_name'])
	
	df['discoverydate'] = pd.to_datetime(df['discoverydate'],format='%Y-%m-%d %H:%M:%S.%f')
	# df['discoverydate'] = pd.to_datetime(df['discoverydate'], format='%m/%d/%y')
	df['field_date'] = pd.to_datetime(df['field'].str[-10:], format='%Y-%m-%d')
	df['days_since_discovery'] = (df['field_date'] - df['discoverydate']).dt.days


	return df

def assign_roles(group):
	group = group.copy()
	group["role"] = "ref"

	# if any mag has identified a change
	sci_idx = -1
	# print(group)
	for idx,item in group.iterrows():
		# if any magnitude method registered a detection
		# print(item)
		def detection(key):
			S_CM = float(item.get(key))
			return S_CM!=1 and not math.isnan(S_CM)
			
		if detection('autoS_CM') or detection('aperS_CM') or detection('psfS_CM'):
			# if this is a new detection, or closer detection to TNS
			if abs(item['days_since_discovery']) < 20:
				if sci_idx == -1 or item['autoMag'] < group.loc[sci_idx, 'autoMag']: # brighter in this image
					sci_idx = idx
	
	group.loc[sci_idx, "role"] = "sci"

	return group

def fetch_prime_data(transient_file = "~/PycharmProjects/prime-photometry-fiona/tns_objects_11_15.csv"):
	folders = 0
	
	id_tags= list_dirs(BASE)
	if not id_tags:
		print("No directories found under /mnt/photometry.")
		return

	candidates = pd.read_csv(transient_file)

	# add in discovery dates
	dates = pd.read_csv('~/PycharmProjects/prime-photometry-fiona/full_tns_objects_11_15.csv')
	del candidates['discoverydate']
	candidates = candidates.merge(dates[['objid', 'discoverydate']], on='objid', how='left')

	
	print(candidates['discoverydate'])
	candidate_ids = list(candidates['full_name'])

	candidates['status']=''
	candidates['GRB']=''
	triplets = pd.read_csv("data_results/triplets_recent.csv")
	triplet_names = list(triplets["full_name"])
	
	print(triplet_names)
	
	obs_rows = []
	
	for i,id_tag in enumerate(candidate_ids):
		candidates.loc[i, 'GRB'] = GRB_SN_Match(candidates, i) # clunky
					
		if id_tag not in id_tags: 
			continue
			
		id_path = f"{BASE}/{id_tag}"
		print(f"\n ID Tag: {id_tag}")

		fields = list_dirs(id_path)
		for field in fields:
			field_path = f"{id_path}/{field}"
			bands = list_dirs(field_path)

			for band in bands:
				stack_path = f"{field_path}/{band}/stack"
				try:
					files = os.listdir(stack_path)
				except IOError:
					continue

				
				folders += 1
				# get files
				cutout_files = [f for f in files if contains('Cutout', f) and f.endswith(".fits")]
				cutout_chip_file, cutout_chip = chip_search(cutout_files)

				ecsv_files = [f for f in files if f.startswith('GRB') and f.endswith('.ecsv')]
				ecsv_chip_file, ecsv_chip = chip_search(ecsv_files)

				coadd_files = [f for f in files if f.startswith("coadd") and f.endswith('.fits')]
				
				if len(coadd_files) ==0:
					print('no coadd found, data not fully processed!')
					continue
				
				if len(cutout_files) == 0:
					print('no cutout found, data not fully processed!')
					continue

				
				


				# enforce chip match if GRB cutout is labelled with a chip
				if cutout_chip:
					coadd_files = [f for f in coadd_files if contains(cutout_chip, f)]
				
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
					"chip": cutout_chip,
					"coadd_path": f"{stack_path}/{coadd_file}",
					# "has_ecsv": len(ecsv_files) > 0,
					"ecsv_path": f"{stack_path}/{ecsv}",
					"discoverydate":candidates.loc[i, 'discoverydate']

					
					# "n_cat": count_cat_lines(f"{field_path}/{band}", candidates, i),
				}


				# coadd_ecsv = [f for f in files if f.startswith('coadd') and f.endswith('.ecsv')]
				
				# if len(coadd_ecsv) >1:
				# 	keys = ["MAG_AUTO"]
				# 	row.update({key: header[key] for key in header_keys})


				

				# add in data from header

				fits_files = glob.glob(f"{field_path}/{band}/sky/*.fits")
				if len(fits_files) == 0:
					print("no sky fits files found!")
					continue
					
				fits_file = fits_files[0] # select the first one
				header = fits.getheader(fits_file)
				header_keys = ["EXPTIMEE", "DITHTYP", "DITHRAD", "DITH_REP", "NINT", "ROTOFF"]
				
				row.update({key: header[key] for key in header_keys})

				row["detection"] = []
				# add in data from ecsv
				if len(ecsv_files)>0:
					df = pd.read_csv(row["ecsv_path"], sep='\s+', comment='#')
					ecsv_data = df.to_dict(orient='records')  # list of dicts, one per row
					for ecsv_dict in ecsv_data:
						row.update(ecsv_dict)
						# print(ecsv_dict.get('autoS_CM'), type(ecsv_dict.get('autoS_CM')))

						
				print(row["detection"])
				
				obs_rows.append(row)
				print("added row!")
				year = id_tag[4:6]
				# print(year)

				
				if id_tag in triplet_names and (year=="24" or year=="25"): 

					print( "running: photometrus", "process", "stack"
								    " -astrom_only"
								    " -stack", stack_path,
								    "-chip", f"{header['CHIP']}")
					try:
						subprocess.run([
								    "photometrus", "process", "stack",
								    "-astrom_only",
								    "-stack", stack_path,
								    "-chip", f"{header['CHIP']}"
								], check=True)
					except Exception as e:
							print(f"Error: {e}")
					
		
						
					try:
						print("running photometrus photometry")
						

						subprocess.run([
							"photometrus", "photometry",
							"-stackpath",stack_path,
							"-band", band,
							"-chip", f"{header['CHIP']}",
							"-grb_ra", f"{candidates.loc[i, 'ra']}",
							"-grb_dec", f"{candidates.loc[i, 'declination']}",
							"-grb_radius", "1.0"
						])

					except Exception as e:
						print(f"Error: {e}")
						# continue

					
	# sftp.close()
	# ssh.close()

	# print(obs_rows)
	obs_df = pd.DataFrame(obs_rows)
	# print(obs_df.head())
	obs_df = combine_mag(obs_df)
	# print(obs_df.head())

	obs_df = fix_dates(obs_df)
	# print(obs_df.head())

	
	print(folders)
	obs_df.to_csv("data_results/prime_observations_recent.csv", index=False)
	candidates.to_csv("data_results/transient_candidates_results_default_recent.csv", index=False)

	return candidates, obs_df, folders

def label_edge_transients(df):
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
			if id_tag == 'AT2025xfk':
				print('saw AT2025xfk:', (id_tag, field, band))

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

	df = df[df['days_since_discovery'] < 8]
	print(len(df), "detections within 7 day window of ZTF measurement")

	return df


def make_triplet_df(candidates, obs_df, folders):

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

	# columns={ # should just automatically append _ref to these probably 
			# "coadd_path": "ref_coadd_paths",
			# "field": "ref_fields", # should all be on same grid field, this is for the date 
			# "chip":"ref_chips",
			# "has_ecsv":"detection",
			# "autoMag":"ref_auto_mags",
			# "autoMag_Err":"ref_auto_mag_err",
			# "psfMag":"ref_psf_mags",
			# "psfMag_Err":"ref_psf_mag_err",
			# "aperMag":"ref_aper_mags",
			# "aperMag_Err":"ref_aper_mag_err",
			# "days_since_discovery":"ref_days_since_discovery",
	
	# print(scis.loc[0, candidate_groups])
	# print(refs_agg.loc[0, candidate_groups])
	

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
	scis = scis.merge(
		candidates,
		how='left',
		on=["full_name"],
	)

	print(scis.keys())

	sci_only = find_single_detections(scis)
	# print("AT2025xfk in scis 2:", "AT2025xfk" in list(sci_only['full_name']))
	sci_only_minus_edge = label_edge_transients(sci_only)
	# print("AT2025xfk in scis 3:", "AT2025xfk" in list(sci_only_minus_edge['full_name']))
	# print("AT2025xfk in scis 4:", "AT2025xfk" in list(sci_only['full_name']))
	# print("sci:", len(sci_only), list(sci_only["full_name"]))	


	
	# print(sci_only.head())
	# save files 
	timestamp = time.time()
	# candidates.to_csv(f'data_results/transient_candidates_results_{timestamp}.csv', index=False)
	
	sci_only_minus_edge.to_csv(f'data_results/to_observe_{timestamp}.csv', index=False)
	print("sci only:", len(sci_only), list(sci_only["full_name"]))	


	# print(triplets["full_name"].head(5))
	# print(triplets["coadd_path"].head(5))
	# print(triplets["ref_coadd_paths"].head(5))

	def constrain_dt(df):
		df['discoverydate'] = pd.to_datetime(df['discoverydate'])
		print(len(df), "total triplets")
		start_date = datetime(2024, 1, 1, 0, 0, 0)
		df = df[df['discoverydate'] > start_date]
		print(len(df), "detections after 2023")
		df['days_since_discovery_abs'] = abs(df['days_since_discovery'])
		df = df[df['days_since_discovery_abs'] < 8]
		print(len(df), "detections within 7 day window of ZTF measurement")
		return df

	triplets_subset = constrain_dt(triplets) #UNCOMMENT
	triplets.to_csv("data_results/triplets_small_radius.csv", index=False)

	
	triplet = len(triplets)
	sci = len(scis)
	ref = len(refs)
	print("distinct transients:",len(set(list(triplets["full_name"]))))
	print(len(set(list(triplets["full_name"]))))
	print(f"triplets: {triplet}, only sci: {sci-triplet}, only ref: {ref-triplet}, partially processed: {folders-sci-ref+triplet}")
	print("\n Complete.")

	return candidates, triplets

def save_summary_statistics(candidates):

	# pie chart for data results
	
	status_counts = pd.Series(candidates['status'], name="status").value_counts()
	print(status_counts)
	status_counts = status_counts[status_counts < 10000]
	# status_counts.plot(kind='pie', autopct='%1.1f%%')


	# plt.savefig(f"data_results/status_pie_{timestamp}.png")

	# histogram for cat file length
	# cat_counts = pd.Series(candidates['cat_length'], name="cat_length").value_counts()
	# print(cat_counts)
	# cat_counts = cat_counts[cat_counts < 10000]

	candidates = candidates.where(candidates['cat_length'] < 700000, 0)

	plt.hist(candidates[candidates['cat_length'] > 0]['cat_length'] , bins=100)
	plt.legend()

	plt.savefig(f"data_results/cat_hist_{timestamp}.png")

	
def get_relevant_info(id_tag, df):


	row = df[df["full_name"]==id_tag]
	
	def get_item(key):
		return row[key].iloc[0]
		
	ra, dec = get_item('ra'), get_item('declination')
	
	
	FITS_SCI = str(row["coadd_path"].iloc[0])
	FITS_REF = str(row["ref_coadd_path"].iloc[0][0])

	sci_stack = os.path.dirname(FITS_SCI)
	ref_stack = os.path.dirname(FITS_REF)

	do_subtractions = False

	diff_path = glob.glob(f"/mnt/photometry/{id_tag}/*.fits")	
	diff_path = [f for f in diff_path if f.startswith("sub") or f.endswith(".sfftdiff.fits")]
	cat = ""
	for path in diff_path:
		cat_path = glob.glob(f"{path[:-4]}*.cat")
		if len(cat_path)>0:
			FITS_DIFF = path
			cat = cat_path[0]

	print()
	
	if id_tag=="AT2025xfk":# getting stuck
		return ""
		
	if cat=="" or do_subtractions:
		# return ""
		print("doing subtractions, extracting sources...")
		FITS_SCI, FITS_REF, FITS_DIFF, PixA_DIFF, SFFTPrepDict, sci_bkg_path, ref_bkg_path =  sfft_util.SFFT(FITS_SCI, FITS_REF)
		rms_ref, rms_sci, rms_diff_predicted, rms_diff = sfft_util.get_diff_stats(SFFTPrepDict, PixA_DIFF)
		sources = sfft_util.source_extract(FITS_DIFF)

	else: 
		print("reading in diff image and diff source catalog")
		print("catalog:", cat)
		sci_bkg_path, ref_bkg_path = glob.glob(f"{sci_stack}/background.png")[0], glob.glob(f"{ref_stack}/background.png")[0]

		# sci_bkg, sci_bkg_path = sfft_util.get_bkg(FITS_SCI)
		# ref_bkg, ref_bkg_path = sfft_util.get_bkg(FITS_REF)
		PixA_DIFF = fits.getdata(FITS_DIFF)

		rms_ref, rms_sci, rms_diff = sfft_util.rms(fits.getdata(FITS_SCI), "sci"), sfft_util.rms(fits.getdata(FITS_REF), "ref"), sfft_util.rms(PixA_DIFF, "diff")
		sources = Table.read(cat, hdu=2)


	print("sci, ref, dif:",FITS_SCI, FITS_REF, FITS_DIFF)
	diff_dir = os.path.dirname(FITS_DIFF)


	html_content = ""

	def prnt(string_content):
		return "<p>"+string_content+"</p>"
		
	
	html_content += prnt(f"\nSTATS for {id_tag} (discovery mag: {get_item('discoverymag')}, sci ref difference: {round(get_item('sci_ref_dif'),3)}, {get_item('band')} band, ra={round(ra,4)}, dec={round(dec,4)})")

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


	def print_mags(img_str, header):
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

	def GRB_png_html_full(file_path):
		with open(file_path, "rb") as f:
			img_data = base64.b64encode(f.read()).decode("utf-8")
			html_content = f'<img src="data:image/png;base64,{img_data}" width="200" height="200" alt="cutout">'
			html_content += f"<i>{prnt(os.path.basename(file_path))}</i>"
			return html_content
			

	def GRB_png_html(directory, filename):
		file_path = f"{directory}/{filename}.png"
		return GRB_png_html_full(file_path)


	
	def diff(sci_stack, ref_stack, ra, dec):
		html_content = ""
		
		good_ztf_sources = sfft_util.ztf_cuts(sources, PixA_DIFF)

		source_savename, source_threshname, closest_ra, closest_dec = sfft_util.closest_source(good_ztf_sources, FITS_DIFF, ra, dec)
		distance = util.distance(ra, dec, closest_ra, closest_dec)

		html_content += f"AB mag: {util.forced_photometry(ra, dec, FITS_DIFF)}"
		html_content += f"<p>distance to closest source: {distance} arcsec (ra={round(closest_ra,4)}, dec={round(closest_dec,4)})</p>"
		html_content += f"<p>RMS: sci={round(rms_sci,3)}, ref={round(rms_ref,3)}, diff={round(rms_diff,3)}</p>" # pred={round(rms_diff_predicted,3)}

		return html_content, GRB_png_html(diff_dir, source_savename)

	def source_info():
		html=''
		for i,item in enumerate(row['RA']):
			prime_ra, prime_dec = row['RA'].iloc[i], row['DEC'].iloc[i]
			html+= f"distance to source {util.distance(ra, dec, prime_ra, prime_dec)}"
		return html



	sci_savename, sci_threshname = util.make_cutout(FITS_SCI, ra, dec, png=True, display_file=False, name_ext="sci")
	ref_savename, ref_threshname = util.make_cutout(FITS_REF, ra, dec, png=True, display_file=False, name_ext="ref")
	diff_savename, diff_threshname = util.make_cutout(FITS_DIFF, ra, dec, png=True, display_file=False, name_ext="diff")

	diff_html_content, source_png_html = diff(sci_stack, ref_stack, ra, dec)

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
	    {append_table_row(["SCI", "REF", "DIFF"])}
	  </thead>
	  <tbody>
	    {append_table_row([GRB_png_html(sci_stack, sci_savename)+source_info(), GRB_png_html(ref_stack, ref_savename), GRB_png_html(diff_dir,diff_savename)])}
		{append_table_row([GRB_png_html_full(sci_bkg_path), GRB_png_html_full(ref_bkg_path), source_png_html])}
	    {append_table_row([print_mags("sci",fits.getheader(FITS_SCI)) , print_mags("ref",fits.getheader(FITS_REF)), diff_html_content])}
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

	df["mag_dif"]= 10
	df["sci_ref_dif"] = 10
	df["average_sci_mag"] = 10
	for i, row in df.iterrows():
		df.loc[i,"mag_dif"] = abs(float(row['psfMag']) - float(row['autoMag']))
		mag_keys = ["autoMag", "aperMag"] # not psf?
		sci_ref_difs = [float(row[key]) - np.nanmean(row[f'ref_{key}']) for key in mag_keys]
		df.loc[i,"sci_ref_dif"] = np.nanmean(sci_ref_difs)
		df.loc[i,"average_sci_mag"] = np.nanmean([float(row[key]) for key in mag_keys])

	# df = df[df["mag_dif"] <= 0.5 ] #1] for ZTF #UNCOMMENT

	return df

	
	
def process_ecsv_data(triplet_df):

	# triplet_df = combine_mag(triplet_df)
	# combine mag into one col?
	print(triplet_df.keys())
	print("num triplets with small dt",len(triplet_df))

	best_candidates = apply_ZTF_cuts(triplet_df)
	# best_candidates = best_candidates[best_candidates["full_name"].isin(["SN2025adcv", "SN2024xyu", "AT2025lbx","AT2025uvl", "AT2025aayv", "AT2024xwh", "SN2025adcv", "AT2025bob"])]
	# best_candidates = best_candidates.head(3)
	
	best_candidates['days_since_discovery_abs'] = abs(best_candidates['days_since_discovery'])
	best_candidates = best_candidates[best_candidates['days_since_discovery_abs']<20]
	# best_candidates = best_candidates.sort_values('average_sci_mag', ascending=False).head(100)
	best_candidates = best_candidates.sort_values('sci_ref_dif', ascending=True)

	# detections=["AT2025bnx", "SN2024aeqs", "AT2025gwl", "AT2025bte", "AT2024ahg", "AT2024acdq", "AT2025gex", "SN2025ihl","AT2025sox", "AT2025ably", "AT2025rpv"]
	# neg_residuals=["AT2024aevn", "AT2025bnu", "AT2024dfx", "AT2025pwa", "AT2025bnv", "AT2025fdc", "AT2024svc", "AT2024zoz", "AT2025vjl", "AT2024fje", "AT2024xwh", "AT2025wog", "AT2024yfr", "AT2024yfg"]
	
	# best_candidates = best_candidates[best_candidates["full_name"].isin(detections)]
	

	print(best_candidates[['full_name', 'days_since_discovery', 'sci_ref_dif']])
	# best_candidates = pd.concat([best_candidates,best_candidates_picked])

	
	# print(best_candidates[['full_name', 'SNR', 'Elongation', 'Mag', 'apMag2', 'autoMag', 'psfMag', 'coadd_path', 'ref_coadd_path']].head(20))


	html_content = """
	<!DOCTYPE html>
	<html>
	<head>
		<link rel="stylesheet" href="styles.css">
	</head>
	<style>
	body {background-color: powderblue;}
	h1   {color: black;}
	p    {color: black;}
	
	#wrapper {width: 1100px; border: 1px solid black, margin;15px}
	#first {width: 300px; float:left; padding:15px;}
	#second {width: 300px; float:left; padding:15px;}
	#third {padding:15px;}
	</style>
	
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

	filename= "detections"
	
	with open(f'/home/alex/PycharmProjects/prime-photometry-fiona/triplet_{filename}.html', 'w') as f:
	    f.write(html_content)
	    print(f"HTML file '{filename}' generated successfully.")

	
if __name__ == "__main__":
	candidates, obs_df, folders = fetch_prime_data()
	# candidates = pd.read_csv("data_results/transient_candidates_results_default_recent.csv")
	# obs_df = pd.read_csv("data_results/prime_observations_recent.csv")
	# folders=6435
	timestamp = time.time()


	candidates, triplet_df = make_triplet_df(candidates, obs_df, folders)
	process_ecsv_data(triplet_df)



	
	# save_summary_statistics(candidates, triplet_df)
	# generate_plots(triplet_df)




