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

from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u
import astropy.time
import astropy.coordinates

import photometrus.photometry.photometry as photo
import photometrus_utils as util
import sfft_util


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
		df = df[df['days_since_discovery_abs'] < 20]
		print(len(df), "detections within 20 day window of ZTF measurement")
		return df


def add_galactic_coords(df):

	coords = SkyCoord(
		ra=df['ra'].values * u.deg,
	   dec=df['declination'].values * u.deg,
	   frame='icrs'
	)
	
	df['l'] = coords.galactic.l.deg
	df['b'] = coords.galactic.b.deg
	return df


	# print(df.loc[0,'field'], df.loc[0,'field'][-10:])
	# df['date'] = [field[-10:] for field in df['field']]
	# for i, row in df.iterrows():
	#	coord = astropy.coordinates.TETE(
	#		ra=row['ra'] * u.deg,
	#		dec=row['declination'] * u.deg,
	#		obstime=astropy.time.Time(row['date'])
	#	)
		
	#	galactic_coord = coord.transform_to(astropy.coordinates.Galactic())
	#	df.loc[i,'l'] = galactic_coord.l.deg
	#	df.loc[i,'b'] = galactic_coord.b.deg
	
	return df


def get_latitude(ra, dec, field):
	date = field[-10:]

	coord = astropy.coordinates.TETE(
		ra=ra * u.deg,
		dec=dec * u.deg,
		obstime=astropy.time.Time(date)
		)
		
	galactic_coord = coord.transform_to(astropy.coordinates.Galactic())
	
	return galactic_coord.b.deg

	


def GRB_png_html_full(file_path):
	with open(file_path, "rb") as f:
		img_data = base64.b64encode(f.read()).decode("utf-8")
		html_content = f'<img src="data:image/png;base64,{img_data}" width="200" height="200" alt="cutout">'
		html_content += f"<i><p>{(os.path.basename(file_path))}</p></i>"
		return html_content
			

def GRB_png_html(directory, filename):
	file_path = f"{directory}/{filename}.png"
	return GRB_png_html_full(file_path)




######## unused #######
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