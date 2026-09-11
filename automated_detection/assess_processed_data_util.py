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
	"""Determine if a string contains """
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
	"""
	Calculate galactic coordinates from ra and declination
	"""
	coords = SkyCoord(
		ra=df['ra'].to_numpy() * u.deg,
		dec=df['declination'].to_numpy() * u.deg,
		frame='icrs'
	)

	# coords = SkyCoord(
	# 	ra=df['ra'].values * u.deg,
	#	dec=df['declination'].values * u.deg,
	#	frame='icrs'
	# )
	
	df['l'] = coords.galactic.l.deg
	df['b'] = coords.galactic.b.deg
	return df


	


def GRB_png_html_full(file_path):
	"""
	Create a row of pngs for dashboard from a full image filapath
	"""
	with open(file_path, "rb") as f:
		img_data = base64.b64encode(f.read()).decode("utf-8")
		html_content = f'<img src="data:image/png;base64,{img_data}" width="200" height="200" alt="cutout">'
		html_content += f"<i><p>{(os.path.basename(file_path))}</p></i>"
		return html_content
			

def GRB_png_html(directory, filename):
	"""
	Create a row of pngs for dashboard from a directory name
	"""
	file_path = f"{directory}/{filename}.png"
	return GRB_png_html_full(file_path)




######## unused #######
def count_cat_lines(band_path, candidates, i):
	"""
	Count lines in all catalog files
	"""
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

