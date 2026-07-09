import matplotlib.pyplot as plt
import math
import numpy as np
import pandas as pd
import os
from pathlib import Path
import time
import glob

from astropy.io import fits
from astropy.table import Table

import photometrus_utils as util
import sfft_util
import assess_processed_data_util as gen_util
from assess_processed_data_util import list_dirs, list_files, contains, combine_mag, constrain_dt

def prnt(string_content):
	"""
	create a string within a paragraph
	"""
	return "<p>"+string_content+"</p>"
		

def print_lim_mag(header, img_str):
	"""
	print limiting magnitudes from header
	"""

	html_content = ""

	try:
		html_content += prnt(f"{img_str} (psf, auto) lim mag: {float(header['HIERARCH lim_mag_psf'])}, {float(header['HIERARCH lim_mag_auto'])}")
	except KeyError as e:
		try:
			html_content += prnt(f"{img_str} (psf, auto) lim mag: {float(header['LM_PSF'])}, {float(header['LM_AUTO'])}")
		except KeyError as e2:
			html_content += prnt(f"Error:{e}, {e2}")
	return html_content


def append_table_row(contents):
	""" 
	Add a row of contents to the table displacing a certain transient
	"""
	html_content = "<tr>"
	for item in contents:
		html_content +=f"<td>{item}</td>"
	html_content += "<tr>"
	return html_content


def diff(sci_stack, ref_stack, FITS_DIFF, PixA_DIFF, ra, dec, sources, rms_sci, rms_ref, rms_diff, px):
	"""
	Get HTML information for the difference image
	"""
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
	
	good_ztf_sources = sfft_util.ztf_cuts(sources, PixA_DIFF, crop=300)
	source_savename, source_threshname, closest_ra, closest_dec = sfft_util.closest_source(sources, FITS_DIFF, ra, dec)
	distance = util.distance(ra, dec, closest_ra, closest_dec)
	diff_savename, diff_threshname = util.make_cutout(FITS_DIFF, ra, dec, png=True, display_file=False, name_ext=name_ext)

	
		
	html_content += f"<p>AB mag: {mag}</p><p>AB mag flipped: {flipped_mag}</p>"

	html_content += f"<p>distance to closest source: {distance} arcsec (ra={round(closest_ra,4)}, dec={round(closest_dec,4)})</p>"
	html_content += f"<p>RMS: sci={round(rms_sci,3)}, ref={round(rms_ref,3)}, diff={round(rms_diff,3)}</p>" # pred={round(rms_diff_predicted,3)}
	# html_content += f"<p>&mu;x={udx}, &sigma;x={sdx}, &mu;x={udy}, &sigma;x={sdy}</p>"
	html_content += f"{len(good_ztf_sources)} good / {len(sources)} raw sources in diff image"
	html_content += f"<p>pixels from edge: {px}</p>"

	return html_content, diff_savename


		
def get_transient_html(id_tag, df):
	"""
	Retrieve one row of relevant information about a triplet for transient HTML webpage
	"""

	row = df[df["full_name"]==id_tag]
	
	def get_item(key):
		return row[key].iloc[0]
		
	ra, dec = get_item('ra'), get_item('declination')
	
	
	FITS_SCI = str(row["coadd_path"].iloc[0])
	FITS_REF = str(row["ref_coadd_path"].iloc[0][2:-2])

	sci_stack = os.path.dirname(FITS_SCI)
	ref_stack = os.path.dirname(FITS_REF)

	sci_cat, sci_cat_png, fwhm_sci = sfft_util.get_cat_old(FITS_SCI)
	ref_cat, ref_cat_png, fwhm_ref = sfft_util.get_cat_old(FITS_REF)


	do_subtractions = False
		
	diff_path = glob.glob(f"/mnt/photometry/{id_tag}/*.fits")	


	diff_path = [f for f in diff_path if os.path.basename(f).startswith("sub_6_24")] # or f.endswith(".sfftdiff.fits")]
	cat = ""
	for path in diff_path:
		cat_path = glob.glob(f"{path[:-4]}*.cat")
		if len(cat_path)>0:
			FITS_DIFF = path
			cat = cat_path[0]

	
	if id_tag=="AT2025xfk":# getting stuck
		return ""
		
	if do_subtractions or cat=="":

		# print('SKIPPING UNSUBTRACTED FIELD')
		# return ""
		
		print("doing subtractions, extracting sources...")
		FITS_SCI, FITS_REF, FITS_DIFF, PixA_DIFF, SFFTPrepDict, sci_bkg_path, ref_bkg_path =  sfft_util.SFFT(FITS_SCI, FITS_REF)
		rms_ref, rms_sci, rms_diff_predicted, rms_diff = sfft_util.get_diff_stats(SFFTPrepDict, PixA_DIFF)
		udx, sdx, udy, sdy = sfft_util.sfft_source_reg_gen(SFFTPrepDict)
		sources = sfft_util.source_extract(FITS_DIFF)

	else: 
		print()

		print("reading in diff image and diff source catalog")
		print("catalog:", cat)
		
		PixA_DIFF = fits.getdata(FITS_DIFF)

		rms_sci, rms_ref, rms_diff = sfft_util.rms(fits.getdata(FITS_SCI), "sci"), sfft_util.rms(fits.getdata(FITS_REF), "ref"), sfft_util.rms(PixA_DIFF, "diff")
		sources = Table.read(cat, hdu=2)


	print("sci, ref, dif:",FITS_SCI, FITS_REF, FITS_DIFF)
	diff_dir = os.path.dirname(FITS_DIFF)


	html_content = ""

	html_content += prnt(f"\nSTATS for {id_tag} (discovery mag: {get_item('discoverymag')}, sci ref difference: {round(get_item('sci_ref_dif'),3)}, {get_item('band')} band, ra={round(ra,4)}, dec={round(dec,4)}, chip={get_item('chip')})")



	def print_sci_ref_info(img_str, image):
		"""
		Print data features for sci and ref images 
		"""
		header = fits.getheader(image)
		html_content = ""
		def get_item(key):
			if img_str == "sci":
				return row[key].iloc[0] # get instead of idx?
			else:
				return row[f"ref_{key}"].iloc[0][1:-1]
			
		try:

			CM_list = [get_item('autoS_CM'), get_item('aperS_CM'), get_item('psfS_CM')]
			CM_arr = np.array(CM_list).flatten()
			color = "yellow" if 2 in CM_arr else "green"
			color = "cyan" if 3 in CM_arr else color
			
			open_span, close_span = f"<span style=\"background-color:{color};\">", "</span>"
			html_content += prnt(f"{open_span}Survey Crossmatch (auto, aper, psf): {get_item('autoS_CM')}, {get_item('aperS_CM')}, {get_item('psfS_CM')}{close_span}")
			
			mag_arr = np.array([float(get_item('autoMag')),float(get_item('aperMag')),float(get_item('psfMag'))]).flatten()
			avg_mag = round(np.average(mag_arr),3)
			stdev = round(np.std(mag_arr),3)
			html_content += prnt(f"<b>{img_str} mag: &mu;={avg_mag}, &sigma;={stdev}</b> \n(auto: {get_item('autoMag')}+-{get_item('autoMag_Err')}, aper:  {get_item('aperMag')} +- {get_item('aperMag_Err')}, psf: {get_item('psfMag')} +- {get_item('psfMag_Err')})")

			html_content += print_lim_mag(header, img_str)
			html_content += prnt(f"{img_str}: {get_item('field')} | {get_item('days_since_discovery')} days since discovery")
			

		except KeyError as e:
			html_content += prnt(e)

		html_content += prnt("")
		return html_content


	def source_info():
		"""
		Print distance to sources in images
		"""
		html='distance to source: '
		for i,item in enumerate(row['RA']):
			prime_ra, prime_dec = row['RA'].iloc[i], row['DEC'].iloc[i]
			html+= f"{util.distance(ra, dec, prime_ra, prime_dec)},"
		return html

	px = round(get_item('px_from_edge'), 3)
	diff_html_content, diff_savename = diff(sci_stack, ref_stack, FITS_DIFF, PixA_DIFF, ra, dec, sources, rms_sci, rms_ref, rms_diff, px)

	
	sci_savename, sci_threshname = util.make_cutout(FITS_SCI, ra, dec, png=True, display_file=False, name_ext="sci")
	ref_savename, ref_threshname = util.make_cutout(FITS_REF, ra, dec, png=True, display_file=False, name_ext="ref")
	# diff_savename, diff_threshname = util.make_cutout(FITS_DIFF, ra, dec, png=True, display_file=False, name_ext="diff")
	diff_bigger_savename, diff_bigger_threshname = util.make_cutout(FITS_DIFF, ra, dec, png=True, display_file=False, name_ext="diff_bigger", photoDistThresh=30.0)

	sci_whole_savename = util.display(FITS_SCI, png=False, show=False, savename="sci")
	ref_whole_savename = util.display(FITS_REF, png=False, show=False, savename="ref")
	diff_whole_savename = util.display(FITS_DIFF, png=False, show=False, savename="diff")

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
		{append_table_row([gen_util.GRB_png_html(sci_stack, sci_savename)+source_info(), gen_util.GRB_png_html(ref_stack, ref_savename), gen_util.GRB_png_html(diff_dir,diff_savename)])}

		{append_table_row([gen_util.GRB_png_html_full(sci_whole_savename), gen_util.GRB_png_html_full(ref_whole_savename), gen_util.GRB_png_html_full(diff_whole_savename)])}
		
		{append_table_row([gen_util.GRB_png_html_full(sci_cat_png)+f"fwhm:{fwhm_sci}", gen_util.GRB_png_html_full(ref_cat_png)+f"fwhm:{fwhm_ref}", gen_util.GRB_png_html(diff_dir, diff_bigger_savename)])}
		
		{append_table_row([print_sci_ref_info("sci",FITS_SCI) , print_sci_ref_info("ref",FITS_REF), diff_html_content])}
	  </tbody>
	</table>
	"""	

	return html_content
		

	
def get_transient_html_sci_only(id_tag, df):
	"""
	Retrieve one row of relevant information about a single epoch for transient HTML webpage
	"""

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

			html_content += print_lim_mag(header, img_str)
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




def get_html_header():
	"""
	Get header for html file
	"""
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
	
	
def make_html(df, triplet=True):
	"""
	Select transients of interest, then iterate through them (triplets or single epochs) to make a HTML website 
	"""

	print(df.keys())
	print("num triplets with small dt",len(df))
	

	if triplet:	
		ids=["AT2025lnp", "SN2024aeqs", "AT2025aayv", "AT2025lbx", "AT2024acdq", "AT2025zee"]
		best_candidates = df[df["full_name"].isin(ids)]
		best_candidates = best_candidates.sort_values('sci_ref_dif', ascending=True)
		print(best_candidates[['full_name', 'days_since_discovery', 'sci_ref_dif']])
		# print(best_candidates[['full_name', 'SNR', 'Elongation', 'Mag', 'apMag2', 'autoMag', 'psfMag', 'coadd_path', 'ref_coadd_path']].head(20))

	else:
		good_sci = ["SN2025wzh", "AT2025umy", "AT2025wdt"]
		# good_sci = ["SN2025wzh", "AT2025umy", "AT2025wdt", "SN2025abgi", "AT2025abcc", "AT2025aahv", "AT2025zwq", "SN2025yrt", "AT2025zgx","AT2025yzj", "SN2025vzq", "AT2025vma", "AT2025vdr", "AT2025uwp", "AT2025ukj",  "AT2025szq", "AT2025swv", "AT2025rkt", "AT2025smc", "AT2025rxk", "AT2025rnf", "AT2025qsr", "AT2025nqj", "AT2025gee", "SN2025hkm", "AT2025ika", "AT2025ihi", "AT2025dil", "SN2024qje", "AT2024xwo", "AT2024vuk", "AT2024vpy", "AT2024vlx", "AT2024uwv", "AT2024uwy", "AT2024uwk", "SN2024sfg", "AT2024sjl", "AT2024scc", "AT2024rwy", "AT2024qir", "AT2024dtu", "AT2024eel", "AT2024bcj", "AT2024axt", "AT2024ase"]
		
		best_candidates = df[df["full_name"].isin(good_sci)]



	html_content = get_html_header() 
	html_content += """
	<body>
		<h1>Sci-Ref pairs currently processed with updated photometry methods</h1>
	"""
	
	for name in best_candidates['full_name']:
		# try:
		if triplet:
			html_content += get_transient_html(name, df)
		else: # sci only
			html_content += get_transient_html_sci_only(name, df)
		# except Exception as e:
		#     print(f"Error: {e}")
			

	html_content+="""
		</body>
		</html>
		"""

	filename= "7_2_test"

	name_str = "triplet" if triplet else "sci"
	
	with open(f'/home/alex/PycharmProjects/prime-photometry-fiona/{name_str}_{filename}.html', 'w') as f:
		f.write(html_content)
		print(f"HTML file '{filename}' generated successfully.")


	
	
if __name__ == "__main__":

	triplet_df = pd.read_csv("data_results/triplets.csv")
	sci_df = pd.read_csv("data_results/to_observe_7_7.csv")
	
	# make_html(triplet_df)
	make_html(sci_df, triplet=False)
	