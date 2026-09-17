import numpy as np
import math
import glob
import pandas as pd
import torch.nn.functional as F
import ast 
import os
import re
import os.path as pa
from pathlib import Path

import matplotlib.pyplot as plt

import astropy
from astropy.io import fits, ascii
from astropy.table import Table, vstack
from astropy.coordinates import Angle, SkyCoord
from astropy.wcs import WCS, utils
from astropy.stats import sigma_clip, sigma_clipped_stats
import astropy.units as u
from astropy.nddata import Cutout2D

from regions import CircleSkyRegion


import sfft_util
import photometrus_utils as util
from photometrus.astrom.astrom_img_sub import multi_epoch_astrom


def augment_image(images, labels):
	"""
	Augment training data by flipping and rotating images
	"""
	def flip_vertical(image):
		return np.stack([np.flipud(ch) for ch in image])

	def flip_horizontal(image):
		return np.stack([np.fliplr(ch) for ch in image])
	
	augmented_images = []
	augmented_labels = []
	for i, myimage in enumerate(images):
		for k in range(4):
			rotated = np.stack([np.rot90(ch, k=k) for ch in myimage])
			augmented_images.append(rotated)
			augmented_images.append(flip_vertical(rotated))  # vertical flip
			augmented_images.append(flip_horizontal(rotated))  # horizontal flip
			augmented_images.append(flip_horizontal(flip_vertical(rotated)))  # both flips
		augmented_labels.extend([labels[i]] * 16)
	return augmented_images, augmented_labels


def pad_image(tensor_batch, px=63, resize=False, pad=True):
	"""
	Currently unused functionality for padding images
	"""
	if resize:
		resized_tensor = F.interpolate(tensor_batch, size=(px, px), mode='bilinear', align_corners=False)
		tensor_batch = resized_tensor.squeeze(0)
	if pad:
		print ("Shape before padding:", tensor_batch.shape)
		d1, d2, d3, d4 = tensor_batch.shape
		def get_px(d):
			half_px = (px-d)/2
			return math.floor(half_px), math.ceil(half_px)
		tensor_batch = F.pad(tensor_batch, list(get_px(d3) + get_px(d4)), value = 10**-9)
	return tensor_batch


def cutout(file_list, ra, dec, id_tag, photoDistThresh=8):
	"""
	Read in image data for deep learning
	"""
	cutouts = []
	for file in file_list:
		coords = SkyCoord(ra=[ra], dec=[dec], frame='icrs', unit='degree')
		threshname, savename, cutout = grb_cutout(file, coords, id_tag, photoDistThresh=8)
		cutouts.append(cutout.data)
	return cutouts


def grb_cutout(imageName, GRBcoords, id_tag, photoDistThresh, regprimename=None, regsurvname=None):
	"""
	Return a cutout of an image with specified coordinates and radius
	"""
	imgdata = fits.getdata(imageName)
	img = fits.open(imageName)
	head = img[0].header
	w = WCS(head)

	field = re.search(r'field\d\d\d\d\d-\d\d\d\d-\d\d-\d\d', imageName)
	if not field:
		field = re.search(r'field\d\d\d\d-\d\d\d\d-\d\d-\d\d', imageName)
	if field:
		field = field.group()
	elif re.search('dif', imageName):
		field = 'dif'
		
	save_dir = f'image_data/{id_tag}' #'%s_%s_Cutout_%s_%s' % (grb_name, band, survey, num)
	os.makedirs(save_dir, exist_ok=True)
	
	savename = f'{save_dir}/{field}_Cutout'

	threshname = 'GRB_query_thresh.reg' #%s_query_thresh.reg' % grb_name
	
	size = 4 * photoDistThresh * u.arcsec
	try:
		cutout = Cutout2D(imgdata, GRBcoords, size, wcs=w, copy=True)
		region = CircleSkyRegion(center=GRBcoords[0], radius=Angle(photoDistThresh, unit='arcsec'))
		pix_region = region.to_pixel(cutout.wcs)
	except astropy.nddata.utils.NoOverlapError:
		print(' Area of GRB threshold not found within image, cannot generate cutout!')
		return savename, threshname, None
	
	mean, median, sigma_cut = sigma_clipped_stats(cutout.data)
	plt.figure(10, figsize=(8, 8))
	plt.imshow(cutout.data, vmin=median - 3 * sigma_cut, vmax=median + 3 * sigma_cut, origin='lower', cmap='viridis')
	pix_region.plot(color='cyan', ls='--', label='Input GRB threshold')
	
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
	plt.savefig(savename + '.png', dpi=300)
	plt.clf()
	
	fits.writeto(savename + '.fits', cutout.data, cutout.wcs.to_header(), overwrite=True)
	return savename, threshname, cutout



def coadd_to_sources(triplet):
	"""
	Read in all the sources from a single triplet to pass through classifier, evaulating images as we would during automated detection. 
	"""
	print(triplet['ref_coadd_path'])
	print(ast.literal_eval(triplet['ref_coadd_path']))
	ref = ast.literal_eval(triplet['ref_coadd_path'])[0]
	sci = triplet['coadd_path']

	diff_file = triplet["sfft_diff"]
	diff_data = fits.getdata(diff_file)
	diff_cat_name = triplet["diff_cat"]
	
	diff_cat = Table.read(diff_cat_name, hdu=2)
	cat_good = diff_cat[sfft_util.get_good_mask(diff_cat)]

	cat_good_df = sfft_util.make_df(cat_good)
	band = triplet['band']

	# print("raw cat:",len(diff_cat), "good cat",len(cat_good))

	# crop
	xdim, ydim = np.shape(diff_data)
	cat_good_cropped = sfft_util.crop_mask(cat_good_df, xdim, ydim, crop=1000)
	cat_good_df = cat_good_cropped[cat_good_cropped["SNR_WIN"]>10] 
	
	print("raw cat:",len(diff_cat), "good cat",len(cat_good), "cropped cat", len(cat_good_cropped), "snr > 10 ", len(cat_good_df))

	images = []
	for i, source in cat_good_df.iterrows():
		print("found source at", source['ALPHA_J2000'], source['DELTA_J2000'])
		ra, dec = source['ALPHA_J2000'], source['DELTA_J2000']
		triplet_imgs = cutout([triplet['coadd_path'], ast.literal_eval(triplet['ref_coadd_path'])[0], triplet['sfft_diff']], ra, dec, triplet['full_name'])
		images.append(np.nan_to_num(triplet_imgs))

	return images


def n_neg(path, ra, dec):
	data = fits.getdata(path)
	wcs = WCS(fits.getheader(path))
	xpx, ypx = wcs.world_to_pixel(SkyCoord(ra=ra*u.deg, dec=dec*u.deg))
	x, y = round(xpx.item()), round(ypx.item())
	stamp = data[y-3:y+3, x-3:x+3]	
	n_neg = np.sum(stamp < 0)
	return n_neg


	
def add_analysis_fields(triplet_df):
	"""
	handle some rows with nans, some rows with real values?
	"""

	def nan_test(path):
		try:
			nan = math.isnan(float(path))
		except:
			nan=False
		return nan
	
	for i, row in triplet_df.iterrows():
	# for ra, dec, fits_diff in zip(triplet_df['ra'], triplet_df['declination'], triplet_df['sfft_diff']):
		nan = nan_test(row['sfft_diff'])
		if nan:
			print("no sfft diff path?")
			triplet_df.loc[i,'diff_mag'] = math.nan
			triplet_df.loc[i,'n_neg'] = math.nan
		else: 
			try:
				triplet_df.loc[i,'diff_mag'] = util.forced_photometry(row['ra'], row['declination'], row['sfft_diff'])
				triplet_df.loc[i,'n_neg'] = n_neg(row['sfft_diff'], row['ra'], row['declination'])
			
				sci_mag = np.nanmedian(row[['autoMag', 'aperMag', 'psfMag']])
			except ValueError:
				print("Issue with nan or inf positions, wcs error?")
				print(row['ra'], row['declination'], row['sfft_diff'])
	
			nan_sci = nan_test(sci_mag)
			if nan_sci: 
				sci_mag = util.forced_photometry(row['ra'], row['declination'], row['coadd_path'])
			triplet_df.loc[i,'sci_mag'] = sci_mag

			# print(row['sfft_diff'])
   #		  print("dif, nneg, sci mag", triplet_df.loc[i,'diff_mag'], triplet_df.loc[i,'n_neg'], sci_mag)
   #		  util.make_cutout(row['sfft_diff'], row["RA"], row["DEC"], png=True, display_file=True)
			
	triplet_df['psf_minus_aper'] = triplet_df['psfMag'] - triplet_df['aperMag']
	triplet_df = triplet_df[triplet_df['sfft_diff']].isnan()
	
	return triplet_df



def display_confusion_images(all_preds, all_labels, all_images):
	"""
	Display sample images which as TP, FP, TN, and FN
	"""
	all_preds = np.array(all_preds)
	all_labels = np.array(all_labels)
	categories = {
		"True Real (0,0)": np.where((all_labels == 0) & (all_preds == 0))[0],
		"False Bogus (0,1)": np.where((all_labels == 0) & (all_preds == 1))[0],
		"False Real (1,0)": np.where((all_labels == 1) & (all_preds == 0))[0],
		"True Bogus (1,1)": np.where((all_labels == 1) & (all_preds == 1))[0],
	}
	channel_titles = ['Science', 'Reference', 'Difference'] 

	rows_per_category = [min(len(idxs), max_examples_per_category) for idxs in categories.values()]
	total_rows = sum(max(r, 1) for r in rows_per_category)  # reserve 1 row even if empty, for the label

	fig, axes = plt.subplots(total_rows, 3, figsize=(9, 3 * total_rows))
	if total_rows == 1:
		axes = axes.reshape(1, 3)

	current_row = 0
	for title, idxs in categories.items():
		n_shown = min(len(idxs), max_examples_per_category)

		if n_shown == 0:
			for col in range(3):
				axes[current_row, col].axis('off')
			axes[current_row, 0].set_ylabel(f"{title}\n(none found)", rotation=0,
											 labelpad=60, fontsize=10, va='center')
			current_row += 1
			continue

		chosen_idxs = list(idxs[:n_shown])
		for i, idx in enumerate(chosen_idxs):
			img = all_images[idx]
			plot_row = current_row + i
			for col in range(3):
				axes[plot_row, col].imshow(img[col].numpy(), cmap='gray')
				axes[plot_row, col].set_xticks([])
				axes[plot_row, col].set_yticks([])
				if plot_row == 0 or i == 0:
					axes[plot_row, col].set_title(channel_titles[col])
			# only label the first row of this category's block
			axes[plot_row, 0].set_ylabel(title if i == 0 else "", rotation=0,
										  labelpad=60, fontsize=10, va='center')
		current_row += n_shown

	plt.tight_layout()
	plt.show()


def display_triplet_images(triplet_images):
	"""
	Display read in image data
	"""
	fig, axes = plt.subplots(len(triplet_images), 3, figsize=(9, 3 * len(triplet_images)))
	
	for i, img in enumerate(triplet_images):
		for col in range(3):
			axes[i, col].imshow(img[col], cmap='gray')
			axes[i, col].set_xticks([])
			axes[i, col].set_yticks([])
	
		title = triplet_df.iloc[i]["full_name"] + str(triplet_df.iloc[i]["bogus"])
		axes[i, 0].set_ylabel(title, rotation=0,
									  labelpad=60, fontsize=10, va='center')
	
	plt.tight_layout()
	plt.show()



def ML_prep(triplets):
	"""
	Add relevant analysis fields and enforce detection in real examples
	"""
	bogus_df = triplets[triplets["bogus"]==1].copy()
	real_df = real_df[real_df['distance'] < 4]# pd.read_csv("data_results/triplets_diff_detected.csv")
	print("bogus triplets, real triplets:",len(bogus_df), len(real_df))
	
	triplets = pd.concat([real_df, bogus_df])
	triplet_df = add_analysis_fields(triplets.copy()) 

	return triplet_df





	
	
## FROM ML ANALYSIS: source distance stuff. ADD TO HERE AND SAVE!!?
	 # triplets[triplets["bogus"]==0].copy()

# diff_distance = []
# for real in real_df.itertuples():
#	 print("finding min distance for source", real.full_name)
#	 ra, dec = real.ra, real.declination
#	 try:
#		 sources = Table.read(real.diff_cat, hdu=2)
#	 except:
#		 print("error with table format")
#		 diff_distance.append(100)

#		 continue
#	 # print(len(sources))	
#	 sources["distance"] = [util.distance(source['ALPHA_J2000'], source['DELTA_J2000'], ra, dec) for source in sources]	
#	 # min_source_idx = sources["distance"].argmin()
#	 # source = sources[min_source_idx]
#	 diff_distance.append(min(sources["distance"]))
#	 print("min_distance:", min(sources["distance"]))


# real_df["diff_distance"] = diff_distance
# real_df = real_df[real_df["diff_distance"]<4]




# real_df["diff_distance"] = diff_distance
# real_df = real_df[real_df["diff_distance"]<4]


# # dubious_reals = ["AT2025ably", "AT2025adoj","AT2025admu", "AT2025admo", "AT2025adkn", "AT2025adgu", "SN2025tpu", "AT2025ssh", "AT2025kqj", 
# #				  "AT2025gwd", "AT2025fof", "AT2025efu","AT2024aht","AT2025dbp","TDE2025chm","AT2024ahec","AT2025bti","AT2025btg",
# #				  "SN2024aeqs","AT2025admd"]
# # real_df = real_df[~real_df["full_name"].isin(dubious_reals)]

# print("bogus triplets, real triplets:",len(bogus_df), len(real_df))
# bogus_df.head()

# triplets = pd.concat([real_df, bogus_df])


# print(list(triplets.keys()))

# real_df.to_csv("data_results/triplets_diff_detected.csv")

# triplets.head()



# triplet_df['diff'] = [get_diff(id_tag) for id_tag in triplet_df['full_name']]
# # triplet['diff_mag'] = [util.forced_photometry(triplet['ra'], triplet['declination'], triplet['sfft_diff']) for i, triplet in triplets.iterrows()]
# triplet_df['n_neg'] = [n_neg(fits_diff, ra, dec) for ra, dec, fits_diff in zip(triplets['ra'], triplets['declination'], triplets['sfft_diff'])]

# triplet_df['sci_rms'] = [bkg_rms(fits.getdata(image_path)) for image_path in triplets["coadd_path"]]
# triplet_df['ref_rms'] = [bkg_rms(fits.getdata(image_path)) for image_path in triplets["ref_coadd_path"]]

# triplet_df['sci_cat_len'] = [len(get_cat(image_path)) for image_path in triplets["coadd_path"]]
# triplet_df['ref_cat_len'] = [len(get_cat(image_path)) for image_path in triplets["ref_coadd_path"]]

# triplet_df['sci_seeing'] = [sfft_util.seeing(get_cat(image_path)) for image_path in triplets["coadd_path"]]
# triplet_df['ref_seeing'] = [sfft_util.seeing(get_cat(image_path)) for image_path in triplets["ref_coadd_path"]]

# triplet_df['ref_mag'] = triplet_df[['ref_autoMag', 'ref_aperMag', 'ref_psfMag']].mean(axis=1) # not right data structure
# triplets['lim_mag'] = [bogus.limiting_mag(fits.getheader(path)) for path in triplets['coadd_path']]





# def image_heatmap(model)

# 	print("Heatmap for input image")
# 	model = net.eval()
# 	cam_extractor = GradCAM(model, target_layer='conv2')
# 	myinput, label = next(iter(test_loader))
	
# 	out = model(myinput)
# 	class_idx = out[0].argmax().item()
	
# 	activation_map = cam_extractor(class_idx, out)
# 	cam = activation_map[0].sum(dim=0)
# 	cam = F.relu(cam)
# 	cam = cam / cam.max()
	
# 	cam_np = cam.detach().cpu().numpy()
# 	plt.imshow(cam_np, cmap='jet')
# 	plt.axis('off')
# 	plt.show()
	
	
# 	# img_np = myinput[0].detach().cpu().numpy()
# 	# plt.imshow(img_np, cmap='jet')
# 	# plt.axis('off')
# 	# plt.show()
	
	
# 	# img_pil = to_pil_image(myinput[0].detach().cpu())
# 	# cam_pil = to_pil_image(cam.detach().cpu(), mode='F')
# 	# result = overlay_mask(img_pil, cam_pil, alpha=0.5)
	
# 	# plt.imshow(result)
# 	# plt.axis('off')
# 	# plt.show()






### OTHER NET ARCH



# class WinterNet(nn.Module):
#	 def __init__(self):
#		 super().__init__()
#		 self.kernel = 3
#		 self.padding = 0 # reduces each side by 1 during Conv

#		 self.conv1 = nn.Conv2d(3, 16, self.kernel, padding=self.padding)
#		 self.conv2 = nn.Conv2d(16, 16, self.kernel, padding=self.padding)
#		 self.pool1 = nn.MaxPool2d((2, 2), stride=2)

#		 self.conv3 = nn.Conv2d(16, 32, self.kernel, padding=self.padding)
#		 self.conv4 = nn.Conv2d(32, 32, self.kernel, padding=self.padding)
#		 self.pool2 = nn.MaxPool2d((2, 2), stride=4)

#		 self.fc1 = nn.Linear(32 * 6 * 6, 256)
#		 self.fc2 = nn.Linear(256, 2) # Theirs is 1

#	 def forward(self, x):
#		 # input per image [3, 63, 63]
#		 x = F.relu(self.conv1(x)) # [16, 61, 61]
#		 x = self.pool1(F.relu(self.conv2(x))) # [16, 59, 59] -> [16, 29, 29]
#		 x = F.dropout(x, p=0.25)

#		 x = F.relu(self.conv3(x)) # [32, 27, 27]
#		 x = self.pool2(F.relu(self.conv4(x))) # [32, 25, 25] -> [32, 6, 6]
#		 x = F.dropout(x, p=0.25)

#		 x = torch.flatten(x, 1) # flatten all dimensions except batch [32 x 6 x 6]
#		 x = F.relu(self.fc1(x)) # [256]
#		 x = F.relu(self.fc2(x)) #[2]
#		 x = F.dropout(x, p=0.5)

#		 return x

# class ZTFNet(nn.Module):
#	 # their inputs are 63 x 63 x 3
#	 def __init__(self):
#		 super().__init__()
#		 self.kernel = 5
#		 self.padding = 2

#		 self.conv1 = nn.Conv2d(3, 32, self.kernel, padding=self.padding)
#		 self.conv2 = nn.Conv2d(32, 32, self.kernel, padding=self.padding)
#		 self.pool1 = nn.MaxPool2d(2, 2)

#		 self.conv3 = nn.Conv2d(32, 64, self.kernel, padding=self.padding)
#		 self.conv4 = nn.Conv2d(64, 64, self.kernel, padding=self.padding)
#		 self.pool2 = nn.MaxPool2d(4, 4)

#		 self.fc1 = nn.Linear(64 * 7 * 7, 8)
#		 self.fc2 = nn.Linear(8, 2)

#	 def forward(self, x):

#		 # input for each image: [3, 63, 63]
#		 x = F.relu(self.conv1(x)) # [32, 63, 63]
#		 x = self.pool1(F.relu(self.conv2(x))) # [32, 31, 31]
#		 # x = F.dropout(x, p=0.5)

#		 x = F.relu(self.conv3(x)) # [64, 31, 31]
#		 x = self.pool2(F.relu(self.conv4(x))) # [64, 7, 7]
#		 # x = F.dropout(x, p=0.55)

#		 # ZTF has metadata branch, add LATER

#		 x = torch.flatten(x, 1) # flatten all dimensions except batch [64 x 7 x 7]
#		 x = F.relu(self.fc1(x)) # [8]
#		 x = F.dropout(x, p=0.2)
#		 x = F.relu(self.fc2(x)) #[2] one for each class

#		 return x
