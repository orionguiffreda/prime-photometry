import pandas as pd
import numpy as np
import subprocess
import ast
from collections import Counter
from datetime import datetime, timedelta
from pathlib import Path
import matplotlib.pyplot as plt




import astropy.units as u
import astropy.time
import astropy.coordinates
from astropy.io import fits

import assess_processed_data_util as gen_util



def select_candidates(tns_csv='data_results/tns_8_20.csv', combo_file='data_results/combo_commands_2026-09-03'):

	# read in, reshape TNS data

	df = pd.read_csv(tns_csv)
	
	df['discovery_date'] = pd.to_datetime(df['discoverydate'],format='%m/%d/%y')
	
	mydf = df[df['declination'] < 30].copy()
	mydf['discoverydate_datetime'] = pd.to_datetime(df['discoverydate'],format='%m/%d/%y')
	mydf = mydf[mydf['discoverydate_datetime'] > pd.to_datetime('01/01/24')]
	print(len(mydf), "relevant transients in TNS")
	print()
	
	threshold_deg = 0.001
	data_list = []
	
	with open(combo_file, 'r') as f:
	    for line in f:
	        line = line.replace("np.float64(", "").replace(")", "")
	        line = ast.literal_eval(line)
	        data_list.append(line)
	
	df_prime = pd.DataFrame(data_list)
	df_prime['datetime'] = pd.to_datetime(df_prime['date'], format='%Y%m%d')
	df_prime['field'] = [Path(mydict['parentdir']).parts[-2] for idx, mydict in df_prime.iterrows()]
	df_prime['grid_field'] = [mydict['field'][:-11] for idx, mydict in df_prime.iterrows()]
	
	
	print(len(df_prime), "identified target fields overlapping with these detections")
	

	
	df['dt'] = [[] for i in range(len(df))]
	df['field'] = [[] for i in range(len(df))]
	df['grid_field'] = [[] for i in range(len(df))]
	df['band'] = [[] for i in range(len(df))]
	
	df['command'] = [[] for i in range(len(df))]
	df['discovery_date'] = pd.to_datetime(df['discoverydate'],format='%m/%d/%y')
	
	for idx, prime in df_prime.iterrows(): # for potential observation in prime
	
	    # for each idx where ra and dec of transient was close to observation
	    idxs = df[(abs(df['ra'] - prime['grb_ra']) < threshold_deg) & (abs(df['declination'] - prime['grb_dec']) < threshold_deg)].index.tolist()
	    for idx in idxs:
	
	        # add dt (between prime obs and our tns detection) to the dataframe
	        dt = (df.loc[idx, 'discovery_date']-prime['datetime']).days
	        df.loc[idx, 'dt'].append(dt)
	        df.loc[idx, 'field'].append(prime['field'])
	        df.loc[idx, 'grid_field'].append(prime['grid_field'])
	        df.loc[idx, 'band'].append(prime['band'])
	        
	        prime_original = prime.copy()
	        del prime_original['datetime'] 
	        del prime_original['field']
	        del prime_original['grid_field']
	
	        df.loc[idx, 'command'].append(prime_original)
	
	
	
	
	    
	df_nonempty=df.copy()
	# df_nonempty['min_dt'] = df_nonempty['dt'].apply(lambda x: min(map(abs, x)) if x else float('inf'))
	df_nonempty['min_dt'] = df_nonempty['dt'].apply(lambda x: min(x, key=abs) if x else None)
	
	
	sorted_ids = df_nonempty.sort_values('min_dt')
	# plt.hist(sorted_ids['min_dt'])
	# plt.show()
	
	close_obs = sorted_ids[np.abs(sorted_ids['min_dt'])<21]
	print(len(close_obs), "obs within a 20 days of a TNS obs")
	
	start_date = datetime(2024, 1, 1, 0, 0, 0)
	close_obs = close_obs[close_obs['discovery_date']>start_date]
	print(len(close_obs), "obs after 2023")
	
	close_obs = gen_util.add_galactic_coords(close_obs)
	close_obs = close_obs[abs(close_obs['b'])>16]
	print(len(close_obs), "obs in sparse locations")
	
	# close_obs = close_obs[close_obs['discoverymag']<19.5]
	# print(len(close_obs), "obs brighter than 19.5")
	
	def count_refs(row):
	    if row['min_dt']>8:
	        return 0
	    
	    min_idx = np.argmin(list(map(abs, row['dt'])))
	    sci_key = (row['grid_field'][min_idx], row['band'][min_idx])
	    all_keys = list(zip(row['grid_field'], row['band']))
	
	    return sum(k == sci_key for i, k in enumerate(all_keys) if i != min_idx)
	
	close_obs['ref_match_count'] = close_obs.apply(count_refs, axis=1)
	del close_obs['grid_field']
	
	
	solo_field = close_obs[close_obs['ref_match_count'] == 0]
	multi_field = close_obs[close_obs['ref_match_count'] > 0]
	
	print(len(multi_field), "obs with a ref of matching field and band")
	print(f"({len(solo_field)} solo obs)")
	
	
	top_hits_df = close_obs 
	
	
	list(multi_field["full_name"][100:])
	recent_obs = pd.read_csv("data_results/single_detection_8_20.csv")
	old_triplet = pd.read_csv("data_results/full_data.csv")
	
	multi_field_new = multi_field[multi_field["full_name"].isin(recent_obs["full_name"])]
	multi_field_old = multi_field[multi_field["full_name"].isin(old_triplet["full_name"])]
	print(len(multi_field_new), "new obs, ", len(multi_field_old), "old obs,",len(multi_field), "total")
	
	
	multi_field_target = multi_field #pd.concat([multi_field_new, multi_field_old])
	
	
	
	by_obs = []
	
	multi_field_target['truncated_field'] = [[field[:-11] for field in row['field']] for i, row in multi_field_target.iterrows()]
	
	for j, row in multi_field_target.iterrows():
	    # if row['min_dt']>8:
	        # continue
	    # print(row['dt'])
	        
	    sci_idx = np.argmin(list(map(abs, row['dt'])))
	    sci_field = row['truncated_field'][sci_idx]
	    sci_band = row['band'][sci_idx]
	
	    # REF = obs with matching field and band that has the smallest dt
	    
	    
	    matching_ref_idx = []
	    ref_idxs = np.argsort(list(map(abs, row['dt'])))
	    ref_idx = -1
	    for i in ref_idxs:
	        sci = i == sci_idx
	
	        # change this to include only one ref 
	        ref = row['truncated_field'][i] == sci_field and row['band'][i] == sci_band and i!=sci_idx 
	        if sci or ref:
	            try:
	                new_row = row.copy()
	                new_row['field'] = row['field'][i]
	                new_row['dt'] = row['dt'][i]
	                new_row['command'] = row['command'][i]
	                by_obs.append(new_row)
	                if ref: 
	                    ref_idx = i
	            except:
	                print("ERROR")
	            
	# finding sci even if no ref         
	by_obs = pd.DataFrame(by_obs)
	
	new_file_path = combo_file + '-top-candidates'
	
	with open(new_file_path, 'w') as f_out:
	    for command_list in by_obs['command']:
	            command_dict = dict(command_list)
	            
	            # print(mydict)
	            f_out.write(str(command_dict) + '\n')
	
	print("wrote", len(by_obs), "good transients to file", new_file_path)
	
	return new_file_path
	            


def reduce_combo(combo_file):

	base_path = "/home/alex/PycharmProjects/prime-photometry-fiona/automated_detection"
	combo_out_file = f"{base_path}/data_results/new_combo_file"
	
	result = subprocess.run(
	    [
	        "photometrus", "archive",
	        "-coord_ra_field", "ra",
	        "-coord_dec_field", "declination",
	        "-coord_file_object_field", "full_name",
	        "-coord_file_sep", ",",
	        "-coord_file", "abc.csv",
			"-combo_out_file", combo_out_file,
	        "-combo_dict_file", combo_file,

	    ],
	    # check=True,
	    # text=True,
	    # capture_output=True,
	)
	
	print(result.stdout)
	print(result.stderr)
	return combo_out_file

def identify_observations(coord_file):

	base_path = "/home/alex/PycharmProjects/prime-photometry-fiona/automated_detection"
	combo_out_file = f"{base_path}/data_results/new_combo_file"
	
	result = subprocess.run(
	    [
	        "photometrus", "archive",
	        "-coord_ra_field", "ra",
	        "-coord_dec_field", "declination",
	        "-coord_file_object_field", "full_name",
	        "-coord_file_sep", ",",
	        "-coord_file",  f"{base_path}/{coord_file}",
	        "-combo_out_file",  combo_out_file,
			"-no_reduce"
			
	    ],
	    # check=True,
	    # text=True,
	    # capture_output=True,
	)
	
	print(result.stdout)
	print(result.stderr)
	return combo_out_file



	

