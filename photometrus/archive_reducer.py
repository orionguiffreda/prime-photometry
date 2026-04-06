import datetime
import os
import traceback

import argparse
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
import ast
# from getfiles import get_log_file
from photometrus.multi_combo import combo
from photometrus.getfiles import get_log_file
from photometrus.settings import PIPELINE_DEFAULT_DIR, gen_config_file_name

# PIPELINE_DEFAULT_DIR = '/mnt/photometry'
# _grid_df = pd.read_csv('obsable_all_sky_grid.csv')
_defaults = dict(
	coord_file_sep='\s+', coord_file_object_field='ObjectName', coord_ra_field='RA', coord_dec_field='DEC', frame='icrs',
	unit=u.degree,
	grid_file=gen_config_file_name('obsable_all_sky_grid.csv'), parent=PIPELINE_DEFAULT_DIR,
	rot_val=48, no_shift=False, astromnet=False,
	sky_override_path=False, removal=False, no_get_files=False, no_download=False, no_mflat=False,
	survey=None
)
_defaults['grid_df'] = pd.read_csv(_defaults['grid_file'])
chip_dict = {
	'11': 3,
	'10': 1,
	'01': 4,
	'00': 2,
}

def get_dates(start_date=datetime.date(day=11, month=10, year=2022), end_date=datetime.date.today()):
	num_days = (end_date-start_date).days
	date_list = [(end_date - datetime.timedelta(days=x)).strftime('%Y-%m-%d') for x in range(num_days)]
	return tuple(date_list)


def get_all_date_logs(date_list=get_dates(), return_missing=False):
	logs = []
	nologs = []
	for date in date_list:
		try:
			log = get_log_file(date)
			log['date'] = date
			logs.append(log)
		except FileNotFoundError:
			print(date, 'not found')
			nologs.append(date)
	if return_missing:
		return pd.concat(logs, ignore_index=True), nologs
	else:
		return pd.concat(logs, ignore_index=True)


def find_objname_commands(objname, logs):
	obj_log = logs[logs.OBJNAME == objname]  # .date.unique()
	commands = ('photometrus full -target ' + obj_log.OBJNAME + ' -date ' + obj_log.date + ' -band ' + obj_log.FILTER1
				+ obj_log.FILTER2)
	commands = [cmd.replace('Open', '') for cmd in commands]
	obj_log['command'] = commands
	return obj_log.drop_duplicates(subset=['command'])


def get_chip_df(df):
	ra_sign = df['ra_offsets'] >= 0
	dec_sign = df['dec_offsets'] >= 0
	signs = ra_sign.astype(int).astype(str) + dec_sign.astype(int).astype(str)
	chips = signs.map(chip_dict)
	return chips


def calculate_distance_all(target, grid=_defaults['grid_df']):
	grid_coords = SkyCoord(grid['RA'], grid['DEC'], unit=('hourangle', 'degree'))
	grid['distance'] = target.separation(grid_coords).arcmin
	min_dist = np.min(grid['distance'])
	index = np.where(grid['distance'] == min_dist)
	dra, ddec = target.spherical_offsets_to(grid_coords)
	grid['ra_offsets'], grid['dec_offsets'] = dra.arcmin, ddec.arcmin
	grid['ra_degrees'] = grid_coords.icrs.ra.deg
	grid['dec_degrees'] = grid_coords.icrs.dec.deg
	grid['chip'] = get_chip_df(grid)
	return index, min_dist, grid.sort_values('distance')


def obtain_point_source_grid_df(
		dataframe, point_source_radius=3/60, dither_radius=-1, min_ra_dec_offset=3.5, max_ra_dec_offset=37.5
):
	new_df = dataframe.copy()
	print(new_df[['distance', 'ra_offsets', 'dec_offsets', 'chip']].head(10))
	# print(new_df['ra_offsets'] > settings.MIN_RA_DEC_OFFSET_ARCMIN)
	ra_abs = new_df['ra_offsets'].abs()
	dec_abs = new_df['dec_offsets'].abs()
	min_offset = min_ra_dec_offset + point_source_radius + dither_radius
	max_offset = max_ra_dec_offset - (point_source_radius + dither_radius)
	print('min_offset', min_offset)
	print('max_offset', max_offset)
	new_df = new_df.loc[(ra_abs > min_offset) & (ra_abs < max_offset)]
	new_df = new_df.loc[(dec_abs > min_offset) & (dec_abs < max_offset)]
	print(new_df)
	current_rows = new_df.shape[0]
	if current_rows == 0:
		raise ValueError('no on grid location for this object')
	return new_df


def get_grid_locations(coords, grid_df=_defaults['grid_df']):
	"""
	return a list of grids matching each coordinate
	"""
	grids = []
	# print("RA DEC:",grid_df['RA'], grid_df['DEC'])
	
	coord_series = SkyCoord(
		ra=grid_df['RA'].values,
		dec=grid_df['DEC'].values,
		unit=('hourangle', 'degree')
		)	
	for coord in coords:
		# print(grid_df.keys())
		# print(type(grid_df['RA'][0]),type(coord.icrs.ra))

		delta_ra = (coord_series.icrs.ra - coord.icrs.ra).wrap_at(180 * u.deg).deg
		delta_ra_corrected = np.abs(delta_ra * np.cos(np.pi / 180 * coord.icrs.dec.deg))
		delta_dec = np.abs((coord_series.icrs.dec - coord.icrs.dec).deg)
		cropped_df = grid_df[(delta_ra_corrected < 2) & (delta_dec < 2)]
		
		index, min_dist, grid = calculate_distance_all(coord, cropped_df.reset_index(drop=True))
		
		try:
			ongrid_df = obtain_point_source_grid_df(grid)
			grids.append(ongrid_df)
		except ValueError:
			grids.append(None)
	return grids


def get_coords(
		coordinates_file, coord_file_sep=_defaults['coord_file_sep'], coord_ra_field=_defaults['coord_ra_field'],
		coord_dec_field=_defaults['coord_dec_field'], frame=_defaults['frame'], unit=_defaults['unit']
):
	coord_df = pd.read_csv(coordinates_file, sep=coord_file_sep)
	coords = SkyCoord(coord_df[coord_ra_field], coord_df[coord_dec_field], unit=unit, frame=frame)
	coord_df['coords'] = coords
	return coord_df


# def combo(target, date, band, chip=None, parentdir=False, rot_val=48, no_shift=False, astromnet=False,
#		   sky_override_path=False, removal=False, no_get_files=False, no_download=False, no_mflat=False, survey=None,
#		   grb_ra=None, grb_dec=None, grb_coordlist=None, grb_radius=None, auto_mode=False, input_ramp_lists=None
#		   ):
#	 pass

def execute_combo_dict_file(combo_dict_file, no_reduce):
	data_list = []
	with open(combo_dict_file, 'r') as f:
		for line in f:
			line = line.replace("np.float64(", "").replace(")", "")
			line = ast.literal_eval(line)
			data_list.append(line)

			print("combo dict:", line, type(line))
			try:	   
				if not no_reduce:
					combo(**line)
				else:
					print('no_reduce')
			except Exception:
				tb = traceback.format_exc()
				print(line)
				print(tb)

def get_GRB_follow_up(logs):

	mask = (~logs.OBJNAME.str.lower().str.startswith("field", na=False)) & (logs.OBSERVER == 'NASA')
	obj_log = logs[mask].copy()
	commands = ('photometrus full -target ' + obj_log.OBJNAME + ' -date ' + obj_log.date + ' -band ' + obj_log.FILTER1
				+ obj_log.FILTER2)
	commands = [cmd.replace('Open', '') for cmd in commands]
	obj_log['command'] = commands
	return obj_log.drop_duplicates(subset=['command'])

# def make_combo_commands(object_dict, chip=None):

def archive_reducer(
	coord_file, coord_file_sep=_defaults['coord_file_sep'],
	coord_file_object_field=_defaults['coord_file_object_field'],
	coord_ra_field=_defaults['coord_ra_field'],
	coord_dec_field=_defaults['coord_dec_field'],
	frame=_defaults['frame'], unit=_defaults['unit'], grid_file=_defaults['grid_file'], parent=_defaults['parent'], rot_val=_defaults['rot_val'], no_shift=_defaults['no_shift'], astromnet=_defaults['astromnet'], sky_override_path=_defaults['sky_override_path'], removal=_defaults['removal'], no_get_files=_defaults['no_get_files'], no_download=_defaults['no_download'], no_mflat=_defaults['no_mflat'], survey=_defaults['survey'], no_reduce=False, combo_dict_file=None
):
	if combo_dict_file is not None:
		execute_combo_dict_file(combo_dict_file, no_reduce)
		return
 
	archive_command_filename = 'combo_commands_{}'.format(datetime.datetime.utcnow().isoformat())
	archive_command_filename = os.path.join(os.getcwd(), 'combo_commands', archive_command_filename)
	coord_df = get_coords(
		coord_file, coord_file_sep=coord_file_sep, coord_ra_field=coord_ra_field, coord_dec_field=coord_dec_field,
		frame=frame, unit=unit
		)

	
	archive = get_all_date_logs()
	archive_filename = os.path.join(parent, 'prime_log_archive_{}'.format(datetime.date.today().strftime('%Y-%m-%d')))
	archive.to_csv(os.path.join(archive_filename))
	print('output file to: ', archive_filename)

	use_all_sky_grid = True
	grids = None
	if use_all_sky_grid:
		grid_df = pd.read_csv(grid_file)
		grids = get_grid_locations(coord_df.coords, grid_df)

	"""
	Crossmatch inputted coordinates with previously observed locations (typically: objects on grid) from archive
	
	coord_df: dataframe storing coordinate information inputted by user
	grids: list of grid tiles that go with each coordinate in coord_df
	archive: dataframe of archival data about PRIME observations
	"""
	grid_search=True
	GRB_follow_up=False
		
	if GRB_follow_up:
		
		for i, coord_dict in coord_df.iterrows():
				if i>0:
					break

				object_df = get_GRB_follow_up(archive) # pass smth else here
				object_dicts = object_df.to_dict(orient='records')
			
				print("coord dict",coord_dict)		
				print("obj dict 1", object_dicts[0])
				for object_dict in object_dicts:
					bandpass = (object_dict['FILTER1']+object_dict['FILTER2']).replace('Open', '')
						
					combo_dict = dict(
						target=object_dict['OBJNAME'], date=object_dict['date'].replace('-', ''),
						band=bandpass,
						chip=str(object_dict['CHIP']), # take from the archive when its not from grid
						parentdir=os.path.join(
							parent, coord_dict[coord_file_object_field], # name these what they are called in archive?
							'{}-{}'.format(object_dict['OBJNAME'], object_dict['date']), bandpass
						),
						grb_ra=coord_dict['coords'].icrs.ra.deg, grb_dec=coord_dict['coords'].icrs.dec.deg, # make these object dict
						rot_val=rot_val, no_shift=no_shift, astromnet=astromnet,
						sky_override_path=sky_override_path, removal=removal, no_get_files=no_get_files,
						no_download=no_download, no_mflat=no_mflat, survey=survey,
					)
					print('combo command dict:')
					print(combo_dict)
					with open(archive_command_filename, 'a') as f:
						f.write('{}\n'.format(combo_dict))
					try:
						if not no_reduce:
							combo(**combo_dict)
						else:
							print('no_reduce')
					except Exception:
						tb = traceback.format_exc()
						print(combo_dict)
						print(tb)


	
	if grid_search:
		for i, coord_dict in coord_df.iterrows():
			print("coord dict",coord_dict)
			# objname = coord_dict.get('full_name', None)
			# if objname is None:
			# 	continue
			coord_dict = coord_dict.to_dict()
			if grids[i] is None:
				continue
			
			grid_tiles = grids[i].to_dict(orient='records')
			print(grid_tiles)
			for grid_tile in grid_tiles:
				print("tile:",grid_tile)
				object_df = find_objname_commands(grid_tile['ObjectName'], archive)
		
				# object_df_coords = find_objname_commands(objname, archive) # from coords not from grid?


				
				object_dicts = object_df.to_dict(orient='records')	
				# print("OBJ DICT",object_dict)
				
				for object_dict in object_dicts:
					# is_nasa_non_grid_field = (not object_dict['OBJNAME'].lower().startswith("field")) & (object_dict['OBSERVER'] == 'NASA')
		
					# if not is_nasa_non_grid_field: # if not (is_nasa_non_field or is_grid):
					# 	continue
					# access additional metadata from grids?
		
					bandpass = (object_dict['FILTER1']+object_dict['FILTER2']).replace('Open', '')
						
					combo_dict = dict(
						target=object_dict['OBJNAME'], date=object_dict['date'].replace('-', ''),
						band=bandpass,
						chip=str(grid_tile['chip']), # should chip come from grid_tile, not archive? 
						parentdir=os.path.join(
							parent, coord_dict[coord_file_object_field],
							'{}-{}'.format(object_dict['OBJNAME'], object_dict['date']), bandpass
						),
						grb_ra=coord_dict['coords'].icrs.ra.deg, grb_dec=coord_dict['coords'].icrs.dec.deg,
						rot_val=rot_val, no_shift=no_shift, astromnet=astromnet,
						sky_override_path=sky_override_path, removal=removal, no_get_files=no_get_files,
						no_download=no_download, no_mflat=no_mflat, survey=survey,
					)
					print('combo command dict:')
					print(combo_dict)
					with open(archive_command_filename, 'a') as f:
						f.write('{}\n'.format(combo_dict))
					try:
						if not no_reduce:
							combo(**combo_dict)
						else:
							print('no_reduce')
					except Exception:
						tb = traceback.format_exc()
						print(combo_dict)
						print(tb)
			

	
def test_main():
	coord_file = '10_G4Jy_sources_in_fields_completed_by_PRIME_J_10.6Hmag16.5_K.txt'
	coord_ra_field = 'centroid_RAJ2000_1'
	coord_dec_field = 'centroid_DEJ2000_1'
	parent_dir = os.path.join(PIPELINE_DEFAULT_DIR, 'G4Jy')
	archive_reducer(coord_file, coord_ra_field=coord_ra_field, coord_dec_field=coord_dec_field, parent=parent_dir)


def main():
	parser = argparse.ArgumentParser(description='Use to process all previous observations for a given set of coordinates.')
	parser.add_argument('-coord_file', type=str, help='[str] Required. Path to file containing coordinates and target names in a character separated format with a header row.')
	parser.add_argument('-coord_file_sep', type=str, default=_defaults['coord_file_sep'], help='[str] Optional. Separator between columns. Default is any whitespace.')
	parser.add_argument(
		'-coord_file_object_field', type=str,
		default=_defaults['coord_file_object_field'], help='[str] Optional. Header for Object Name. Default is {}'.format(
			_defaults['coord_file_object_field'])
	)
	parser.add_argument(
		'-coord_ra_field', type=str,
		default=_defaults['coord_ra_field'], help='[str] Optional. Header for Coordinate RA. Default is {}'.format(
			_defaults['coord_ra_field'])
	)
	parser.add_argument(
		'-coord_dec_field', type=str,
		default=_defaults['coord_dec_field'], help='[str] Optional. Header for Coordinate Dec. Default is {}'.format(
			_defaults['coord_dec_field'])
	)
	parser.add_argument(
		'-frame', type=str,
		default=_defaults['frame'], help='[str] Optional. Coordinate frame. Default is {}'.format(
			_defaults['frame'])
	)
	parser.add_argument(
		'-unit', type=str,
		default=_defaults['unit'], help='[str] Optional. Coordinate unit. Default is {}'.format(
			_defaults['unit'])
	)
	parser.add_argument(
		'-grid_file', type=str,
		default=_defaults['grid_file'], help='[str] Optional. Grid file used for archival search. Default is {}'.format(
			_defaults['grid_file'])
	)
	parser.add_argument('-parent', type=str, help='[str] *NOW OPTIONAL* specify parent directory to '
												  'store all data products, otherwise it will automatically generate @ '
												  'the default directory w/ the format "/target_date/band/"', default=_defaults['parent'])
	parser.add_argument('-band', type=str, help='[str] Optional, filter, ex. "J". Default is to reduce all bands.', default=None)  # TODO: implement this filter
	parser.add_argument('-no_get_files', action='store_true', help='optional flag, use if you want to download '
																	 'through old method (scp), new method passes file '
																	 'paths (new saves space and time)')
	parser.add_argument('-no_download', action='store_true', help='optional flag, use if you *ALREADY* have the data'
																  'downloaded, *NOT* to use new file path method')
	parser.add_argument('-no_shift', action='store_true', help='optional flag, DO NOT use astrometric shift'
											' script in place of astrom.net, will not use either (shift is default)')
	parser.add_argument('-astromnet', action='store_true', help='optional flag, use astrom.net to reinforce astrometry')
	parser.add_argument('-removal', action='store_true',
						help='optional flag, used to remove intermediate subdirectories and data, leaving only the '
							 'stacks & skies; intended for space saving in large nights of observation')
	parser.add_argument('-rot_val', type=float, help='[float] optional, put in your rot angle in deg,'
													 ' if you had a non-default rotation angle in your obs'
													 ' (default = 48 deg or 172800")', default=_defaults['rot_val'])
	parser.add_argument('-sky_override', type=str, help='[str], Optional path to specify already generated '
														'sky to use in sky sub, skipping sky gen. Input full file path.',
						default=_defaults['sky_override_path'])
	parser.add_argument('-no_mflat', action='store_true', help='optional flag, use if you *DO NOT* want to'
															   ' automatically generate mflats for this night if none'
															   ' exist')
	parser.add_argument('-survey', type=str, help='Specify specific survey to query for photometry (default'
												  ' picks for you), see photometrus single_photometry -h for list of '
												  'available surveys')
	parser.add_argument('-no_reduce', action='store_true', help='optional flag, use if you *DO NOT* want to'
															   ' actually reduce the data, and just want the list of '
															   ' coommands')
	parser.add_argument('-combo_dict_file', type=str, help='[str] optional flag, use combo dict file instead of coord_file') 

	args, unknown = parser.parse_known_args()
	archive_reducer(
		coord_file=args.coord_file, coord_file_sep=args.coord_file_sep, coord_file_object_field=args.coord_file_object_field, coord_ra_field=args.coord_ra_field, coord_dec_field=args.coord_dec_field,
		frame=args.frame, unit=args.unit, grid_file=args.grid_file, parent=args.parent, rot_val=args.rot_val, no_shift=args.no_shift, astromnet=args.astromnet, sky_override_path=args.sky_override,
		removal=args.removal, no_get_files=args.no_get_files, no_download=args.no_download, no_mflat=args.no_mflat, survey=args.survey, no_reduce=args.no_reduce, combo_dict_file=args.combo_dict_file
	)
if __name__ == '__main__':
	main()
