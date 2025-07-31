from astropy import units as u

from photomitrus.settings import PIPELINE_DEFAULT_DIR, gen_config_file_name, PHOTOMETRY_MAG_UPPER_LIMIT, PHOTOMETRY_MAG_LOWER_LIMIT

PROCESSING_DEFAULTS = dict(
    coord_file_sep='\s+', coord_file_object_field=None, coord_ra_field='RA', coord_dec_field='DEC', frame='icrs',
    unit=u.degree,
    grid_file=gen_config_file_name('obsable_all_sky_grid.csv'),
    parent=PIPELINE_DEFAULT_DIR, band='J', date=None,
    target=None, chip='1,2,3,4', rot_val=48, no_shift=False, astromnet=False,
    sky_override_path=False, removal=False, no_get_files=False, no_download=False, no_mflat=False, sigma=4, bulge=False,
    sex=False, compress=False, ramplist=None, automode=False,
    stackpath=None, survey=None, grb_ra=None, grb_dec=None, grb_coordlist=None, grb_radius='4.0', exp_query=False,
    no_plots=False, keep=False, grb_only=False, filepath=None, crop=300, sigma_photom=3, catalog=None,
    mag_low=PHOTOMETRY_MAG_LOWER_LIMIT, mag_high=PHOTOMETRY_MAG_UPPER_LIMIT, int_cal=False
)

FILE_DEFAULTS = dict(
    save_dir='.', redownload=False, overwrite=False, ftype='ramp', ip=None, user=None,
    password=None, objname=None, objtype=None, observer=None, chip=(1,2,3,4), filter1=None, filter2=None,
    n_retry=3,
    funpack_fz=False

)

ASTROM_DEFAULTS = dict(range=3,length=100,num=15,stdev=1,segstd=2)
ASTROM_DEFAULTS_NEW = dict(length=125,num=400,thresh_high=0.4,thresh_low=0.1,iters=10)
