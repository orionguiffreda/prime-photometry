# PRIME Photometry Pipeline - Astrometry Scripts

There are two scripts in the pipeline that deal with astrometry for input images.  _gen_astrometry.py_ generates initial astrometry for ramp images using astrometry.net.  _astrometry.py_ utilizes sextractor and scamp to generate improved astrometry for the purpose of stacking processed images.  Below, usage and info for both scripts is shown.

## Gen_astrometry.py 

This script utilizes a conda package of astrometry.net on input images, generating astrometry and creating new images with this astrometry in the header.  Before using this script, make sure you possess the right index files (as explained in the _master.py_ documentation).  In addition, if you cloned this repo, you'll have to change or remove '--backend-config' field from the functions.  The directory in this field should be the same as where 'astrometry.cfg' is located.  

### Gen_astrometry - Usage and output files

Running the --help flag reveals the formatting is:

> gen_astrometry.py [-h] [-list] [-output OUTPUT] [-input INPUT]

-_-list_ is an optional flag, for use only if you're using a list of image paths as input.  Currently, this really isn't used.  However, if the ramp and raw image directories are mounted to PC-1, then this will likely see use. 

-_output_ is the output path for the images w/ new astrometry.  

-_input_ is the input field.  This is most likely a directory with ramp FITS images (currently).    

The output files are new FITS files w/ WCS in the header, with the format '.ramp.new'.  It also outputs separate wcs files as '.ramp.wcs.'

## Astrometry.py

This script utilizes sextractor, generating source catalogues in FITS_LDAC format.  These catalogues have specific parameters which are then read and used by scamp, which in turn generates FITS header keywords in external header files, containing updated astrometric and photometric information.  The purpose of these external headers are to be utilized with swarp in image stacking, a later step.  _astrometry.py_ is located in a different folder (_/astrom/_) than the previous script, as it deals with a later step in the pipeline.  Keep this in mind when calling it from command line.

### Astrometry - Usage and output files

Running the --help flag reveals the formatting is:

> astrometry.py [-h] [-sex] [-scamp] [-all] [-path PATH]

-_-sex_ is an optional flag, use this if you want to only use sextractor on input images.

-_-scamp_ is an optional flag, use this if you want to only use scamp.  This step requires .cat sextractor catalogues to work.

-_-all_ is an optional flag, use this if you want to use sextractor, then scamp right after.  This will be the one you use most often, and is the flag used in _master.py_.  

There are two types of output files.  '.cat' files are the sextractor catalogues, while '.head' files are the external headers.
