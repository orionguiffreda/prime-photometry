# PRIME Photometry Pipeline - Astrometry Scripts

There are two scripts in the pipeline that deal with astrometry for input images.  _gen_astrometry.py_ generates initial astrometry for ramp images using astrometry.net.  _astrometry.py_ utilizes sextractor and scamp to generate improved astrometry for the purpose of stacking processed images.  Below, usage and info for both scripts is shown.

## Gen_astrometry.py 

Before using this script, make sure you possess the right index files (as explained in the _master.py_ documentation).  In addition, if you cloned this repo, you'll have to change or remove '--backend-config' field from the functions.  The directory in this field should be the same as where 'astrometry.cfg' is located.  

### Gen_astrometry - Usage and output files

Running the --help flag reveals the formatting is:

> gen_astrometry.py [-h] [-list] [-output OUTPUT] [-input INPUT]

-_-list_ is an optional flag, for use only if you're using a list of image paths as input.  Currently, this really isn't used.  However, if the ramp and raw image directories are mounted to PC-1, then this will likely see use. 

-_output_ is the output path for the images w/ new astrometry.  

-_input_ is the input field.  This is most likely a directory with ramp FITS images (currently).    

The output files are new FITS files w/ WCS in the header, with the format '.ramp.new'.  It also outputs separate wcs files as '.ramp.wcs.'

## Astrometry.py

This script is located in a different folder (_/astrom/_) than the previous script.
