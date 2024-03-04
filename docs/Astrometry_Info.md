# PRIME Photometry Pipeline - Astrometry Scripts

There are two scripts in the pipeline that deal with astrometry for input images.  _gen_astrometry.py_ generates initial astrometry for ramp images using astrometry.net.  _astrometry.py_ utilizes sextractor and scamp to generate improved astrometry for the purpose of stacking processed images.  Below, usage and info for both scripts is shown.

## Gen_astrometry.py 

Before using this script, make sure you possess the right index files (as explained in the _master.py_ documentation).  In addition, if you cloned this repo, you'll have to change or remove '--backend-config' field from the functions.  

### Gen_astrometry - Usage

Running the --help flag reveals the formatting is:

> gen_astrometry.py [-h] [-list] [-dir] [f ...]
