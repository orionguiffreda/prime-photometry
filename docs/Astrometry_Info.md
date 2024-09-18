# PRIME Photometry Pipeline - Astrometry Scripts

There are four scripts in the pipeline that deal with astrometry for input images.  

For initial astrometry, one can choose either _astrom_angle.py_ or _gen_astrometry.py_.  The former uses the pointing of the telescope and default rotation offset to generate quicker astrometric solutions (this is the normal pipeline choice), while the latter generates initial astrometry for ramp images using astrometry.net.  

If _astrom_angle.py_ is used, then a further improvement is run using _astrom_shift.py_, which corrects for translation error in the initial step.

_astrometry.py_ is the final refinement, utilizing sextractor and scamp to generate improved astrometry for the purpose of stacking processed images.  Below, usage and info for both scripts is shown.

## Gen_astrometry.py 

This script utilizes a conda package of astrometry.net on input images, generating astrometry and creating new images with this astrometry in the header.  Before using this script, make sure you possess the right index files (as explained in the _master.py_ documentation).  In addition, if you cloned this repo, you'll have to change or remove '--backend-config' field from the functions.  The directory in this field should be the same as where 'astrometry.cfg' is located.  

### Gen_astrometry - Usage and output files

Running the --help flag reveals the formatting is:

> gen_astrometry.py [-h] [-list] [-hard] [-output OUTPUT] [-input INPUT]

- _-list_ is an optional flag, for use only if you're using a list of image paths as input.  Currently, this really isn't used.  However, if the ramp and raw image directories are mounted to PC-1, then this will likely see use.

- _-hard_ is an optional flag that uses a more specific command for astrometry.net, in an attempt to get a faster solve.  It includes variables such as pixel scale range, center RA and DEC, and search radius.  This is automatically picked for use in the pipeline.

- _output_ is the output path for the images w/ new astrometry.  

- _input_ is the input field.  This is most likely a directory with ramp FITS images (currently).    

The output files are new FITS files w/ WCS in the header, with the format '.ramp.new'.  It also outputs separate wcs files as '.ramp.wcs.'

## Astrom_shift.py

This script uses sextractor (w/ psfex) and astroquery to compare the source positions in a known catalog to source positions from a given image, correcting for a translation offset between them.  

For the sake of quickness, this comparison is done within a box of 9-12' at the center of the image.  Both the catalog & image sources are sorted by magnitude, then euclidean distances are calculated between each of the first specified # of sources in each group are calculated.  If distances agree w/in a specified # of pixels and are less than a specified # of pixels in length, they are a considered pair.  The X and Y pixel distances (distance between the catalog and image source) are examined for these considered pairs.  A median is taken for all the considered X and Y distances respectively, and any pair of sources that falls greater than a specified # of stdevs away from this median are pruned out.  Finally, a median is taken across the X and Y distances of the remaining pairs, and is used for the final X shift value and Y shift value.  These values are then written into the _CRPIX_ fields of the header of the given image.  

In the pipeline, this solution is applied to all images in a given observation.

### Astrom_shift - Usage and output files

Running the --help flag reveals the formatting is: 

> astrom_shift.py [-h] [-pipeline] [-remove] [-segment] [-dir DIR] [-imagename IMAGENAME] [-filter FILTER] [-range RANGE] [-length LENGTH] [-num NUM] [-stdev STDEV] [-segstd SEGSTD]

Let's start with the positional arguments:

- _-dir_ is the directory the image (or images) are located in.

- _-imagename_ is the filename of the image this script will be run on.

- _-filter_ is the filter used during the observation.

The optional arguments are:

- _-pipeline_ applies the solution found to all images in the directory and generates new images w/ the same name, moving the old images to a subdirectory '/old/' (this is done to combat a disk error on the computer).  This, obviously, is the default on the pipline itself.

- _-remove_ removes intermediate catalogs created by sextractor and psfex, this is also default on the pipeline.

- _-segment_ increases the accuracy of this script through splitting the box at the center of the image into 4 quadrants.  A final shift is generated for each quadrant, where a quadrant's final values will be pruned out if they do not agree w/in a specified # of stdevs (_-segstd_) from the median value.  This introduces redundancy, where if 1 quadrant is for some reason incorrect, it won't affect the final value.

The arguments below have default values and don't need to be specified in the command unless the user desires it.

- _-range_ is the acceptable range for calculated euclidean distances to agree, default = 3 pix.

- _-length_ is the maximum length of the euclidean distance that can be considered, default = 100 pix.

- _-num_ is the number of sources from both the catalog and image to have euclidean distances calculated, default = 10.

- _-stdev_ is the number of stdevs away from the median, above which considered pairs are pruned out, default = 1.

- _-segstd_ is used with the _-segment_ flag, and specifies the # of stdevs away from the median above which final values are pruned.  The median is taken across all 4 quadrants' final X and Y values.  Default = 2.

Here is example usage of this command, in the setting of normal pipeline usage (in the pipeline, the first image is automatically selected from the image directory):

    python ./astrom/astrom_shift.py -remove -segment -dir /obs_dir/C1_sub/ -imagename /first_image.fits -filter J


## Astrometry.py

This script utilizes sextractor, generating source catalogues in FITS_LDAC format.  These catalogues have specific parameters which are then read and used by scamp, which in turn generates FITS header keywords in external header files, containing updated astrometric and photometric information.  The purpose of these external headers are to be utilized with swarp in image stacking, a later step.  _astrometry.py_ is located in a different folder (_/astrom/_) than the previous script, as it deals with a later step in the pipeline.  Keep this in mind when calling it from command line.

### Astrometry - Usage and output files

Running the --help flag reveals the formatting is:

> astrometry.py [-h] [-sex] [-scamp] [-all] [-path PATH]

- _-sex_ is an optional flag, use this if you want to only use sextractor on input images.

- _-scamp_ is an optional flag, use this if you want to only use scamp.  This step requires .cat sextractor catalogues to work.

- _-all_ is an optional flag, use this if you want to use sextractor, then scamp right after.  This will be the one you use most often, and is the flag used in _master.py_.

- _-path_ is the directory where this script will run.

There are two types of output files.  '.cat' files are the sextractor catalogues, while '.head' files are the external headers.
