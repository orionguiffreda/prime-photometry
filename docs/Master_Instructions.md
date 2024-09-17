# PRIME Photometry Pipeline - Master Instruction & Basic Usage

## Before using Master.py

Before trying to use, make sure of course everything is installed correctly.  In addition, make sure your astrometry.net index files are downloaded.

## Using Master.py

Once the pipeline is installed, make sure to add the _/photomitrus_ directory to your PYTHONPATH.  Once you've made sure to build the conda environment from the provided .yaml file, cd to _/photomitrus_ and try running the command:

	python ./master.py -h

Master.py should run with no errors and present text that details the format of the command line arguments.  Currently, the usage is as follows:

> python ./master.py [-h] [-angle] [-fpack] [-shift] [-skygen_start] [-parent PARENT] [-chip CHIP] [-filter FILTER] [-sigma SIGMA]

We will learn what all these arguments mean, starting with the positional arguments.  
- _-parent_ is a string that is the pathname to the directory where your pre-processed ramp FITS files are stored.  If the data was downloaded using _getdata.py_ then a good directory structure should look similar to this:
	```
	├── J_Band
	│   ├── C1
	│   │   ├── **C1.ramp.fits...
	│   ├── ramp_fit_log_****-**-**.dat
	│   ├── ramp_fit_log_****-**-**.clean.dat
	```
	Here we've previously downloaded only C1's data for J Band.  _-parent_ should, in this case, be the path for _J_Band_

- _-chip_ is an integer that tells the pipeline which detector is being processed.  In the previous example, we would simply use '1' for this argument.  Currently, this master file only processes 1 chip at a time, though perhaps improvements will be coming soon!

- _sigma_ is an integer that is used in the creation of the sky for sky subtraction (default = 4).  This value is how many &#963; away from the median value will be clipped.  If &#963; = 4, pixels with values >4&#963; will be clipped and replaced with the median value.  If you would like to know more, look to the individual documentation for _gen_sky.py_.  Due to the variation in image quality when comparing ramps from each chip, I usually use a sigma of 4 for C1 and C2, with a sigma of 6 for C3 and C4.

- _skygen_start_ is an optional argument for starting the pipeline on the sky generation step.  This would be used if you had already done the astrometry on all ramp files or for some reason the pipeline errored out after the astrometry step.  There are also similar _...start_ arguments for every step of the pipeline, so one can easily continue their processing progress if a certain step errors out.  These do require updating in the future however, and are outdated.

- _-filter_ is an argument that must be included with the flat-fielding step.  This is simply the letter corresponding to the filter of the images you're processing, ex. "J".

Optional arguments are important and part of the default usage of the pipeline in _multi-master.py_.  Below are the optional arguments currently in use.

- _-FF_ includes a flat-fielding step.  Currently, this extra step happens after the initial astrometry.  The ramp images (w/ astrometry) are divided by an appropriate master flat.  These master flats are stored in their own subdirectory called _mflats_, and are called by a function in _settings.py_ by their filter and chip number.  Currently, this only grabs the appropriate flats that have already been created.  These master flats were created using _gen_flat.py_.  In the future, I hope to have an addition which will read the night's log, determine if flats have been taken for the same filter as the images you're processing, then generate a new flat before processing.  This will combat any temperature adjustments (causing thermal noise) made between when I made these master flats and future observations.

- _-angle_ forgoes initial usage of astrometry.net (due to inconsistency with that software), utilizing the telescope pointing and corners of the detectors to generate a quick astrometric solution using astropy (_astrom-angle.py_).  This solution is not as accurate as astrometry.net, but it is much quicker and far more consistent, so this arg is part of normal usage. NOTE: this step relies on the fact that the rotation offset of the telescope has not been changed from the default.

- _-shift_ further improves astrometry later down the pipeline.  This arg corrects for the inaccuracy in _-angle_ (usually a correction of a translation offset) through using _astrom-shift.py_ and applies the solution to all processed files.

- _-refine_ is not currently part of the normal usage of the pipeline, but this arg will apply astrometry.net instead of the algorithm in _-shift_.  When it works, it works well! But unfortunately it is greatly inconsistent across chips and observations, so it is not currently in use.

- _-fpack_ optionally compresses the final stacked image using fpack.  This can save storage space, especially if you're trying to upload final images to google drive or somewhere else.

I acknowledge it is confusing to have several optional arguments in normal usage of the pipeline, so in future refactoring I may change things.  To combat this complexity, _multi-master.py_ simplifies _master.py_'s normal usage.  A sample argument utilizing normal usage would look like:

	python ./master.py -FF -angle -shift -parent ../../J_Band/ -chip 1 -filter J -sigma 4

## How Master.py Works

This section details how _master.py_ works.  

- To begin, it imports the makedirs function from _settings.py_, which creates directories to store the ramp images with astrometry, the generated sky, the processed (cropped and sky-subtracted) images, and the final stacked image.  
- Then the initial astrometry step is run.  The normal usage involves using _astromangle-new.py_ to generate quick astrometry.  Astrometry.net with _gen-astrometry.py_ is also an option.
- It is at this point, where if you included flat-fielding, the ramp images would be cropped accordingly and divided by the master flat specific to the filter and chip, using _flatfield.py_.
- Another command is then run that calls _gen_sky.py_, which generates a sky image based on the given sigma value.  The sky image is then subtracted from each image using _sky.py_.
- Further astrometric refinement is then run, normal usage involves _astrom-shift.py_ to quickly solve for the translation offset from the initial step.  Astrometry.net is also an option.
- Once each image has been processed, _astrometry.py_ is used to run sextractor and scamp to generate final astrometry for stacking purposes, through generating sextractor .cat and scamp .head files.
- Finally _stack.py_ is utilized, running swarp to co-add all the processed images into 1 final stacked image.  This image sometimes has bad astrometry at the very edges, so astrometry.net can be run again on this stacked image, correcting the edge astrometry.
- If one wants compressed data products, the _-fpack_ flag will fpack the final stacked image optionally.   

## Data Products

There are many outputted data products from this pipeline.  However, _master.py_ sorts them nicely in folders within the parent directory.  We can look at a sample parent directory to see how things will be stored.  Below, is a parent directory after _master.py_ has processed 1 chip's worth of data. 

	├── J_Band (Parent)
	│   ├── C1
	│   │   ├── **C1.ramp.fits...
	│   ├── C1_astrom
 	│   │   ├── **C1.ramp.new...
   	│   ├── C1_sub
    	│   │   ├── old
  	│   │   ├── **C1.sky.flat.fits...
  	│   │   ├── **C1.sky.flat.cat...
  	│   │   ├── **C1.sky.flat.head...
   	│   ├── sky
	│   │   ├── sky.Open-#.12345678-12345678.C#.fits
  	│   ├── stack
	│   │   ├── coadd.Open-#.12345678-12345678.C#.fits
 	│   │   ├── weight.Open-#.12345678-12345678.C#.fits
  	│   ├── ramp_fit_log_****-**-**.dat
	│   ├── ramp_fit_log_****-**-**.clean.dat

We've already gone over C1 (the directory where the ramps are stored).  We'll now go over the rest of the subdirectories and what is stored within.

- C1_astrom is where _gen_astrometry.py_ deposits the new ramps that now have astrometry, hence the name '.ramp.new'.  Astrometry.net also outputs the '.ramp.wcs' wcs files, though since the WCS for the new ramps are written in the header, we don't use these.
- If you used flat-fielding, there would be another directory included called C1_FF.  C1_FF is a directory which stores the flat fielded ramps.  These files have names '.flat.fits'.
- C1_sub is where the sky subtracted images are stored.  These are the 'sky.flat.fits' files.  This subdirectory is also where the _astrometry.py_ places the sextractor catalogues ('sky.flat.cat') and the scamp .head ('.sky.flat.head') files, for better astrometry when swarp is used.  The subdirectory 'old/' is created by _astrom-shift.py_ and is where the processed images w/o the improved solution are moved to (this is currently done to combat a disk issue on GoddardPC01).  
- sky is where _gen_sky.py_ outputs the sky image.  The filename is similar to the final stacked image, look at _photometry.py_ or _stack.py_'s documentation for more info on how this name is generated.
- stack is where the final co-added image for the chip is placed.  
