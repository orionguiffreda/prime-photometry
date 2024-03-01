# PRIME Photometry Pipeline - Master Instruction & Basic Usage

## Before using Master.py

Before trying to use, make sure of course everything is installed correctly.  In addition, make sure your astrometry.net index files are downloaded.  See these links for details about index files:

https://astrometry.net/doc/readme.html

http://data.astrometry.net/

The 4200 series is useful, as I first downloaded those and they've worked (mostly) just fine, being built off of 2MASS.  I recommend downloading the 4100 and some of the 5200 series, though I haven't yet tested if these are better for PRIME than the 4100 and 5200.

## Using Master.py

Once the pipeline is installed, make sure to add the _/photomitrus_ directory to your PYTHONPATH.  Once you've made sure to build the conda environment from the provided .yaml file, cd to _/photomitrus_ and try running the command:

	python ./master.py -h

Master.py should run with no errors and present text that details the format of the command line arguments.  Currently, the usage is as follows:

> python ./master.py [-h] [-skygen_start] [-parent PARENT] [-chip CHIP] [-sigma SIGMA]

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

- _sigma_ is an integer that is used in the creation of the sky for sky subtraction (default = 4).  This value is how many &#963; away from the median value will be clipped.  If &#963; = 4, pixels with values >4&#963; will be clipped and replaced with the median value.  If you would like to know more, look to the individual documentation for _gen_sky.py_

- _skygen_start_ is an optional argument for starting the pipeline on the sky generation step.  This would be used if you had already done the astrometry on all ramp files or for some reason the pipeline errored out after the astrometry step.  In the future, I aim to be able to start the pipeline on any individual step, in case things go wrong halfway through.

Now that we've gone over each argument, here is a sample command utilizing _master.py_:

	python ./master.py -parent ../../J_Band/ -chip 1 -sigma 4

This would process C1's data for J band.

## How Master.py Works

This section details how _master.py_ works.  

- To begin, it imports the makedirs function from _settings.py_, which creates directories to store the ramp images with astrometry, the generated sky, the processed (cropped and sky-subtracted) images, and the final stacked image.  
- Then, using subprocess, it runs a command which calls _gen_astrometry.py_, running astrometry.net to generate astrometry for all the ramp images in the chip's directory.
- Another command is then run that calls _gen_sky.py_, which generates a sky image based on the given sigma value.  The sky image is then subtracted from each image using _sky.py_.
- Once each image has been processed, _astrometry.py_ is used to run sextractor and scamp to generate improved astrometry for stacking purposes.
- Finally _stack.py_ is utilized, running swarp to median-combine all the processed images into 1 final stacked image!   

