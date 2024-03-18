# PRIME Photometry Pipeline - Stacking Instructions

_stack.py_ is this pipeline's script for stacking sets of images.  It utilizes a conda package for Swarp, which utilizes external header files (with improved astrometry) matched with processed images to create a final, median-combined, stacked image. 

## Usage

_stack.py_ has the formatting as followed:

> stack.py [-h] [-mask] [-sub SUB] [-stack STACK] [-parent PARENT] [-chip CHIP]

The arguments you will utilize (currently) are the following:

- _-sub_ is a string that is the path to where your input processed images are stored.  These will be stacked.

- _-stack_ is a string that is the path to where your final stacked image & weight image will be stored. 

### Bad Pixel Mask - *IN DEVELOPMENT*

The arguments below are currently in development, but I've put them in for an initial test.  They technically work if given the right files and paths, but the swarped end product is subpar, so it's not currently part of the pipeline.

- _-mask_ is be an optional flag, utilized when you want to use a bad pixel mask to remove bad pixels from your images.

- _-parent_ is be a string for the path to the greater parent directory, where all the ramps and other images are stored. It's the same path as _-parent_ in _master.py_.

- _-chip_ is an integer for the chip number.

What this part of _stack.py_ does is apply an already created bad pixel mask (there will be one for each detector) that identifies bad pixels in each detector.  This mask is a FITS file, but it is converted to a boolean array and then applied to each image.  This removes bad pixels by converting all the false values into nan values.  The images are then rewritten with these nan values.  These masked images are then inputted into swarp for stacking.

In THEORY, this should just remove bad parts of the detector in the final image.  However, the stacked images I've generated using this do not look very good at all, even with trying with multiple different swarp config settings.  As a result, this isn't part of the current pipeline.  I may experiment more with swarp weight-maps instead in hopes of getting a good stacked image, so this may be removed eventually.

## Output files

There are two output files this script will deliver.  One is the final stacked image, which has the following format:

> coadd.Open-#.12345678-12345678.C#.fits

Here, the # after 'Open' is the letter of the filter, ex. Open-J.  The two 8 digit sequences in the middle are the numbers corresponding to the first and last ramp image processed, showing the range of numbers stacked.  The final # after 'C' is the number of the chip.  

The other output file is the weight-map of the final image.  

> weight.Open-#.12345678-12345678.C#.fits
