# PRIME Photometry Pipeline - Multi_master Instruction & Basic Usage

This .py file extends the functionality of _master.py_ to all 4 detector chips, along with including the data downloading step.  Essentially, through 1 command, one is able to download and process an entire observation in a specific band, producing stacked images for each chip.

## Using Multi_master.py

Before running this file, it's good practice to verify that your observation target has been observed in the date and filter you think it is.  You can verify by checking the ramp log at the link below:

http://www-ir.ess.sci.osaka-u.ac.jp/prime_staff/LOG/Ramp_LOG/

Simply navigate to your date and search for your field and filter.

You can run this file and see the arguments through the format:

    python ./multi_master.py -h

Here is the format:

> python ./multi_master.py [-h] [-parallel] [-parent PARENT] [-target TARGET] [-date DATE] [-filter FILTER] [-chip CHIP]

The necessary positional arguments include _-parent_, _-target_, _-date_, and _-filter_.  We will go over these now.
- _-parent_ is is a string that is the pathname to the parent directory where you want your ramp file folders, along with the folders for the stacked images and intermediate images, to be stored.

- _-target_ is a string that denotes the field you want to download.  In the log, this corresponds to the OBJNAME field.  This is the location on the grid your observation was taken.

- _-date_ is a string that denotes the date the observation was taken.  This tells the downloader which log file to look at.  The format of this argument is: _yyyymmdd_

- _-filter_ is a 1 character string that denotes the filter of the observation you want to download.  'J' or 'H' is currently supported.  In the log, this corresponds to the FILTER2 field.

Optional arguments include _-parallel_ and _-chip_.
- _-parallel_ is an optional argument for running the processing of chips in parallel.  This is currently in development and should not be used.

- _-chip_ is an optional argument.  If you want to only process a single chip's observations, you can specify the chip here.

Now that we've gone over the arguments associated with this script, we'll go over a sample command.  Say we have observed a target in field1111 in J band only, on 06/06/2024.  We first create a parent directory '/field1111/'.  To download and process that target's data, this is how this script would be used.  

    python ./multi_master.py -parent ../../field1111/ -target field1111 date 20240606 -filter J

## How Multi_master.py Works

- To begin, _getdata.py_ is utilized to download a specified observation's data for all 4 chips to the parent directory, along with the accompanying log file.
- Once the data is downloaded, _master.py_ is run on each chip individually, or a single chip if specified.  Refer to the documentation on _master.py_ to see the inner workings of that processing.
- Once the final stacked image is created, _multi_master.py_ is currently set up to automatically correct for any weak astrometry on the edges of the image through running astrometry.net on the stacked image.  This regenerates the stacked image with new astrometry, along with creating a .wcs file for each stack.
- Finally, the new stacked images are compressed into .fz files with fpack.  If you need to uncompress them (for perhaps photometry purposes) use funpack.

## Data Products

There are many, many data products from this pipeline.  Simply enough, it is 4x the amount than a single run of _master.py_, delivering all those data products for each chip.  Below is a sample parent directory using the aforementioned example command.  Below, only the intermediate data products for C1 are shown (for the sake of space).

      ├── field1111
      │   ├── C1
      │   │   ├── **C1.ramp.fits...
      │   ├── C1_astrom
      │   │   ├── **C1.ramp.new...
      │   │   ├── **C1.ramp.wcs...
      │   ├── C1_FF
      │   │   ├── **C1.flat.fits...
      │   ├── C1_sub
      │   │   ├── **C1.sky.flat.fits...
      │   │   ├── **C1.sky.flat.cat...
      │   │   ├── **C1.sky.flat.head...
      │   ├── C2
      │   ├── C2_astrom
      │   ├── C2_FF
      │   ├── C2_sub
      │   ├── C3
      │   ├── C3_astrom
      │   ├── C3_FF
      │   ├── C3_sub
      │   ├── C4
      │   ├── C4_astrom
      │   ├── C4_FF
      │   ├── C4_sub
      │   ├── sky
      │   │   ├── sky.Open-#.12345678-12345678.C1.fits
      │   │   ├── sky.Open-#.12345678-12345678.C2.fits
      │   │   ├── sky.Open-#.12345678-12345678.C3.fits
      │   │   ├── sky.Open-#.12345678-12345678.C4.fits
      │   ├── stack
      │   │   ├── coadd.Open-#.12345678-12345678.C1.fits.fz
      │   │   ├── coadd.Open-#.12345678-12345678.C2.fits.fz
      │   │   ├── coadd.Open-#.12345678-12345678.C3.fits.fz
      │   │   ├── coadd.Open-#.12345678-12345678.C4.fits.fz
      │   │   ├── coadd.Open-#.12345678-12345678.C1.wcs
      │   │   ├── coadd.Open-#.12345678-12345678.C2.wcs
      │   │   ├── coadd.Open-#.12345678-12345678.C3.wcs
      │   │   ├── coadd.Open-#.12345678-12345678.C4.wcs
      │   │   ├── weight.Open-#.12345678-12345678.C1.fits
      │   │   ├── weight.Open-#.12345678-12345678.C2.fits
      │   │   ├── weight.Open-#.12345678-12345678.C3.fits
      │   │   ├── weight.Open-#.12345678-12345678.C4.fits
      │   ├── ramp_fit_log_2024-06-06.dat
      │   ├── ramp_fit_log_2024-06-06.clean.dat

To examine the data products, explore the _master.py_ documentation.  This script simply does the work of _master.py_ across all 4 chips.  
