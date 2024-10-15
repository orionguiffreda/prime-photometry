# PRIME Photometry Pipeline - Multi_master Instruction & Basic Usage

This .py file extends the functionality of _photometrus stack_, or _master.py_, to all 4 detector chips, along with including the data downloading step.  Essentially, through 1 command, one is able to download and process an entire observation in a specific band, producing stacked images for each chip.

## Using Multi_master.py

Before running this file, it's good practice to verify that your observation target has been observed in the date and filter you think it is.  You can verify by checking the ramp log at the link below:

http://www-ir.ess.sci.osaka-u.ac.jp/prime_staff/LOG/Ramp_LOG/

Simply navigate to your date and search for your field and filter.

You can run this file and see the arguments through the format:

    photometrus pipeline -h

Using the default command 'photometrus' with the keywrod 'pipeline' allows utilization of _multi_master.py_.

Here is the format:

> photometrus [-h] [-parent PARENT] [-target TARGET] [-date DATE] [-band BAND] [-chip CHIP] [-rot_val ROT_VAL] [-parallel] [-no_download] [-astromnet]

The necessary positional arguments include _-parent_, _-target_, _-date_, and _-filter_.  We will go over these now.
- _-parent_ is is a string that is the pathname to the parent directory where you want your ramp file folders, along with the folders for the stacked images and intermediate images, to be stored.

- _-target_ is a string that denotes the field you want to download.  In the log, this corresponds to the OBJNAME field.  This is the location on the grid your observation was taken.

- _-date_ is a string that denotes the date the observation was taken.  This tells the downloader which log file to look at.  The format of this argument is: _yyyymmdd_

- _-band_ is a 1 character string that denotes the filter of the observation you want to download.  In the log, this corresponds to the FILTER2 field (in the case of Z band, it is FILTER1).

Optional arguments include _-parallel_, _-no_download_, _-astromnet_, & _-chip_.
- _-chip_ is an optional argument.  If you want to only download and process a single chip or a couple chips' observations, you can specify them here.  The format is comma separated, like _-chip 1,2,3,4_

- _-rot_val_ denotes the rotation offset of the telescope during the observation.  Currently, PRIME images possess an error that prevents the pipeline from correctly reading the telescope's rotation value.  To counteract this, we have this optional arg.  The default for this arg is 48 (degrees), as that is the default telescope rotation offset.  However, if you're processing an observation, be sure to check the _ROToffset_ column in your observation request .csv.  The default should be 172800 (arcsec).  If it is something different, you must convert the value from arcsec to deg and then input that value here.

- _-parallel_ is an optional argument for running the processing of chips in parallel.  This is currently in development and should not be used.

- _-no_download_ is an optional argument to run _multi_master.py_ if you already have the data downloaded.

- _-astromnet_ is an optional argument for the astrometry improvement step, utilizing astrometry.net through _gen_astrometry.py_ in _master.py_ (examine the documentation for either script for more info).  This is much less consistent and slower, but more accurate when it works.

Now that we've gone over the arguments associated with this script, we'll go over a sample command.  Say we have observed a target in field1111 in J band only, on 06/06/2024.  To download and process that target's data, this is how this script would be used in a normal scenario.  

    photometrus pipeline -parent ../../field1111/ -target field1111 date 20240606 -band J

## How Multi_master.py Works

- To begin, the _-parent_ directory is generated (or specified if the directory already exists).
- _getdata.py_ is then utilized to download a specified observation's data for all 4 chips (or specific chips) to the parent directory, along with the accompanying log file.
- Once the data is downloaded, _master.py_ is run on each chip individually, or a single chip if specified.  Refer to the documentation on _master.py_ to see the inner workings of that processing.

## Data Products

There are many, many data products from this pipeline.  Simply enough, it is 4x the amount than a single run of _master.py_, delivering all those data products for each chip (if, of course, you ran this without the -chip optional arg).  Below is a sample parent directory using the aforementioned example command.  Below, only the intermediate data products for C1 are shown (for the sake of space).

      ├── field1111
      │   ├── C1
      │   │   ├── **C1.ramp.fits...
      │   ├── C1_astrom
      │   │   ├── **C1.ramp.new...
      │   ├── C1_FF
      │   │   ├── **C1.flat.fits...
      │   ├── C1_sub
      │   │   ├── old
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
