# PRIME Photometry Pipeline - Photometry Instruction and Basic Usage

## Before using Photometry.py

Before using this script, make sure you already possess a co-added final image, likely created through usage of _master.py_.  Currently, _master.py_ and _photometry.py_ are separate, requiring 1 command to process into a final image and another command to run photometry.  I could (and plan to) combine them into 1 big script that does everything in 1 go.  However, while combining the basic functionality of _photometry.py_ with _master.py_ would be relatively easy, due to the many current (and future) optional arguments for _photometry.py_ I don't want to have a very cluttered single command that would be prone to misinputting, leading to many minutes of lost time.  So, for now, they're separate until I update them.

## Using Photometry.py

Same start point as _master.py_, make sure you've cd'd into _/photomitrus_.  Running the help command: 

    python ./photometry/photometry.py -h

will present helpful text on what each argument is, though the amount of them can be daunting.  The formatting is as follows:

> photometry.py [-h] [-exp_query] [-plots] [-grb] [-dir DIR] [-name NAME] [-filter FILTER] [-survey SURVEY] [-crop CROP] [-RA RA] [-DEC DEC] [-thresh THRESH]

Let us again start with the required positional arguments:

### Required Arguments

- _-dir_ is a string that denotes the parent directory where the co-added final image is stored.

- _-name_ is a string that denotes the filename of the co-added final image.  Note, this is not the full path, but just the filename.  If you utilized _stack.py_ to create the image, the format of the name should be (see _stack.py_'s documentation for more info):

  > coadd.Open-#.12345678-12345678.C#.fits

- _-filter_ is a single character string which denotes the filter of the co-added final image.  This should be clear from the filename of the image.  For instance, 'Open-J' is J band, so we would put 'J' for this argument.

- _-survey_ is a string which denotes which survey _photometry.py_ will be using for calibration.  Currently, there are 2 surveys supported: 2MASS and the VISTA Hemisphere Survey (VHS).  VHS is a deeper survey and generally a better pick for photometry, but it is limited in area to targets in the southern hemisphere, so keep this in mind or you may get an error.  2MASS is all-sky, so it should be fine anywhere.  To utilize 2MASS, put '2MASS', to utilize VISTA, put 'VHS'.

- _-crop_ in an integer which denotes how far in pixels from the edge of the stacked image the script should crop out sources.  This is used to combat the empty edges present in stacked images, along with bad edge artifacts which would result in many false positive sources.  The default is 300 pixels.

### Optional Arguments

The previous were all required arguments, now we will discuss the optional arguments.

- _-exp_query_ is an optional flag that will make _photometry.py_ export an .ecsv file of the survey query results.  This file will contain that survey's information about any sources that were found, excluding the cropped edges.  This can be useful when comparing PRIME's calculated magnitudes to 2MASS or VHS, beyond the scope of what the script already does.

- _-plots_ is an optional flag that will currently generate 2 plots.  The first is a magnitude comparison plot.  It plots the calculated PRIME mags and survey mags of crossmatched sources against eachother.  This can be a good, quick, visual representation of how well the photometry is doing.  For a more in depth look, the 2nd plot is a residuals plot.  A linear fit is run on the crossmatched mag data and then residuals are calculated and plotted.  The residual plot contains useful information, such as lines representing &#963; values, fitted slope & intercept values and error, along with R<sup>2</sup> and RSS values.

- _-grb_ is an optional flag to be used when you're trying to find a specific source.  With this flag, you input the RA and Dec of the GRB you wish to find, along with a search threshold diameter.  _photometry.py_ will search for a sextracted source within that area.  If it finds one, it will print information about the source to the command line, along with generating an .ecsv file with the same information.  This info includes the source's RA & Dec, calculated magnitude & error, 50% flux radius, and SNR.  Below are the optional arguments that must be inputted when using this flag.
    - _-RA_: GRB RA
    - _-DEC_: GRB Dec
    - _-thresh_: Search diameter in arcsec, default = 4"

Now that all of that is out of the way, here is a sample command using all the arguments:

    python ./photometry/photometry.py -plots -grb -dir /../../J_Band/stack/ -name coadd.Open-J.00654321-00123456.C1.fits -filter J -survey 2MASS -crop 400 -RA 221.248375 -DEC -60.697972 -thresh 2.0

## How Photometry.py Works

- To begin, _photometry.py_ reads the WCS of the co-added final image, in order to find the RA and Dec at the center of the image, along with the width and height of the image in arcmin.  (Though, due to the empty edges of the stacked images having no direct coordinate data, these width and height values can be inaccurate and rough, so currently we rely on the edge cropping of sources later on).  
- After getting the size and center of the image in sky coordinates, it then uses astroquery to remotely run a box query the chosen catalog for crossmatching and calibration.
- Sextractor is then run on the co-added image in order to initially extract sources.  Once this catalog is generated, PSFex is run on the catalog to create a PSF model for sources in the catalog.  Then, this PSF model is run back through sextractor, performing model-fitting on sources and calculating PSF-fit magnitudes.
- Once this final sextractor catalog is generated, both it and the previous survey query catalogue have sources on the edges cropped.  The two catalogues are then crossmatched, with a crossmatch being defined as sources in both catalogues that are within 1" of eachother.  _photometry.py_ prints the amount of crossmatches to the command line.
  - If you included _-exp_query_, it is at this point that the cropped survey query catalogue is saved.
- Each crossmatched source then has their sextractor calculated mags subtracted from the catalog source mags.  These offsets are used to calculate the zero point.  This is done by calculating the median and standard deviation of the offset data, excluding sources > 3&#963; from the median.
- For each source, the zero point is then added to the sextractor mags to generate the final filter mag.  The error in this mag is calculated by combining the sextractor mag error and offset data stdev in quadrature.  These mags are then written to columns and added as an addition to the sextractor columns, before being written to a new .ecsv file.  (The _VIGNET_ column is removed from the sextractor catalogue before this writing)
    - After this point, your photometry is done! If you've included _-grb_, then after this file is written, it is read back in and an attempt at a crossmatch is made using your inputted RA, Dec, and thresh.  If it finds a source within these parameters, it will give source information and write it to an .ecsv file.  If not, then it'll let you know.
    - If you included _-plots_, then it will also read in the final sextractor .ecsv and the cropped survey query in order to generate the mag comparison and residual plots.
 
## Data Products

After running this script, the main product that will be outputted will be the .ecsv file with the corrected mags and sextractor columns.  If we call the filename of the co-added image: **COADD** (just for simplicity), then the naming format of the .ecsv would be:

> COADD.fits._survey_.ecsv

Here, _survey_ is the chosen survey name.
There are several other catalogues that are outputted, including the initial sextractor catalogue, the psf model file, and the final model-fitted sextractor catalogue (look at 'How Photometry.py Works' section for more info).  The format for these files, respectively, is:

> COADD.fits.cat,
> COADD.fits.psf,
> COADD.psf.cat

As we don't really utilize these catalogues once the photometry is done, so I could have an option to delete these after the fact.  But for now, they stay.

In addition, the optional args can output more products.  In the case for _-exp_query_, the outputted query catalogue will be (where C# is the chip number):

> _survey_-C#-query.ecsv

For _-plots_, the outputted mag comparison and residual plots will be PNGS with the respective names:

> _survey_-mag_comp_plot.png,
> _survey_-residual_plot.png

Finally, for _-grb_, the .ecsv with the found source's information will be named:

> GRB_Data.ecsv
