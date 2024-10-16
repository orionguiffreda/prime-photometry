# PRIME Photometry Pipeline - Photometry Instruction and Basic Usage

## Before using Photometry (Photometry.py)

Before using this script, make sure you already possess a co-added final image, likely created through usage of _photometrus pipeline_.  Currently, _photometrus pipeline_ (_multi_master.py_) and _photometrus photometry_ (_photometry.py_) are separate, requiring 1 command to process into a final image and another command to run photometry.  I could (and plan to) combine them into 1 big script that does everything in 1 go.  However, while combining the basic functionality of both would be relatively easy, due to the many current (and future) optional arguments for _photometrus photometry_ I don't want to have a very cluttered single command that would be prone to misinputting, leading to many minutes of lost time.  So, for now, they're separate until I update them.

## Using Photometry

Running the help command: 

    photometrus photometry -h

will present helpful text on what each argument is, though the amount of them can be daunting.  The formatting is as follows:

> photometry.py [-h] [-filepath FILEPATH] [-band BAND] [-survey SURVEY] [-crop CROP] [-grb] [-grb_only] [-grb_ra GRB_RA] [-grb_dec GRB_DEC] [-grb_radius GRB_RADIUS] [-grb_coordlist ...] [-exp_query] [-exp_query_only] [-no_plots] [-plots_only] [-keep] 

Let us again start with the required positional arguments:

### Required Arguments

- _-filepath_ is a string that denotes the full filepath to the stacked image.  One can also simply put the filename if they are cd'd into the directory with that file.  An example is below:

  > /obs/J/stack/coadd.Open-#.12345678-12345678.C#.fits

- _-band_ is a single character string which denotes the filter of the co-added final image.  This should be clear from the filename of the image.  For instance, 'Open-J' is J band, so we would put 'J' for this argument.

- _-survey_ is a string which denotes which survey _photometry.py_ will be using for calibration.  Currently, there are 6 surveys supported: 2MASS, VISTA Hemisphere Survey (VHS), VIKING, Skymapper, SDSS, and UKIDSS.  VHS is a deeper survey and generally a better pick for photometry, but it is limited in area to targets in the southern hemisphere, so keep this in mind or you may get an error.  In addition, while VHS has nearly all S. hemi. coverage in J (and Ks, but PRIME doesn't have that filter), its coverage in H, Y, and Z is far less.  2MASS is all-sky, so it should be fine anywhere in J and H (though it is shallow).

    To supplement the gaps in both coverage and filters from VHS, we have the other surveys.  VIKING covers some gaps in VHS in both J and H.  Skymapper has Z band coverage in nearly all of the Southern Hemisphere, though is shallower than the equivalent VHS limiting mag in J.  SDSS is much deeper in Z, but covers less area.  Finally, UKIDSS provides Y band coverage, but only in certain areas.  

    Currently, you have to manually specify one of these surveys.  In the future, I plan to create a more automatic query system.  The links below show coverage maps for VISTA and other surveys in each filter.  To utilize 2MASS, put '2MASS', to utilize VISTA, put 'VHS', for the others, put the exact names I mentioned earlier.

    - http://casu.ast.cam.ac.uk/vistasp/overview/
    - https://skymapper.anu.edu.au/data-release/
    - https://www.sdss4.org/science/data-release-publications
    - http://www.ukidss.org/surveys/surveys.html

- _-crop_ in an integer which denotes how far in pixels from the edge of the stacked image the script should crop out sources.  This is used to combat the empty edges present in stacked images, along with bad edge artifacts which would result in many false positive sources.  The default is 300 pixels.

### Optional Arguments

The previous were all required arguments, now we will discuss the optional arguments.

- _-grb_ is an optional flag to be used when you're trying to find a specific source.  With this flag, you input the RA and Dec of the GRB you wish to find, along with a search threshold diameter.  _photometry.py_ will search for a sextracted source within that area.  If it finds one, it will print information about the source to the command line, along with generating an .ecsv file with the same information.  If it finds multiple sources, it will let you know in the command line and log info for all found sources in the .ecsv.  This info includes the source's RA & Dec, calculated magnitude & error, 50% flux radius, SNR, and distance from inputted coords.  Below are the optional arguments that must be inputted when using this flag.
    - _-grb_ra_: GRB RA
    - _-grb_dec_: GRB Dec
    - _-grb_radius_: Search diameter in arcsec, default = 4"
 
- _-grb_coordlist_ is an optional flag to be used with _-grb_.  If a user has many areas they want to search for sources in, instead of running photometry several times with the _-grb_ flag, they can use this field.  Input RA and DEC of each area with a comma between them, and put spaces between each pair of coordinates.  If sources are found in any of the search areas, they will be output to seperate .ecsvs.  For instance the usage for 2 search areas is below.  
    *NOTE* When using this flag, you must still specify _-grb_radius_.  Currently, only 1 thresh is supported for all search areas.   

    -coordlist 123,-45 321,-54

- _-grb_only_ is an optional flag to be used if you want to run the grb functionality again.

- _-exp_query_ is an optional flag that will make _photometry.py_ export an .ecsv file of the survey query results.  This file will contain that survey's information about any sources that were found, excluding the cropped edges.  This can be useful when comparing PRIME's calculated magnitudes to 2MASS or VHS, beyond the scope of what the script already does.

- _-exp_query_only_ is an optional flag which will only generate the survey query results, and not complete photometry.

- _-no_plots_ is an optional flag that will **stop the generation** of several plots.  The _default_ behavior is as follows.  The first plot is a magnitude comparison plot.  It plots the calculated PRIME mags and survey mags of crossmatched sources against eachother.  This can be a good, quick, visual representation of how well the photometry is doing.  For a more in depth look, the 2nd plot is a residuals plot.  A linear fit is run on the crossmatched mag data (y int is forced to zero) and then residuals are calculated and plotted.  The residual plot contains useful information, such as lines representing &#963; values, fitted slope & intercept values and error, along with R<sup>2</sup> and RSS values.  Another residual plot is generated, taking into account the y int.  Finally, a limiting magnitude plot is generated.  This bins all sextracted sources by 0.1 mag, from 12 to 25 mag, plotting number of sources per mag against mag.  The half-max point of this plot (the limiting mag) is printed to command line and presented in the plot.  **Remember, using this flag will STOP this plot generation**

- _-plots_only_ is an optional flag that you can use if photometry has already been generated once.  This will quickly generate the plots only.

- _-keep_ is an optional flag that will **keep** intermediate catalogs that have no use after photometry.  The default behavior is to remove these, as these are only used to generate the photometry and never again.

Now that all of that is out of the way, here is a sample command using all the arguments:

    photometrus photometry -grb -filepath /../../J_Band/stack/coadd.Open-#.12345678-12345678.C#.fits -band J -survey 2MASS -crop 400 -grb_ra 221.248375 -grb_dec -60.697972 -grb_radius 2.0

## How Photometry Works

- To begin, _photometry.py_ reads the WCS of the co-added final image, in order to find the RA and Dec at the center of the image, along with the width and height of the image in arcmin.  (Though, due to the empty edges of the stacked images having no direct coordinate data, these width and height values can be inaccurate and rough, so currently we rely on the edge cropping of sources later on).  
- After getting the size and center of the image in sky coordinates, it then uses astroquery to remotely run a box query the chosen catalog for crossmatching and calibration.
- Sextractor is then run on the co-added image in order to initially extract sources.  Once this catalog is generated, PSFex is run on the catalog to create a PSF model for sources in the catalog.  Then, this PSF model is run back through sextractor, performing model-fitting on sources and calculating PSF-fit magnitudes.
- Once this final sextractor catalog is generated, both it and the previous survey query catalogue have sources on the edges cropped.  The two catalogues are then crossmatched, with a crossmatch being defined as sources in both catalogues that are within 1" of eachother.  _photometry.py_ prints the amount of crossmatches to the command line.
  - If you included _-exp_query_, it is at this point that the cropped survey query catalogue is saved.
- Each crossmatched source then has their sextractor calculated mags subtracted from the catalog source mags.  These offsets are used to calculate the zero point.  This is done by calculating the median and standard deviation of the offset data, excluding sources > 3&#963; from the median.
- For each source, the zero point is then added to the sextractor mags to generate the final filter mag.  The error in this mag is calculated by combining the sextractor mag error and offset data stdev in quadrature.  These mags are then written to columns and added as an addition to the sextractor columns, before being written to a new .ecsv file.  (The _VIGNET_ column is removed from the sextractor catalogue before this writing)
    - After this point, your photometry is done! If you've included _-grb_, then after this file is written, it is read back in and an attempt at a crossmatch is made using your inputted RA, Dec, and thresh.  If it finds a source within these parameters, it will give source information and write it to an .ecsv file.  If not, then it'll let you know.
- Unless you included _-no_plots_, then it will also read in the final sextractor .ecsv and the cropped survey query in order to generate the mag comparison, residual plots, and lim mag plots.
 
## Data Products

After running this script, the main product that will be outputted will be the .ecsv file with the corrected mags and sextractor columns.  If we call the filename of the co-added image: **COADD** (just for simplicity), then the naming format of the .ecsv would be:

> COADD.fits._survey_.ecsv

Here, _survey_ is the chosen survey name.
There are several other catalogues that are outputted, including the initial sextractor catalogue, the psf model file, and the final model-fitted sextractor catalogue (look at 'How Photometry.py Works' section for more info).  The format for these files, respectively, is:

> COADD.fits.cat,
> COADD.fits.psf,
> COADD.psf.cat

We don't really utilize these catalogues once the photometry is done, the _-keep_ flag will keep these if you like.

In addition, the optional args can output more products.  In the case for _-exp_query_, the outputted query catalogue will be (where C# is the chip number):

> _survey_-C#-query.ecsv

For _-plots_, the outputted mag comparison and residual plots will be PNGS with the respective names:

> _survey_-mag_comp_plot.png,
> _survey_-residual_plot.png
> _survey_-residual_plot_int.png
> _survey_chip_-lim_mag_plot.png

Finally, for _-grb_, the .ecsv with the found source's (or sources') information will be named:

> GRB-_filter_-Data-_survey_.ecsv or GRB_Multisource-_filter_-Data-_survey_.ecsv

If you use _-coordlist_ for multiple search areas, the .ecsvs will be output for each area with a found source.  

> GRB-_filter_-Data-_survey_-loc_#.ecsv or GRB_Multisource-_filter_-Data-_survey_-loc_#.ecsv

Where the loc_# is the number corresponding to the pair of coordinates you inputted.  For instance, loc_2 would correspond to the second pair you put in the command.
