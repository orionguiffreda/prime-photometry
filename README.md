![Photometrus logo image](https://github.com/Oriohno/prime-photometry/blob/refactoring_cmd/docs/images/PRIME_logo_v2.png)
# Photometrus PRIME

An image processing and photometry pipeline for PRIME telescope by [Orion Guiffreda](https://github.com/Oriohno) and [Joe Durbak](https://github.com/joedurbak)

# Welcome!

Here is the wiki of the PRIME telescope photometry pipeline!  Here will be pages of information on how to use the pipeline's scripts, how these scripts function, and other various things, such as recommended index and config files.  I'll try my best to explain everything about the pipeline, and will be updating new and current pages over time as functionality improves! 

If you have access to the PRIME computer and want to jump right into how to utilize the pipeline, look at the Quick Start Guide below!  Then, for further information, examine the _photometrus stack_ (_master.py_) and _photometrus pipeline_ (_multi_master.py_) scripts section (as explained below!)

**NOTE!** *A key script for downloading data is not included in the repository, due to security.  As a result (and because it is still in development) the pipeline's greater functions will not operate.*  
If you want to download and set up the pipeline (once it's ready), look at the Installation Guide below!

## Quick Start Guide

### Setting up Access

To begin, if you're using the main pipeline on the PRIME computer, you'll need to access the computer via SSH or TigerVNC.  TigerVNC allows remote desktop access but as we're currently working through some issues with it, it's capabilities are limited.  SSH works as normal.  There is a Google Doc that details how to set up the VPN and access goddardpc01 (I won't link it here as it contains sensitive logins, but if you're in the PRIME ToO slack, it should be a bookmark under 'PRIME data retrieval').  Remember, the pipeline is stored on goddardpc01 specifically, so access that computer only.

### Confirming an Observation

Confirm that you can access goddardpc01.  Before we move to utilizing the pipeline itself, we should first confirm the observation we want to process.  It's good practice to verify that your observation target has been observed in the date and filter you think it is.  You can verify by checking either the PRIME ramp log or the new online log at the links below:

http://www-ir.ess.sci.osaka-u.ac.jp/prime_staff/LOG/Ramp_LOG/

http://www-ir.ess.sci.osaka-u.ac.jp/prime_staff/Online_Log/

Simply navigate to your date and search for your field and filter.

For the purposes of this guide, let us assume a transient has been observed.  We'll use the recent observation of GRB250309B as an example.  This was observed by PRIME on March 10, 2025.  We have the date, now we need the field number (OBJNAME in the ramp log).  This corresponds to where on the PRIME observing grid the target was observed.  When you learn to create observation CSVs to schedule observation, this corresponds to the *ObjectName* column.  Currently, observation CSVs are submitted via a Google submission portal, so they are not easily accessible by people other than the submitter.  A good piece of advice would be to check the ramp log for the field (or fields) on your certain night which were observed by 'NASA' in the *Observer* column.  In the case of this transient, the field number is *no_grid*, as the transient's position caused us to shift off the grid to cover the localization radius.  Most fields will have the format: *field12345*.  

Let's just try and process J band for the filter.  Most PRIME observations are taken in J and or H, sometimes with Y and or Z.  So now we have the necessary date, field, and filter.  Look at the corresponding log to determine if the information is accurate (the observation should be there!).  Once you've confirmed the observation was taken, we can move on.

### Utilizing the Pipeline (Processing & Photometry)

Now we can finally start using the pipeline!  SSH into goddardpc01.  

Currently, to get the pipeline ready for use, let's begin by activating the correct conda environment:

    conda activate prime-photometrus

To utilize all the scripts in this pipeline, you'll call the main command _photometrus_.  There are many, many scripts associated with this pipeline, and _photometrus_ allows access to nearly all of them.  Specifically, we'll begin by using the script most often utilized in the pipeline: _multi_combo.py_.  We can call this script through calling the command _photometrus full_.  This allows the processing, stacking, and photometric analysis of a target field for 1 or more detectors (chips) all from a single command.  This script takes several fields (such as date, band, etc.) as input, if you want detailed explanations of every argument, go to the Scripts section.

We first need to determine what parent directory all the processing will take place in.  Preferably, it should be a new directory in the _/mnt/photometry/TransientEvents/_ path (this is where most GRB observations are stored).  For the sake of this guide, let's make the directory:
_/mnt/photometry/TransientEvents/pipeline_demo/J/_.

This script has many optional args.  For the sake of this quickstart guide, I won't go over them, but explore them through the documentation or running: 

    photometrus full -h

Let's run the pipeline only on chip 1 in J band, as this is where the target landed.  To run the pipeline on this observation, we should utilize the command:

    photometrus full -parent /mnt/photometry/TransientEvents/pipeline_demo/J/ -date 20250310 -target no_grid -band J -chip 1 -grb_ra 210.80129 -grb_dec -8.50302 -grb_radius 4.0

This is a pretty long command, so let's go over a few things briefly.  Some args are self-explanatory (-parent is the parent directory, -date is the date, -band is the filter, etc.).  
Remember, this is a transient observation, so not only are we generating photometry for the whole image, we're also looking for a specific source.  The RA, Dec, and error radius of this transient is inputted at -grb_ra, -grb_dec, and -grb_radius.  We got the coordinates and error radius from the GCN network, specifically from Swift through the GCN below (I increased the error radius slightly):

https://gcn.nasa.gov/circulars/39649

Once you run this command, you'll notice the script is quite verbose.  This is good for monitoring progress, as it details exactly what is occuring on every step.  

_Ideally_, this command should run without issue, producing many subdirectories (the pipeline currently keeps all intermediate data products by default, useful for troubleshooting errors).  For a detailed overview of each of these subdirectories and data products, examine the _photometrus stack_ (_master.py_) documentation.  But quickly, below should be the format:

    ├── J (Parent)
    │   ├── C1_astrom  -  Subdirectory for storage of ramps w/ basic astrometry
    │   ├── C1_sub  -  Subdirectory for storage of processed images w/ improved astrometry
    │   ├── sky  -  Subdirectory for storage of sky image
    │   ├── stack  -  Subdirectory for storage of final stacked image
    │   ├── ramp_fit_log_****-**-**.dat  -  PRIME observation logs for the night
    │   ├── ramp_fit_log_****-**-**.clean.dat

The final stacked image and it's photometric information is the one we're interested in.  If you're using TigerVNC you can open Dolphin, the file viewer, then open the image in DS9 and examine it.  Here, the corrected mounted drive should be titled _wsldata_.  If you're just SSH'd in, you'll have to scp it.  Hopefully it looks acceptable! (no star streaking or blurriness).  Within the _stack_ subdirectory, there should be many files corresponding to photometry.

    ├── stack
    │   ├── coadd.Open-J.02182248-02182886.C1.fits - final stacked image
    │   ├── weight.Open-J.02182248-02182886.C1.fits - weightmap for stacked image
    │   ├── coadd.Open-J.02182248-02182886.C1.fits.VHS.ecsv - full photometric catalog
    │   ├── coadd.Open-J.02182248-02182886.C1.wcs - external WCS header generated by astrometry.net
    │   ├── PSF.Open-J.02182248-02182886.C1.fits - PSF model fits file
    │   ├── Resid_3-sig_Data_J_C1_VHS.ecsv - residual statistics file
    │   ├── GRB_J_Data_VHS.ecsv - transient target information file
    │   ├── VHS_C1_*.png... - Many check-plots for photometry

We won't get into each data product (check the _photometrus photometry_ (_photometry.py_) documentation for more info), but the one we're interested in is the file titled _GRB_J_Data_VHS.ecsv_ is the one we're interested in.  If you open it, you'll see it contains various info on a source at the inputted coordinates and threshold, such as magnitude, radius, SNR, etc.  Congrats! This is our grb! 

## Sections

### Scripts: 
Most of the wiki will be about the scripts, how they function, what they produce, and how to use them.  *NOTE*: Not every utilized script currently has complete documentation, but it will be updated over time!  

* Photometrus Stack & Photometrus Pipeline: The pipeline itself will be run from either _master.py_ or _multi_master.py_ (called by _stack_ and _pipeline_ respectively).  For most use case scenarios and normal observations _multi_master.py_ is sufficient, though specific scenarios may require the greater number of knobs to turn given by _master.py_.  In any case, the sections on these & the photometry scripts are the most important, though other script documentation is useful for knowing how the master scripts run.  To begin to understand how to run things, I recommend first reading through _multi-master.py_ (utilized the most), then _master.py_, and finally the other scripts.

* Photometrus Photometry & Photometrus Single_Photometry: The photometry will be run either from _photometry.py_ or _multi_photom.py_ (called by _single_photometry_ and _photometry_ respectively).  The former runs photometry on all or specific chips in a given observations, while the latter runs on a single chip, but can offer more knobs to turn.

* Photometrus Astrom: This section details each of the 4 current astrometry scripts that can be used in the pipeline.

* Photometrus Process: This section will include details on processing scripts responsible for specific steps, such as stacking.

Please refer to the Wiki tab for the above pages!


## Installation Guide for PRIME Pipeline

### Clone the repo

1) Navigate to the main page of this repo: https://github.com/Oriohno/prime-photometry.  Then click the Code button and copy the HTML or SSH link.
2) Open command line and change your current working directory (cd) to the place you wish to clone the repo.
3) Finally, use the command 'git clone' as below to clone locally.

        git clone https://github.com/Oriohno/prime-photometry.git

### Conda environment setup

1) Begin by navigating to the .yml file included in this repo.  cd to to _/prime-photometry/_ for ease of using the command below.
2) To create the conda environment, run the command below.

        conda env create -f prime-photometry.yml
   
4) You can then activate the environment by using:

        conda activate prime-photometry

## Index file installation

### About index files

To utilize a part of the pipeline, astrometry.net, you will need index files from which the package can pull from.  Astrometry.net is frustratingly inconsistent, and the pipeline normally only uses it for a single final (technically optional) astrometry step.  Thus, it is possible to utilize the pipeline without it.  However, if you want to use the pipeline, you'll need to download the index files.  See these links for key details about and where to download index files:

https://astrometry.net/doc/readme.html

http://data.astrometry.net/

For PRIME I found the 4200 series is useful, as I first downloaded those and they've worked (mostly) just fine, being built off of 2MASS.  Astrometry.net recommends to download index files with sky-marks 10% to 100% the size of your images.  For PRIME, individual detector images are in the 30-40' size range, thus I should download series _08_ to _01_.  I've done this for some areas, but have yet to complete coverage.  

If you want to save time (and delay the inevitable), you can download appropriate index files only for the area where your images are located.  You can use these maps to determine what numbers correspond to what index files: 

https://github.com/dstndstn/astrometry.net/blob/master/util/hp.png

https://github.com/dstndstn/astrometry.net/blob/master/util/hp2.png

It is also recommended to download the 4100 and some of the 5200 series, though I haven't yet tested if these are better for PRIME than the 4200.

### Installation placement

We must place the index files in its own directory, named whatever you like.  Find the 'astrometry.cfg' file, likely located where the astrometry python package is installed.  It should look like this below:

https://github.com/dstndstn/astrometry.net/blob/main/etc/astrometry.cfg

Under the comment "# In which directories should we search for indices?", add your index file directory after "add_path", like the example below:

        # In which directories should we search for indices?
        add_path ../../../indexes

Now, astrometry.net should be able to use your index files.  
