# Welcome!

Here is the wiki of the PRIME telescope photometry pipeline!  Here will be pages of information on how to use the pipeline's scripts, how these scripts function, and other various things, such as recommended index and config files.  I'll try my best to explain everything about the pipeline, and will be updating new and current pages over time as functionality improves! 

To jump right into how to utilize the pipeline, look at the Quick Start Guide below!  Then, for further information, examine the _photometrus stack_ (_master.py_) and _photometrus pipeline_ (_multi_master.py_) scripts section (as explained below!)

## Quick Start Guide

### Setting up Access

To begin, if you're using the main pipeline on the PRIME computer, you'll need to access the computer via SSH or Microsoft Remote Desktop.  I personally use Remote Desktop, as it makes it easy for me to examine scripts, files, and results visually.  There is a Google Doc that details how to set up the VPN and access goddardpc01 (I won't link it here as it contains sensitive logins, but if you're in the PRIME ToO slack, it should be a bookmark under 'PRIME data retrieval').  Remember, the pipeline is stored on goddardpc01 specifically, so access that computer only.

### Confirming an Observation

Confirm that you can access goddardpc01.  Before we move to utilizing the pipeline itself, we should first confirm the observation we want to process.  It's good practice to verify that your observation target has been observed in the date and filter you think it is.  You can verify by checking the PRIME ramp log at the link below:

http://www-ir.ess.sci.osaka-u.ac.jp/prime_staff/LOG/Ramp_LOG/

Simply navigate to your date and search for your field and filter.

For the purposes of this guide, let us assume a transient has been observed.  We'll use GRB240825a as an example.  This was observed by PRIME on 8/25/2024.  We have the date, now we need the field number (OBJNAME in the log).  If you're in the PRIME discord, you should be able to navigate to the _obs_request_ channel.  This is where observation request csvs are sent to the observers.  You can scroll until you find the csv titled 'GRB240825a.csv'.  Within this csv, you can find the *ObjectName* column.  This corresponds to the field on the observation grid that should be observed, and is the field number we need.  If you aren't currently in the PRIME discord, then the field number for this GRB is _field16359_.  

Let's just try and process J band for the filter.  So now we have the necessary date, field, and filter.  Look at the corresponding log to determine if the information is accurate (the observation should be there!).  Once you've confirmed the observation was taken, we can move on.

### Utilizing the Pipeline

Now we can finally start using the pipeline!  Remote into goddardpc01 and open up WSL (it should be on the taskbar).  

Currently, to get the pipeline ready for use, let's begin by activating the correct conda environment:

    conda activate prime-photometrus

To utilize all the scripts in this pipeline, you'll call the main command _photometrus_.  Specifically, we'll begin by using the script most often utilized in the pipeline: _multi_master.py_.  We can call this script through calling the command _photometrus pipeline_.  This allows the data download, processing, and stacking for 1 or more detectors (chips) all from a single command.  This script takes several fields (such as date, band, etc.) as input, if you want detailed explanations of every argument, go to the Scripts section.

We first need to determine what parent directory all the processing will take place in.  Preferably, it should be a new directory in the _/mnt/d/PRIME_photometry_test_files/_ path (this is where most observations are stored).  For the sake of this guide, let's make the directory:
_/mnt/d/PRIME_photometry_test_files/pipeline_demo/_.

This script has many optional args.  For the sake of this quickstart guide, I won't go over them, but explore them through the documentation or running: 

    photometrus pipeline -h

Let's run the pipeline only on chip 2 in J band, to save time and disk space.  To run the pipeline on this observation, we should utilize the command:

    photometrus pipeline -parent /mnt/d/PRIME_photometry_test_files/pipeline_demo/ -target field16539 -date 20240825 -band J -chip 2

Once you run this command, you'll notice the script is quite verbose in WSL.  This is good for monitoring progress, as it details exactly what is occuring on every step.  

_Ideally_, this command should run without issue, producing many subdirectories (the pipeline currently keeps all intermediate data products, useful for troubleshooting errors).  For a detailed overview of each of these subdirectories and data products, examine the _photometrus stack_ (_master.py_) documentation.  But quickly, below should be the format:

    ├── J_Band (Parent)
    │   ├── C1  -  Subdirectory for storage of input ramp images
    │   ├── C1_astrom  -  Subdirectory for storage of ramps w/ basic astrometry
    │   ├── C1_sub  -  Subdirectory for storage of processed images w/ improved astrometry
    │   ├── sky  -  Subdirectory for storage of sky image
    │   ├── stack  -  Subdirectory for storage of final stacked image
    │   ├── ramp_fit_log_****-**-**.dat  -  PRIME observation logs for the night
    │   ├── ramp_fit_log_****-**-**.clean.dat

The final stacked image is the one we're interested in.  Open the image in DS9 and examine it.  Hopefully it looks acceptable! (no star streaking or blurriness).  Once we've confirmed the image is of good quality, let's move onto photometry.

### Utilizing Photometry

We have our image, now it's time to get photometric information from it.  To utilize _photometry.py_, we will call photometrus through the command _photometrus photometry_ (self explanatory).  For this observation, we'll utilize the command:

    photometrus photometry -filepath /mnt/d/PRIME_photometry_test_files/pipeline_demo/stack/coadd.Open-J.01599131-01599329.C2.fits -band J -survey VHS -grb_ra 344.57200 -grb_dec 1.02675 -grb_radius 5.0

This is a pretty long command, so let's go over a few things briefly.  Some args are self-explanatory (-filepath being the path to the stacked image, -band being the filter).  

-survey is the catalog that is used for source matching.  We match PRIME sources in the image to known catalog sources to generate accurate magnitudes.  

Remember, this is a transient observation, so not only are we generating photometry for the whole image, we're looking for a specific source.  The RA, Dec, and error radius of this transient is inputted at -grb_ra, -grb_dec, and -grb_radius.

There are many optional args with this script.  To explore them and how the script works in general, look into the documentation and or run the command:

    photometrus photometry -h

Anyway, once the long command I gave earlier is run, you'll see some info printed to WSL, and several data products created.  As we're looking for a transient, the file titled _GRB_J_Data_VHS.ecsv_ is the one we're interested in.  If you open it in Notepad, you'll see it contains various info on a source at the inputted coordinates and threshold, such as magnitude, radius, SNR, etc.  Congrats! This is our grb!  

## Sections

### Scripts: 
Most of the wiki will be about the scripts, how they function, what they produce, and how to use them.  *NOTE*: Not every utilized script currently has complete documentation, but it will be updated over time!  

* Photometrus Stack & Photometrus Pipeline: The pipeline itself will be run from either _master.py_ or _multi-master.py_ (called by _stack_ and _pipeline_ respectively).  For most use case scenarios and normal observations _multi-master.py_ is sufficient, though specific scenarios may require the greater number of knobs to turn given by _master.py_.  In any case, the sections on these scripts are the most important, though other script documentation is useful for knowing how the master scripts run.  To begin to understand how to run things, I recommend first reading through _multi-master.py_ (utilized the most), then _master.py_, and finally the other scripts.

* Astrometry: This section details each of the 4 current astrometry scripts that can be used in the pipeline.

* Photometry: This section details how photometry is generated for a chip's data, and what data products are produced.

* Stack: This section details how the final pipeline step, stacking, works.


## Installation instructions for PRIME Pipeline (taken from initial readme)

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

To utilize a key part of the pipeline, astrometry.net, you will need index files from which the package can pull from.  See these links for key details about and where to download index files:

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

