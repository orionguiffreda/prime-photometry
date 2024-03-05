# PRIME Image Processing & Photometry Pipeline

This is the photometry pipeline for PRIME telescope.  

## Installation instructions for PRIME Pipeline

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

