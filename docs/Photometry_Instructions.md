# PRIME Photometry Pipeline - Photometry Instruction and Basic Usage

## Before using Photometry.py

Before using this script, make sure you already possess a co-added final image, likely created through usage of _master.py_.  Currently, _master.py_ and _photometry.py_ are separate, requiring 1 command to process into a final image and another command to run photometry.  I could (and plan to) combine them into 1 big script that does everything in 1 go.  However, while combining the basic functionality of _photometry.py_ with _master.py_ would be relatively easy, due to the many current (and future) optional arguments for _photometry.py_ I don't want to have a very cluttered single command that would be prone to misinputting, leading to many minutes of lost time.  So, for now, they're separate until I update them.

## Using Photometry.py

Same start point as _master.py_, make sure you've cd'd into _/photomitrus_.  Running the help command: 

    python ./photometry/photometry.py -h

will present helpful text on what each argument is, though the amount of them can be daunting.  The formatting is as follows:

> photometry.py [-h] [-exp_query] [-plots] [-grb] [-dir DIR] [-name NAME] [-filter FILTER] [-survey SURVEY] [-crop CROP] [-RA RA] [-DEC DEC] [-thresh THRESH]

Let us again start with the required positional arguments:

- _-dir_ is a string that denotes the parent directory where the co-added final image is stored.

- _-name_ is a string that denotes the filename of the co-added final image.  Note, this is not the full path, but just the filename.  If you utilized _stack.py_ to create the image, the format of the name should be (see _stack.py_'s documentation for more info):

  > coadd.Open-*.12345678-12345678.C*.fits

- _-filter_ is a single character string which denotes the filter of the co-added final image.  This should be clear from the filename of the image.  For instance, 'Open-J' is J band, so we would put 'J' for this argument.

- _-survey_ is a string which denotes which survey _photometry.py_ will be using for calibration.  Currently, there are 2 surveys supported: 2MASS and the VISTA Hemisphere Survey (VHS).  VHS is a deeper survey and generally a better pick for photometry, but it is limited in area to targets in the southern hemisphere, so keep this in mind or you may get an error.  2MASS is all-sky, so it should be fine anywhere.  To utilize 
