#!/bin/sh
import sys
from photomitrus.photometry import photometry

if sys.argv[0] == 'photometry':
    photometry.main()

# if sys.argv[1] == 'stack':

