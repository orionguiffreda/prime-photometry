conda create -n prime-photometrus numpy scikit-image scipy pandas astropy statsmodels ipython matplotlib astroquery astromatic-psfex astromatic-scamp astromatic-source-extractor astromatic-swarp astrometry pysftp json5 conda-build -c conda-forge
python build.py
conda develop .
ln -s bin/photometrus $CONDA_PREFIX/bin/photometrus
