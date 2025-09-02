conda create -n prime-photometrus numpy scikit-image scipy pandas astropy statsmodels ipython matplotlib astroquery astromatic-psfex astromatic-scamp astromatic-source-extractor astromatic-swarp astrometry regions pysftp json5 conda-build -c conda-forge
conda activate prime-photometrus
python build.py
conda develop .
CURRENT_WORKING_DIR=$(pwd)
ln -s $CURRENT_WORKING_DIR/bin/photometrus $CONDA_PREFIX/bin/photometrus
