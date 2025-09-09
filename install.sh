conda create -n prime-photometrus numpy scikit-image scipy pandas astropy statsmodels ipython matplotlib astroquery astromatic-psfex astromatic-scamp astromatic-source-extractor astromatic-swarp astrometry regions pysftp json5 cython fitsio conda-build -c conda-forge
conda activate prime-photometrus
cd photometrus/ramp || exit
python setup.py build_ext --inplace
cd ../..
python build.py
conda develop .
ln -s "$(pwd)/bin/photometrus" "$CONDA_PREFIX/bin/photometrus"
