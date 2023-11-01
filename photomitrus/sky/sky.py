"""
Subtracts sky
"""
import os

import numpy as np
from astropy.io import fits
#%%
os.chdir('/Users/orion/Desktop/PRIME/GRB/H_band/H_flat')
img = fits.getdata('/Users/orion/Desktop/PRIME/GRB/H_band/H/00406203C2.cds.fits')
sky = fits.getdata('/Users/orion/Desktop/PRIME/GRB/H_band/H_sky/sky.Open-H.00406203-00406559.C2.fits')
hdl = fits.open('/Users/orion/Desktop/PRIME/GRB/H_band/H_flat/H.flat.C2.fits')
flatdata = hdl[0].data
hdr = hdl[0].header
redflat = flatdata[4:4+4088,4:4+4088]
flat = (redflat-sky)/np.nanmedian(redflat-sky)
#flat = sky/sky.mean()


#%%
reduce = ((img-sky)/flat)

fits.writeto('testimg.fits',reduce,hdr,overwrite=True)


#%%
def subtract_sky_and_normalize(science_data_directory, output_data_dir, sky, flat):
    if not os.path.isdir(output_data_dir):
        os.makedirs(output_data_dir)
    image_fnames = [os.path.join(science_data_directory, f) for f in os.listdir(science_data_directory) if f.endswith('.cds.fits')]
    image_fnames.sort()
    sky = fits.getdata(sky)
    flatdata = fits.getdata(flat)
    flat = flatdata[4:4+4088,4:4+4088]
    for f in image_fnames:
        print(f)
        with fits.open(f) as hdul:
            image = hdul[0].data
            header = hdul[0].header
        reduced_image = (image-sky) / flat
        output_fname = os.path.basename(f)
        output_fname = output_fname.replace('.fits', '.sky.flat.fits')
        output_fname = os.path.join(output_data_dir, output_fname)
        print(output_fname)
        fits.HDUList(fits.PrimaryHDU(header=header, data=reduced_image)).writeto(output_fname, overwrite=True)

#%%

subtract_sky_and_normalize('/Users/orion/Desktop/PRIME/GRB/H_band/H/', '/Users/orion/Desktop/PRIME/GRB/H_band/H_preproctest/',
                           '/Users/orion/Desktop/PRIME/GRB/H_band/sky.3.5.fits', '/Users/orion/Desktop/PRIME/GRB/H_band/H_flat/00384450_00384446C2.fits')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generates sky for given filter and dataset')
    parser.add_argument('paths', nargs='*', type=str, metavar='p', help='Put input path first, output path, sky path,'
                                                                        'then flat path')
    args = parser.parse_args()
    subtract_sky_and_normalize(args.files[0],args.files[1],args.files[2],args.files[3])