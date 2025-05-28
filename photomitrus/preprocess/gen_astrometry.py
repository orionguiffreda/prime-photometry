"""
Runs astrom.net on specified file
"""
# %%
import subprocess
import argparse
import sys
from astropy.io import fits
# sys.path.insert(0,'C:\PycharmProjects\prime-photometry\photomitrus')
from photomitrus.settings import flist
import os


# TODO test this script on PC01 before full implementation (change stuff w/ cfg path)
# %%


def astrom(outpath, inlist, rad, ds):
    for imgpath in inlist:
        img = fits.open(imgpath)
        hdr = img[0].header
        ra = hdr['RA-D']
        dec = hdr['DEC-D']
        try:
            command = (('solve-field '
                        '--backend-config /home/alex/miniconda3/pkgs/astrometry-0.97-py313h139ab80_2/share/astrometry/astrometry.cfg '
                        '--scale-units arcsecperpix --scale-low 0.45 --scale-high 0.55 --ra %s --dec %s --radius %s '
                        '--cpulimit 60 -U none --axy list.axy -S none -M none -R none -B none -O -p -t 4 -z %s -D %s %s') % (
                           ra, dec, rad, ds, outpath, imgpath))
            print('Executing command: %s' % command)
            subprocess.run(command.split(), check=True)
        except subprocess.CalledProcessError as err:
            print('Could not run with exit error %s' % err)
    print('fields solved, all done!')


def astromdir(outpath, directory):
    inpath = sorted(os.listdir(directory))
    for f in inpath:
        try:
            command = ('solve-field '
                       '--backend-config /home/alex/miniconda3/pkgs/astrometry-0.97-py313h139ab80_2/share/astrometry/astrometry.cfg '
                       '--scale-units arcsecperpix --scale-low 0.45 --scale-high 0.55 --no-verify -U none --axy none '
                       '-S none -M none -R none -B none -O -p -t 4 -z 4 -D %s %s') % (
                          outpath, directory + f)
            print('Executing command: %s' % command)
            subprocess.run(command.split(), check=True)
        except subprocess.CalledProcessError as err:
            print('Could not run with exit error %s' % err)
    print('fields solved, all done!')


# astromdir for hard to solve fields, change ra and dec currently


def astromdirhard(outpath, directory, rad, ds):
    inpath = sorted(os.listdir(directory))
    for f in inpath:
        if f.endswith('.fits'):
            img = fits.open(os.path.join(directory, f))
            hdr = img[0].header
            ra = hdr['RA-D']
            dec = hdr['DEC-D']
            try:
                command = (('solve-field '
                            '--backend-config /home/alex/miniconda3/pkgs/astrometry-0.97-py313h139ab80_2/share/astrometry/astrometry.cfg '
                            '--scale-units arcsecperpix --scale-low 0.45 --scale-high 0.55 --ra %s --dec %s --radius %s '
                            '--cpulimit 60 -U none --axy list.axy -S none -M none -R none -B none -O -p -t 4 -z %s -D %s %s') % (
                           ra, dec, rad, ds, outpath, directory+f))
                print('Executing command: %s' % command)
                subprocess.run(command.split(), check=True)
            except subprocess.CalledProcessError as err:
                print('Could not run with exit error %s' % err)
    print('fields solved, all done!')


def gen_astrom(output, input, rad=1, ds=4, soft=False):
    if type(input) is list:
        astrom(output, input, rad, ds)
    elif soft:
        astromdir(output, input)
    else:
        astromdirhard(output, input, rad, ds)


# %%
def main():
    parser = argparse.ArgumentParser(description='Runs astrom.net on specified file')
    parser.add_argument('-soft', action='store_true',
                        help='run on directory w/ looser settings (default has many tightened settings, such as '
                             'center RA and DEC)')
    parser.add_argument('-output', type=str, help='[str] output path for ramps w/ astrom.')
    parser.add_argument('-input', type=str, help='[str] input path for ramps (or list if using that)')
    parser.add_argument('-rad', type=str,
                        help='radius in deg for astrometry searching, default = 1', default=1)
    parser.add_argument('-ds', type=str, help='amount of downsample, default = 4', default=4)
    args, unknown = parser.parse_known_args()

    gen_astrom(args.output, args.input, args.rad, args.ds, args.soft)


if __name__ == "__main__":
    main()
