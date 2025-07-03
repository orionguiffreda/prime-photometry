"""
Runs improved astrom on a given file using sextractor and scamp
"""

import os
import sys
import argparse
import subprocess
import sys
from astropy.io import fits
import threading
from datetime import datetime as dt

# sys.path.insert(0,'C:\PycharmProjects\prime-photometry\photomitrus')
from photomitrus.settings import gen_config_file_name

#%%


def sextract(imgdir, f, sx, ap):
    pre = os.path.splitext(f)[0]
    ext = os.path.splitext(f)[1]
    com = ["sex ", os.path.join(imgdir, pre + ext), ' -c ' + sx, " -CATALOG_NAME " + pre + '.cat', ' -PARAMETERS_NAME ' + ap]
    s0 = ''
    com = s0.join(com)
    out = subprocess.Popen([com], shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    out.wait()
    # print(pre + ext + ' sextracted!')


def sex(imgdir, bulge=False):
    os.chdir(str(imgdir))
    if bulge:
        sx = gen_config_file_name('bulge_new.config')
        ap = gen_config_file_name('tempsource.param')
    else:
        sx = gen_config_file_name('sex.config')
        ap = gen_config_file_name('astrom.param')
    threads = []
    print('Sextracting all images in %s' % imgdir)
    for f in sorted(os.listdir(imgdir)):
        if f.endswith('flat.fits') or f.endswith('.new'):
            thread = threading.Thread(target=sextract, args=(imgdir, f, sx, ap), name=f)
            thread.start()
            threads.append(thread)
    for thread in threads:
        thread.join()
    print('Sextraction of all images complete!')


def sexback(imgdir):
    os.chdir(str(imgdir))
    sx = gen_config_file_name('sex.config')
    ap = gen_config_file_name('astrom.param')
    for f in sorted(os.listdir(str(imgdir))):
        if f.endswith('ramp.fits'):
            pre = os.path.splitext(f)[0]
            ext = os.path.splitext(f)[1]
            com = ["sex ", imgdir + pre + ext, ' -c '+sx, " -CATALOG_NAME " + pre + '.cat', ' -PARAMETERS_NAME '+ap,
                   ' -CHECKIMAGE_NAME '+pre+'.back.fits']
            s0 = ''
            com = s0.join(com)
            out = subprocess.Popen([com], shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            out.wait()
            print(pre + '.back.fits back subbed!')


def scamp(imgdir, distortdeg=None):
    os.chdir(imgdir)
    sc = gen_config_file_name('scamp.conf')
    img_list = [f for f in sorted(os.listdir(imgdir)) if f.endswith('.cat')]
    img_list = [os.path.join(imgdir, f) for f in img_list]
    img_list = ','.join(img_list)
    if distortdeg:
        command = ('scamp %s -c %s -DISTORT_DEGREES %s' % (img_list, sc, distortdeg))
        print('Executing command: %s' % command)
    else:
        command = ('scamp %s -c %s' % (img_list, sc))
        print('Executing command: %s' % command)
    subprocess.run(command.split(), check=True)
    # print(pre + ext + ' scamped!')


def missfits(imgdir):
    os.chdir(imgdir)
    mc = gen_config_file_name('default.missfits')
    img_list = [f for f in sorted(os.listdir(imgdir)) if f.endswith('flat.fits')]
    img_list = [os.path.join(imgdir, f) for f in img_list]
    img_list = ' '.join(img_list)
    command = ('missfits -c %s %s' % (mc, img_list))
    print('Executing command: %s' % command)
    subprocess.run(command.split(), check=True)


#%%

def rename_head(directory):
    fnames = ['.head']
    for f in os.listdir(directory):
        for name in fnames:
            if f.endswith(name):
                f_newname = os.path.splitext(f)[0]+'.2d.head'
                path = os.path.join(directory, f)
                newpath = os.path.join(directory, f_newname)
                try:
                    os.rename(path, newpath)
                except Exception as e:
                    print(f"Error renaming file: {path} - {e}")


def double_astrom(imgdir):
    # beginning from where astrom_shift_bulge solved
    start_time = dt.now()
    print('\nSextracting shift-corrected fits files!')
    sex(imgdir, bulge=True)                     # sextract shift-solved fits files
    print('Running SCAMP w/ 2d distortion polynomial...')
    scamp(imgdir, distortdeg=2)                 # scamp shift-solved cat files w/ 2d poly solve
    print('Adding .head files directly to fits hdrs...')
    missfits(imgdir)                            # add 2d-solved scamp hdrs to shifted fits files
    print('Renaming head files to .2d.head to differentiate from next scamp run...')
    rename_head(imgdir)                         # renames .head files to .2d.head to differentiate betw. later scamp run

    print('\nSextracting 2d scamp-corrected fits files!')
    sex(imgdir, bulge=True)                     # sextract 2d-solved fits files
    print('Running SCAMP w/ 4d distortion polynomial...')
    scamp(imgdir)                               # scamp 2d-solved cat files w/ 4d poly solve
    print('Adding 4d .head files directly to fits hdrs...')
    missfits(imgdir)                            # replace 2d-solved fits file scamp hdrs w/ 4d soln scamp hdrs
    end_time = dt.now()
    print('\ndouble scamp astrometry time:', (end_time - start_time).total_seconds())


#%%


def astrometry(path, run_sex=False, run_scamp=False, run_miss=False, double_solve=False):
    if run_sex:
        sex(path)
    elif run_scamp:
        scamp(path)
    elif run_miss:
        missfits(path)
    elif double_solve:
        double_astrom(path)
    else:
        start_time = dt.now()
        sex(path)
        scamp(path)
        missfits(path)
        end_time = dt.now()
        print('\nnormal astrometry time:', (end_time - start_time).total_seconds())


def main():
    parser = argparse.ArgumentParser(description='runs sextractor and scamp on input imgs; generates LDAC .cat and .head files')
    parser.add_argument('-sex', action='store_true', help='if you want to run JUST sextractor')
    parser.add_argument('-scamp', action='store_true', help='if you want to run JUST scamp')
    parser.add_argument('-missfits', action='store_true', help='if you want to run JUST missfits')
    parser.add_argument('-double_solve', action='store_true', help='uns a 2nd order'
                                                            ' then a 4th order poly scamp solution for max accuracy')
    parser.add_argument('-path', type=str, help='[str] Images path (currently just dumps .cat '
                                                                       '& .head files in same path)')
    args, unknown = parser.parse_known_args()

    astrometry(args.path, args.sex, args.scamp, args.missfits, args.double_solve)


if __name__ == "__main__":
    main()
