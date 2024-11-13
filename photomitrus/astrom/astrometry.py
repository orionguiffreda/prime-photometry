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


def sex(imgdir):
    os.chdir(str(imgdir))
    sx = gen_config_file_name('sex.config')
    ap = gen_config_file_name('astrom.param')
    threads = []
    print('Sextracting all images in %s' % imgdir)
    for f in sorted(os.listdir(imgdir)):
        if f.endswith('.fits') or f.endswith('.new'):
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


def scamp(imgdir):
    os.chdir(imgdir)
    sc = gen_config_file_name('scamp.conf')
    img_list = [f for f in sorted(os.listdir(imgdir)) if f.endswith('.cat')]
    img_list = [os.path.join(imgdir, f) for f in img_list]
    img_list = ','.join(img_list)
    com = ["scamp ", img_list, '-c ' + sc]
    com = ' '.join(com)
    out = subprocess.Popen([com], shell=True)
    out.wait()
    # print(pre + ext + ' scamped!')


#%%


def astrometry(path, run_sex=False, run_scamp=False):
    if run_sex:
        sex(path)
    elif run_scamp:
        scamp(path)
    else:
        sex(path)
        scamp(path)


def main():
    parser = argparse.ArgumentParser(description='runs sextractor and scamp on input imgs; generates LDAC .cat and .head files')
    parser.add_argument('-sex', action='store_true', help='if you want to run JUST sextractor')
    parser.add_argument('-scamp', action='store_true', help='if you want to run JUST scamp')
    parser.add_argument('-path', type=str, help='[str] Images path (currently just dumps .cat '
                                                                       '& .head files in same path)')
    args, unknown = parser.parse_known_args()

    astrometry(args.path, args.sex, args.scamp)


if __name__ == "__main__":
    main()
