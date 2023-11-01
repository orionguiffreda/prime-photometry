"""
Runs improved astrometry on a given file using sextractor and scamp
"""
#%%
import os
import argparse
import subprocess

def sex(imgdir, outdir):
    os.chdir(str(outdir))
    if imgdir.endswith('.fits') or imgdir.endswith('.new'):
        try:
            pre = os.path.splitext(imgdir)[0]
            ext = os.path.splitext(imgdir)[1]
            imgname = str('0'+pre.split('/0')[1])    #don't have directory named with 0 lol
            com = ["sex ",pre + ext, ' -c /Users/orion/Desktop/PRIME/sex.config',
                   " -CATALOG_NAME " + imgname + '.sex'+ ext]  # change url for config file
            s0 = ''
            com = s0.join(com)
            out = subprocess.Popen([com], shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            out.wait()
            print(pre + ext + ' sextracted!')
        except subprocess.CalledProcessError as err:
            print('Could not run sextractor with exit error %s' % err)
#%%
def scamp(imgdir, outdir):
    os.chdir(str(outdir))
    if imgdir.endswith('.fits'): # loop to generate scamp .head headers
        pre = os.path.splitext(imgdir)[0]
        ext = os.path.splitext(imgdir)[1]
        imgname = str('0' + pre.split('/0')[1])
        com = ["scamp ", imgname + ext, ' -c /Users/orion/Desktop/PRIME/scamp.conf']
        s0 = ''
        com = s0.join(com)
        out = subprocess.Popen([com], shell=True)
        out.wait()
        print(pre + ext + ' scamped!')

#%%
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='runs sextractor on input imgs, runs scamp ; generates .head and fits LDAC files')
    parser.add_argument('files', nargs='*', type=str, metavar='f', help='Put input file first and output path last')
    args = parser.parse_args()
    sex(args.files[0],args.files[1])
    scamp(args.files[0],args.files[1])