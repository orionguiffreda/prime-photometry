"""
Runs improved astrometry on a given file using sextractor and scamp
"""
#%%
import os
import argparse

def sex(imgdir, outdir):
    os.chdir(str(outdir))
    for f in os.listdir(str(imgdir)):
        if f.endswith('.fits') or f.endswith('.new'):
            pre = os.path.splitext(f)[0]
            ext = os.path.splitext(f)[1]
            com = ["sex ", imgdir + pre + ext, ' -c /Users/orion/Desktop/PRIME/sex.config',
                   " -CATALOG_NAME " + pre + '.fits']  # change url for config file
            s0 = ''
            com = s0.join(com)
            out = subprocess.Popen([com], shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            out.wait()
            print(pre + ext + ' sextracted!')

def scamp(imgdir, outdir):
    for f in os.listdir(str(outdir)):  # loop to generate scamp .head headers
        if f.endswith('.fits'):
            pre = os.path.splitext(f)[0]
            ext = os.path.splitext(f)[1]
            com = ["scamp ", outdir + pre + ext, ' -c /Users/orion/Desktop/PRIME/scamp.conf']
            s0 = ''
            com = s0.join(com)
            out = subprocess.Popen([com], shell=True,stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            out.wait()
            print(pre + ext + ' scamped!')
    for f in os.listdir(str(outdir)):  # loop to move .heads to img dir, so both are in same dir for swarp to stack, remove if you don't want that
        if f.endswith('.head'):
            pre1 = os.path.splitext(f)[0]
            ext1 = os.path.splitext(f)[1]
            os.rename(outdir + pre1 + ext1, imgdir + pre1 + ext1)
            print(pre1 + ext1 + ' moved!')

#%%
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='runs sextractor on input imgs, runs scamp, then moves scamp ext. hdrs '
                                                 'to img dir.; generates .head and fits LDAC files')
    parser.add_argument('files', nargs='*', type=str, metavar='f', help='Put input file first and output path last')
    args = parser.parse_args()
    sex(args.files[0],args.files[1])
    scamp(args.files[0],args.files[1])