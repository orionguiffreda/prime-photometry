"""
List of sky related functions to use in main pipeline
"""
import os

from photometrus.settings import gen_pipeline_file_name, gen_mflat_file_name, gen_sflat_file_name
from photometrus.sky import gen_sky
from photometrus.sky import sky_sub

from photometrus.utils.defaults import PROCESSING_DEFAULTS as defaults


def sky(FFdir, skydir, sigma, chip, **kwargs):
    # TODO change to support lists
    os.chdir(gen_pipeline_file_name())
    check_for_sky = False
    filelist = None

    if check_for_sky:
        filelist = [f for f in os.listdir(skydir) if f.endswith('.C{}.fits'.format(chip))]
        if filelist:
            print('Previous sky %s found! Skipping sky gen..\n' % filelist[0])
            pass
    if not filelist:
        print('generating sky...')

        print('\nEquivalent argparse cmd: photometrus process gensky -in_path %s -sky_path %s -sigma %s ' % (FFdir, skydir, sigma))

        gen_sky.sky_gen(in_path=FFdir, sky_path=skydir, sigma=sigma)


def skysub(FFdir, subdir, chip, skydir=None, sky_override_path=None, sex=False, bulge=False, **kwargs):
    os.chdir(gen_pipeline_file_name())
    if sky_override_path:
        skyfilepath = sky_override_path
    else:
        if not sex:
            for file in os.listdir(skydir):
                if file.endswith('.C{}.fits'.format(chip)):
                    skyfile = file
                    skyfilepath = os.path.join(skydir, skyfile)
    print('cropping and subtracting sky...')

    if sex or bulge:
        print('\nEquivalent argparse cmd: photometrus process skysub -sex -in_path %s -out_path %s' %
              (FFdir, subdir))

        sky_sub.sky_sub(in_path=FFdir, out_path=subdir, sex=True)
    else:
        print('\nEquivalent argparse cmd: photometrus process skysub -in_path %s -out_path %s -sky_path %s' %
              (FFdir, subdir, skyfilepath))

        sky_sub.sky_sub(in_path=FFdir, out_path=subdir, sky_path=skyfilepath)


SKY_FCTNS = {
    'sky': sky,
    'skysub': skysub,
}


def get_sky_commands(sex=False, sky_override=None, bulge=False, commands=None):
    """Delivers correct command logic to the main sky_functions function"""
    if commands is None:
        commands = list(defaults['default_sky_proc'])
    if sex or sky_override or bulge:
        commands = [c for c in commands if c != 'sky']
    return commands


def sky_fctns(FFdir, commands=None, **kwargs):
    """
    Runs all sky related function steps defined in `commands` (defaults to default_sky_proc).
    kwargs is a superset of everything any step might need, each step function only pulls out what named vars it needs.
    Anything else falls into its own **kwargs and is ignored.
    """
    commands = get_sky_commands(
        sex=kwargs.get('sex', False),
        sky_override=kwargs.get('sky_override', None),
        bulge=kwargs.get('bulge', False),
        commands=commands,
    )
    for cmd in commands:
        SKY_FCTNS[cmd](FFdir, **kwargs)