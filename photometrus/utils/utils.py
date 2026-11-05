import os
import re
import json5
import json
import argparse
import tempfile
from astropy.io import fits
import numpy as np

from photometrus.utils.defaults import PROCESSING_DEFAULTS as defaults



def remove_wcs_headers(header):
    wcs_keywords = [
        'WCSAXES', 'CRPIX*', 'CRVAL*', 'CDELT*', 'CUNIT*', 'CTYPE*',
        'CRDER*', 'CSYSER*', 'CD[1-9]_[1-9]*', 'PC[1-9]_[1-9]*',
        'PV*', 'PS*', 'PZ*', 'WCSPER*', 'SIP*', 'A_*', 'B_*',
        'C_*', 'D_*', 'ONAXIS*'
    ]  # I removed FOC* since it was removing FOCUS. It seems unlikely that FOC* headers will be used at any point
    for pattern in wcs_keywords:
        try:
            del header[pattern]
        except KeyError:
            continue
    return header


def parse_header_file(header_file, remove_comments=True, remove_history=True):
    fix_headers = ['CTYPE1', 'CTYPE2']
    with open(header_file, "r") as fin:
        lines = [line.strip() for line in fin]
    # Don't include END (or later lines)
    end = lines.index('END') if 'END' in lines else len(lines)
    lines = lines[:end]
    # Later pyfits versions changed this to a class method, so you can write
    # pyfits.Card.fromstring(text).  But in older pyfits versions, it was
    # a regular method.  This syntax should work in both cases.
    cards = [fits.Card().fromstring(line) for line in lines]
    header = fits.Header(cards)
    for card in fix_headers:
        header[card] = header[card].replace('TAN', 'TPV')
    if remove_comments:
        del header['COMMENT']
    if remove_history:
        del header['HISTORY']
    return header


def combine_header_and_fits(header_file, fits_file, remove_header_file=False, remove_fits_wcs=True):
    hdr = parse_header_file(header_file)
    with fits.open(fits_file, mode='update') as fin:
        fits_hdr = fin[0].header
        if remove_fits_wcs:
            fits_hdr = remove_wcs_headers(fits_hdr)
        fits_hdr.update(hdr)
        fin[0].header = fits_hdr
    if remove_header_file:
        os.remove(header_file)


def combine_header_and_fits_list(fits_file_paths, remove_header_file=False, header_extension='.head'):
    header_file_paths = [os.path.splitext(f)[0]+header_extension for f in fits_file_paths]
    for header_file_path, fits_file_path in zip(header_file_paths, fits_file_paths):
        combine_header_and_fits(header_file_path, fits_file_path, remove_header_file)


def convert_hdu_to_ldac(hdu):
    """
    Convert an hdu table to a fits_ldac table (format used by astromatic suite)

    Parameters
    ----------
    hdu: `astropy.io.fits.BinTableHDU` or `astropy.io.fits.TableHDU`
        HDUList to convert to fits_ldac HDUList

    Returns
    -------
    tbl1: `astropy.io.fits.BinTableHDU`
        Header info for fits table (LDAC_IMHEAD)
    tbl2: `astropy.io.fits.BinTableHDU`
        Data table (LDAC_OBJECTS)
    """
    tblhdr = np.array([hdu.header.tostring(',')])
    col1 = fits.Column(name='Field Header Card', array=tblhdr, format='13200A')
    cols = fits.ColDefs([col1])
    tbl1 = fits.BinTableHDU.from_columns(cols)
    tbl1.header['TDIM1'] = '(80, {0})'.format(len(hdu.header))
    tbl1.header['EXTNAME'] = 'LDAC_IMHEAD'
    tbl2 = fits.BinTableHDU(hdu.data)
    tbl2.header['EXTNAME'] = 'LDAC_OBJECTS'
    return (tbl1, tbl2)


def convert_table_to_ldac(tbl):
    """
    Convert an astropy table to a fits_ldac

    Parameters
    ----------
    tbl: `astropy.table.Table`
        Table to convert to ldac format
    Returns
    -------
    hdulist: `astropy.io.fits.HDUList`
        FITS_LDAC hdulist that can be read by astromatic software
    """
    for col in tbl.colnames:
        if tbl[col].dtype == np.float64:
            tbl[col] = tbl[col].astype(np.float32)

    f = tempfile.NamedTemporaryFile(suffix='.fits', mode='rb+')
    tbl.write(f, format='fits')
    f.seek(0)
    hdulist = fits.open(f, mode='update')
    tbl1, tbl2 = convert_hdu_to_ldac(hdulist[1])
    new_hdulist = [hdulist[0], tbl1, tbl2]
    new_hdulist = fits.HDUList(new_hdulist)
    return new_hdulist


def parentcreation(target, date, band):
    """
    Autogeneration of parent directories for processing

    Parameters
    ----------
    target: str
        target field, objname in log, ex. "field1234"
    date: str
        date of observation, in yyyymmdd format
    band: str
        filter, ex. "J"
    """

    field_dir_name = '%s_%s' % (target, date)
    field_dir = os.path.join(defaults['parent'], field_dir_name)
    parent_dir = os.path.join(field_dir, band)
    print('\nDefault parent dir: %s' % parent_dir)
    if not os.path.isdir(parent_dir):
        print('Generating parent directory!\n')
        os.makedirs(parent_dir)
    else:
        pass
    return parent_dir, field_dir


def dump_args_to_json(args, unknown=None, output_path=None):
    """
    Dumps input argparse arguments to a json5 file

    Parameters
    ----------
    args
        args that come out of : args, unknown = parser.parse_known_args()
    unknown
        unknown args that come out of above
    output_path : str
        Optionally specify output path to json dump
    """

    args_dict = {k: (str(v) if not isinstance(v, (str, int, float, bool, list, dict, type(None))) else v)
                 for k, v in vars(args).items()}
    # args_dict.pop('output_cmd_file', None)
    args_dict['unknown_args'] = unknown

    if output_path != defaults['parent']:
        json_dump_path = os.path.join(output_path, 'cmd_args.json5')
    else:
        json_dump_path = 'cmd_args.json5'
    with open(json_dump_path, 'w') as f:
        json5.dump(args_dict, f, indent=2)
    print(f'Wrote argparse args to {json_dump_path}')


def load_json_args_file(path):
    """
    Reads in args from input json filepath

    Parameters
    ----------
    path: str
        input filepath for json file of args
    """

    print(f'Input json file specified: {path}')
    with open(path, "r") as f:
        if json5 is not None:
            return json5.load(f)
        return json.load(f)


def parse_args_with_cmd_file(parser):
    """
    Checks for input json cmd file, then loads in args from there as defaults (still overwritten by any CLI args)

    Parameters
    ----------
    parser: Argparse parser
        Initial argparse parser object coming from script
    """

    # Check for input-cmd-file, examine that
    peek_parser = argparse.ArgumentParser(add_help=False)
    peek_parser.add_argument("-input_cmd_file", type=str, default=None)
    peek_args, _ = peek_parser.parse_known_args()

    unknown_from_file = {}

    if peek_args.input_cmd_file:
        file_values = load_json_args_file(peek_args.input_cmd_file)

        # fills keys for real defined args, and handles unknowns
        valid_dests = {action.dest for action in parser._actions}
        known_values = {k: v for k, v in file_values.items() if k in valid_dests}
        unknown_from_file = {k: v for k, v in file_values.items() if k not in valid_dests}

        # File vals become default
        parser.set_defaults(**known_values)

    args, unknown_from_cli = parser.parse_known_args()
    unknown = {"cli": unknown_from_cli, "file": unknown_from_file}
    return args, unknown


def handle_json_input_output(parser, target, date, band):
    """
    Argparse logic for handling json input files and output files

    Parameters
    ----------
    parser: Argparse parser
        Initial argparse parser object coming from script
    target: str
        target field, objname in log, ex. "field1234"
    date: str
        date of observation, in yyyymmdd format
    band: str
        filter, ex. "J"
    """

    args, unknown = parse_args_with_cmd_file(parser=parser)

    if any(var is None for var in (target, date, band)):
        raise ValueError('-target, -date, or -band not specified!  These must be specified to run processing!')

    if args.output_cmd_file:
        if args.parent != defaults['parent']:
            json_outdir = args.parent
        else:
            json_outdir, _ = parentcreation(args.target, args.date, args.band)
        dump_args_to_json(args, unknown, output_path=json_outdir)

    return args, unknown
