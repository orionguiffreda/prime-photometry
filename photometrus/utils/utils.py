import os

from astropy.io import fits


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


def combine_header_and_fits(header_file, fits_file, remove_header_file=False):
    hdr = parse_header_file(header_file)
    with fits.open(fits_file, mode='update') as fin:
        fin[0].header.update(hdr)
    if remove_header_file:
        os.remove(header_file)


def combine_header_and_fits_list(fits_file_paths, remove_header_file=False, header_extension='.head'):
    header_file_paths = [os.path.splitext(f)[0]+header_extension for f in fits_file_paths]
    for header_file_path, fits_file_path in zip(header_file_paths, fits_file_paths):
        combine_header_and_fits(header_file_path, fits_file_path, remove_header_file)
