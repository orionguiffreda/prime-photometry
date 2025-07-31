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
    with fits.open(fits_file) as fin:
        fin[0].header.update(hdr)
    if remove_header_file:
        os.remove(header_file)


def combine_header_and_fits_list(file_list, remove_header_file=False, fits_extension='.new', header_extension='.head'):
    prefixes = [os.path.splitext(f)[0] for f in file_list]
    prefixes = set(prefixes)
    for prefix in prefixes:
        header = prefix + header_extension
        fits_file = prefix + fits_extension
        if os.path.isfile(fits_file) and os.path.isfile(header):
            combine_header_and_fits(header, fits_file, remove_header_file=remove_header_file)
