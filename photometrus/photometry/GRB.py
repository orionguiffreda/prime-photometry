"""
Functions for GRB observations to be used w/ photometry.py, including main grb catalog analysis and new source search
(for radii > 1')
"""
import os
import re

import astropy.nddata.utils
import astropy.units as u
from astropy.coordinates import Angle, SkyCoord
from astropy.wcs import WCS
from astropy.wcs import utils
from astropy.stats import sigma_clipped_stats
from astropy.io import fits
from astropy.io import ascii
from astropy.nddata import Cutout2D
from astropy.table import Column
from astropy.table import Table

import base64
from bs4 import BeautifulSoup
import numpy as np
from regions import CircleSkyRegion
import matplotlib.pyplot as plt

from photometrus.settings import (PHOTOMETRY_MAG_LOWER_LIMIT, PHOTOMETRY_LIM_MAGS, MAGTYPES)

from photometrus.utils.defaults import PROCESSING_DEFAULTS as defaults


def GRB(ra, dec, imageName, survey, band, thresh, massCatCoords, good_cat_stars, directory, chip, coordlist=None,
        grbname=defaults['grb_name'], **kwargs):
    """
    Main GRB functionality.  Searches through photometry .ecsv catalog and crossmatchesv to see if a source was detected
    in a specified area.  If found, crossmatches to input survey to determine if the source is new, then writes many
    statistics to a seperate .ecsv.  Also outputs .reg files for query threshold and sources in this threshold, along
    with a summary .html file.

    Parameters
    ----------
    ra: float
        RA coordinate for search
    dec: float
        Dec coordinate for search
    imageName: str
        Filename of input image
    band: str
        Filter observation was taken in
    thresh: float
        Diameter of threshold to search for sources
    massCatCoords: SkyCoord
        Coordinate data of queried survey sources
    good_cat_stars: Table
        Astropy Table of AB-corrected survey sources (incl. ra, dec, mag, & mag error columns)
    directory: str
        Directory the input image is stored in
    chip: int
        Detector number of input image
    coordlist: list
        Optional, list of multiple different locations to use GRB search on within 1 image
    grbname: str
        Optional, custom name for all GRB data products (default = 'GRB')
    """

    mag_ecsvname = '%s.%s.ecsv' % (imageName, survey)
    mag_ecsvtable = ascii.read(mag_ecsvname)
    mag_ecsvcleanSources = mag_ecsvtable
    mag_ecsvsourceCatCoords = SkyCoord(ra=mag_ecsvcleanSources['ALPHA_J2000'], dec=mag_ecsvcleanSources['DELTA_J2000'],
                                       frame='icrs',
                                       unit='degree')

    if not imageName.startswith(f'coadd.Open-{band}') and not imageName.endswith('.fits'):
        num = 'img'
    else:
        num = imageName[-16:-8]

    # limiting mag for survey
    try:
        survey_lim_mag = PHOTOMETRY_LIM_MAGS[survey]
    except KeyError:
        survey_lim_mag = 19.7
        print(f'Error in finding lim mag for survey catalogs, no matching catalog found? '
              f'Using generous PRIME lim: {survey_lim_mag}')

    # postage stamp cutout fctn (png & fits)
    def grb_cutout(imageName, GRBcoords, photoDistThresh, loc=None, append=False, regprimename=None, regsurvname=None):
        imgdata = fits.getdata(imageName)
        img = fits.open(imageName)
        head = img[0].header
        w = WCS(head)

        if loc:
            savename = '%s_%s_C%i_Cutout_%s_%s_loc_%d' % (grbname, band, chip, survey, num, loc)
            threshname = '%s_queries_thresh.reg' % grbname
        else:
            savename = '%s_%s_C%i_Cutout_%s_%s' % (grbname, band, chip, survey, num)
            threshname = '%s_query_thresh.reg' % grbname

        size = 4 * photoDistThresh * u.arcsec
        try:
            cutout = Cutout2D(imgdata, GRBcoords[0], size, wcs=w, copy=True)
            region = CircleSkyRegion(center=GRBcoords[0], radius=Angle(thresh, unit='arcsec'))
            pix_region = region.to_pixel(cutout.wcs)
        except astropy.nddata.utils.NoOverlapError:
            print(' Area of GRB threshold not found within image, cannot generate cutout!')
            return savename, threshname

        mean, median, sigma_cut = sigma_clipped_stats(cutout.data)
        plt.figure(10, figsize=(8, 8))
        plt.imshow(cutout.data, vmin=median - 3 * sigma_cut, vmax=median + 3 * sigma_cut, origin='lower',
                   cmap='viridis')
        pix_region.plot(color='cyan', ls='--', label='Input GRB threshold')

        if regprimename:
            primeregs = open('%s_PRIME_srcs.reg' % grbname, 'r')
            plt_primeregs = []
            primeallregs = [reg for reg in primeregs if reg != 'fk5\n']
            for reg in primeallregs:
                nums = re.findall(r'[-+]?\d*\.?\d+', reg)
                srcra = float(nums[0])
                srcdec = float(nums[1])
                srcrad = float(nums[2])
                srccoords = SkyCoord(ra=[srcra], dec=[srcdec], frame='icrs', unit='degree')
                srcreg = CircleSkyRegion(center=srccoords[0], radius=Angle(srcrad, unit='arcsec'))
                plt_primeregs.append(srcreg)

            for reg in plt_primeregs:
                pix_reg = reg.to_pixel(cutout.wcs)
                pix_reg.plot(color='red', ls='-', label='PRIME Source')

        if regsurvname:
            survregs = open('%s_%s_srcs.reg' % (grbname, survey), 'r')
            plt_survregs = []
            survallregs = [reg for reg in survregs if reg != 'fk5\n']
            for reg in survallregs:
                nums = re.findall(r'[-+]?\d*\.?\d+', reg)
                srcra = float(nums[0])
                srcdec = float(nums[1])
                srcrad = float(nums[2])
                srccoords = SkyCoord(ra=[srcra], dec=[srcdec], frame='icrs', unit='degree')
                srcreg = CircleSkyRegion(center=srccoords[0], radius=Angle(srcrad, unit='arcsec'))
                plt_survregs.append(srcreg)

            for reg in plt_survregs:
                pix_reg = reg.to_pixel(cutout.wcs)
                pix_reg.plot(color='magenta', ls='-', label='%s Source' % survey)

        handles, labels = plt.gca().get_legend_handles_labels()
        seen = set()
        filtered_handles = []
        filtered_labels = []
        for h, l in zip(handles, labels):
            if l not in seen:
                filtered_handles.append(h)
                filtered_labels.append(l)
                seen.add(l)
        plt.legend(filtered_handles, filtered_labels, loc='best')
        plt.savefig(savename + '.png', dpi=300)
        plt.close()

        fits.writeto(savename + '.fits', cutout.data, cutout.wcs.to_header(), overwrite=True)

        # ds9 regions
        if append:
            newtext = open(threshname, 'a')  # input threshold
            newtext.write(f'\ncircle({ra}, {dec}, {photoDistThresh}") # color=cyan width=2 text={{Query Thresh}}')
        else:
            newtext = open(threshname, 'w+')  # input threshold
            newtext.write('fk5')
            newtext.write(f'\ncircle({ra}, {dec}, {photoDistThresh}") # color=cyan width=2 text={{Query Thresh}}')

        print(' Exported cutout & ds9 region of GRB search area!')

        return savename, threshname

    # source ds9 region writing
    def source_reg_gen(src_ra=0, src_dec=0, rad=2, src_survey=None, append=False):
        if not src_survey:
            name = '%s_PRIME_srcs.reg' % grbname
            color = 'red'
        else:
            name = '%s_%s_srcs.reg' % (grbname, survey)
            color = 'yellow'
        if not append:
            newtext = open(name, 'w+')
            newtext.write('fk5')
            if src_ra != 0 and src_dec != 0:
                newtext.write(f'\ncircle({src_ra}, {src_dec}, {rad}") # color={color}')
        else:
            newtext = open(name, 'a')
            newtext.write(f'\ncircle({src_ra}, {src_dec}, {rad}") # color={color}')

        return name

    # mag diff calc betw. survey and prime for existing source crsmtches
    def mag_diff_calc(survey_cat, prime_cat, survey_idx, prime_idx, d2d, band):
        """calculates diff betw. crossmatched prime source and survey, generating seperation & appropriate flags"""
        # prime & survey mag cols
        all_magtypes = set(MAGTYPES.keys())
        prime_mag_cols = sorted([col for col in prime_cat.colnames if any(mag in col for mag in all_magtypes)
                                 and band in col and f'{band}MAG' in col and 'e_' not in col])
        prime_err_cols = sorted([col for col in prime_cat.colnames if any(mag in col for mag in all_magtypes)
                                 and band in col and f'e_{band}MAG' in col])

        colnames = survey_cat.colnames

        # mag diff calc for all applicable cols
        mag_diff_ar = []
        flag_ar = []

        for prime_mags, prime_errs in zip(prime_mag_cols, prime_err_cols):
            col_magtype = prime_mags.split('_')[-1]

            if len(colnames) > 6:
                magcolname = f'{band}MAG_{col_magtype}'
                magerrcolname = f'e_{band}MAG_{col_magtype}'
            else:
                magcolname = colnames[2]
                magerrcolname = colnames[3]

            if len(prime_idx) > 1:
                if 0 in prime_idx:
                    if not np.isscalar(prime_cat[f'{band}MAG_{col_magtype}']):
                        mag_diff = float(survey_cat[magcolname][survey_idx][0]) - float(
                            prime_cat[f'{band}MAG_{col_magtype}'][prime_idx][0])

                        comb_err = np.sqrt(survey_cat[magerrcolname][survey_idx][0] ** 2 +
                                           prime_cat[f'e_{band}MAG_{col_magtype}'][prime_idx][0] ** 2)
                    else:
                        mag_diff = float(survey_cat[magcolname][survey_idx][0]) - float(prime_cat[f'{band}MAG_{col_magtype}'])

                        comb_err = np.sqrt(survey_cat[magerrcolname][survey_idx][0] ** 2 +
                                           prime_cat[f'e_{band}MAG_{col_magtype}'] ** 2)

                    sep = d2d[0]
                else:
                    mag_diff = float(survey_cat[magcolname][survey_idx]) - float(
                        prime_cat[f'{band}MAG_{col_magtype}'][prime_idx])
                    # errors in quad
                    comb_err = np.sqrt(survey_cat[magerrcolname][survey_idx] ** 2 +
                                       prime_cat[f'e_{band}MAG_{col_magtype}'][prime_idx] ** 2)
                    sep = d2d
            else:
                # survey mag - prime mag
                mag_diff = float(survey_cat[magcolname][survey_idx].item()) - float(prime_cat[f'{band}MAG_{col_magtype}'].item())

                # errors in quad
                comb_err = np.sqrt(survey_cat[magerrcolname][survey_idx] ** 2 +
                                   prime_cat[f'e_{band}MAG_{col_magtype}'] ** 2)

                sep = d2d

            err_3_sig = 3 * comb_err

            if abs(mag_diff) > err_3_sig:
                flg = 2
                # print(' At least 1 source in threshold has a significant mag difference to the catalog!')
            else:
                flg = 1

            if np.size(sep) > 1:
                sep = sep[0]

            sep_dist = sep.to(u.arcsec)
            sep_dist = sep_dist / u.arcsec

            mag_diff_ar.append((mag_diff, col_magtype))
            # mag_diff_ar.append(mag_diff)
            flag_ar.append((np.int16(flg), col_magtype))

        return mag_diff_ar, float(sep_dist), flag_ar

    # mag diff calc betw. survey and prime for existing source crsmtches
    def non_cm_flag_calc(prime_source, band, survey_lim_mag):
        """function to determine if non-crossmatched source is deeper than limiting mag of survey, flag = 3 if so, 0 if not"""
        # prime & survey mag cols
        all_magtypes = set(MAGTYPES.keys())
        prime_mag_cols = sorted([col for col in prime_source.colnames if any(mag in col for mag in all_magtypes)
                                 and band in col and f'{band}MAG' in col and 'e_' not in col])

        # array of flags for each mag col
        flag_ar = []

        for mag_col in prime_mag_cols:
            col_magtype = mag_col.split('_')[-1]

            src_mag = prime_source[mag_col]
            if src_mag > survey_lim_mag:
                flg = 3
            else:
                flg = 0
            flag_ar.append((np.int16(flg), col_magtype))

        return flag_ar

    def mag_and_flag_grabber(mag_diff_ar, flag_ar, output, col_magtype='AUTO'):
        """grabbing correct mag diff & flag values from input arrays"""
        if isinstance(mag_diff_ar[0], tuple):
            # if input mag diff arr is real, return appropriate mag diff and flag vals
            mag_diff_val = [t[0] for t in mag_diff_ar if t[1] == col_magtype]
            flag_val = [t[0] for t in flag_ar if t[1] == col_magtype]
            if output == 'mag':
                return mag_diff_val[0]
            elif output == 'flag':
                return flag_val[0]
            else:
                raise Exception('output option must be either "mag" or "flag"!')
        elif isinstance(flag_ar[0], tuple):
            # if only input flag arr is real, return appropriate flag val
            flag_val = [t[0] for t in flag_ar if t[1] == col_magtype]
            if output == 'mag':
                return 99
            if output == 'flag':
                return flag_val[0]
            else:
                raise Exception('output option must be either "mag" or "flag"!')
        else:
            if output == 'mag':
                return 99
            elif output == 'flag':
                return np.int16(0)
            else:
                raise Exception('output option must be either "mag" or "flag"!')

    # generation of html file
    def html_gen(data, directory, savename, threshname, band, survey, ra, dec, thresh, survname=None, primename=None):
        # building table
        newtbl = data.copy()
        # newtbl.remove_columns([f'{band}apMag', f'{band}apMag_Err'])

        df = newtbl.to_pandas()

        descriptions = {}
        for col in newtbl.colnames:
            desc = newtbl[col].description
            if desc is None or str(desc).strip() == "":
                desc = f"Description of {col}"
            descriptions[col] = desc

        init_html = df.to_html(index=False)

        soup = BeautifulSoup(init_html, "html.parser")

        for th in soup.find_all("th"):
            col_name = th.text.strip()
            if col_name in descriptions:
                th['title'] = descriptions[col_name]

        tbl_html = str(soup)

        # adding img & command
        with open(savename + '.png', "rb") as img_file:
            encoded = base64.b64encode(img_file.read()).decode("utf-8")

        fits_items = [savename + '.fits', threshname, survname, primename]
        fits_items = [os.path.join(directory, f) for f in fits_items if f is not None]

        base = fits_items[0]
        regions = " ".join(f"-regions {f}" for f in fits_items[1:])

        ds9_command = f"ds9 {base} {regions} &"

        # final html gen
        final_html = f"""
        <div style="text-align: center; font-family: Arial, sans-serif;">

            <h2>GRB Information</h2>
            <p style="margin-top: 0; margin-bottom: 10px; font-size: 16px; color: #555;">
                RA = {ra}, Dec = {dec}, threshold = {thresh}"
            </p>

            <img src="data:image/png;base64,{encoded}" alt="GRB Cutout Region" width="600" style="margin-bottom: 10px;">

            <p style="margin-top: 0px; margin-bottom: 10px; font-size: 16px; color: #555;">
                To see source regions, open GRB stamp in DS9 through terminal: 
            </p>
            <p style="margin-top: 0px; margin-bottom: 20px; font-size: 14px; color: #030303;">
                {ds9_command}
            </p>

            <div style="display: inline-block; text-align: left;">
                {tbl_html}
            </div>

        </div>
        """

        with open('%s_Multisource_%s_C%i_Data_%s_%s.html' % (grbname, band, chip, survey, num), "w") as f:
            f.write(final_html)

    # sexigesimal conversion
    try:
        float(ra)
        ra = ra
        dec = dec
    except ValueError:
        coords = ra + ' ' + dec
        print('Sexagesimal RA = %s & Dec = %s' % (ra, dec))

        deci_coords = SkyCoord(coords, frame='icrs', unit=(u.hourangle, u.deg)).to_string()
        deci_coords = deci_coords.split(' ')

        ra = deci_coords[0]
        dec = deci_coords[1]

    photoDistThresh = thresh
    if coordlist:
        for i in range(len(coordlist)):
            print('Checking GRB location %i: %s' % (i, coordlist[i]))
        coordlist = tuple(eval(i) for i in coordlist)
        try:
            idx_GRBpsfdict = {}
            keys = np.arange(0, len(coordlist), 1)
            for i in range(len(coordlist)):
                GRBcoords = SkyCoord(ra=[coordlist[i][0]], dec=[coordlist[i][1]], frame='icrs', unit='degree')
                idx_GRB, idx_GRBcleanpsf, d2d, d3d = mag_ecsvsourceCatCoords.search_around_sky(GRBcoords,
                                                                                               photoDistThresh * u.arcsec)

                # survey source crsmtch
                idx_survey, idx_surveycleanpsf, d2d_surv, d3d_surv = massCatCoords.search_around_sky(GRBcoords,
                                                                                                     photoDistThresh * u.arcsec)
                if len(idx_surveycleanpsf) > 0:
                    print(' %i %s existing sources found within GRB threshold! Writing to DS9 reg files...'
                          % (len(idx_surveycleanpsf), survey))
                    source_reg_gen(src_survey=survey)
                    for idx in idx_surveycleanpsf:
                        src = massCatCoords[idx]
                        source_reg_gen(src.ra.deg, src.dec.deg, src_survey=survey, append=True)
                else:
                    print(' No existing %s sources found within GRB threshold!' % survey)

                if i == 0:
                    grb_cutout(imageName, GRBcoords, photoDistThresh, loc=i)
                else:
                    grb_cutout(imageName, GRBcoords, photoDistThresh, loc=i, append=True)

                idx_GRBpsfdict[keys[i]] = idx_GRBcleanpsf, coordlist[i]
        except NameError:
            print('No Sources found!')

    else:
        print('Checking GRB location %s, %s' % (ra, dec))
        GRBcoords = SkyCoord(ra=[ra], dec=[dec], frame='icrs', unit='degree')
        # prime source crsmtch
        idx_GRB, idx_GRBcleanpsf, d2d, d3d = mag_ecsvsourceCatCoords.search_around_sky(GRBcoords,
                                                                                       photoDistThresh * u.arcsec)
        # survey source crsmtch
        idx_survey, idx_surveycleanpsf, d2d_surv, d3d_surv = massCatCoords.search_around_sky(GRBcoords,
                                                                                             photoDistThresh * u.arcsec)
        if len(idx_surveycleanpsf) > 0:
            print(' %i %s existing sources found within GRB threshold! Writing to DS9 reg files...'
                  % (len(idx_surveycleanpsf), survey))
            source_reg_gen(src_survey=survey)
            for idx in idx_surveycleanpsf:
                src = massCatCoords[idx]
                regsurvname = source_reg_gen(src.ra.deg, src.dec.deg, src_survey=survey, append=True)
        else:
            print(' No existing %s sources found within GRB threshold!' % survey)
            regsurvname = None

        savename, threshname = grb_cutout(imageName, GRBcoords, photoDistThresh)

    # custom col descriptions
    desc = {
        "mag_aper": f'{band} band 2.5" aperture magnitude',
        "mag_aper_err": f'Error in {band} band 2.5" aperture magnitude',
        "mag_err_crsmtch": f'{survey} - PRIME source mag for survey crossmatched source, 99 if no crossmatch',
        "distance": '2D distance (arcsec) betw. PRIME source & input GRB coords',
        "separation": f'2D distance (arcsec) betw. PRIME & crossmatched {survey} source, -1 = no match',
        "crsmtch_flg": (f'Flag for {survey} crossmatch: '
                        f'0 = no match, 1 = match w/ mag diff within 3 sig, '
                        f'2 = match w/ mag diff outside 3 sig, '
                        f'3 = no match, but deeper than survey limiting mag ({survey_lim_mag})')
    }

    mag_col_num = len([col for col in mag_ecsvcleanSources.colnames if any(mag in col for mag in set(MAGTYPES.keys()))
                             and band in col and f'{band}MAG' in col and 'e_' not in col])

    def get_desc(table, colname):
        """Returns col description if col exists"""
        return table[colname].description if colname in table.colnames else None

    # grb table initiation
    grbdata = Table()

    # grb table descriptions, always present
    ra_desc = mag_ecsvcleanSources['ALPHA_J2000'].description
    dec_desc = mag_ecsvcleanSources['DELTA_J2000'].description
    rad_desc = mag_ecsvcleanSources['FLUX_RADIUS'].description
    snr_desc = mag_ecsvcleanSources['SNR_WIN'].description
    elon_desc = mag_ecsvcleanSources['ELONGATION'].description
    dist_desc = desc['distance']
    sep_desc = desc['separation']
    flag_desc = desc['crsmtch_flg']
    mag_diff_desc = desc['mag_err_crsmtch']

    # grb table descriptions, mags
    psf_mag_desc = get_desc(mag_ecsvcleanSources, '%sMAG_PSF' % band)
    psf_mag_err_desc = get_desc(mag_ecsvcleanSources, 'e_%sMAG_PSF' % band)
    aper_mag_desc = get_desc(mag_ecsvcleanSources, '%sMAG_APER' % band)
    aper_mag_err_desc = get_desc(mag_ecsvcleanSources, 'e_%sMAG_APER' % band)
    auto_mag_desc = get_desc(mag_ecsvcleanSources, '%sMAG_AUTO' % band)
    auto_mag_err_desc = get_desc(mag_ecsvcleanSources, 'e_%sMAG_AUTO' % band)

    if coordlist:
        source_reg_gen()
        for key in idx_GRBpsfdict:
            values = idx_GRBpsfdict[key]
            idx_GRBcleanpsf = values[0]
            ra = values[1][0]
            dec = values[1][1]
            print(' idx size = %d for location %d' % (len(idx_GRBcleanpsf), key))
            if len(idx_GRBcleanpsf) == 1:
                print(' GRB source at inputted coords %s and %s, rad = %s arcsec found!' % (ra, dec, photoDistThresh))
                if psf_mag_desc is not None:
                    grb_mag = mag_ecsvcleanSources[idx_GRBcleanpsf]['%sMAG_PSF' % band][0]
                    grb_magerr = mag_ecsvcleanSources[idx_GRBcleanpsf]['e_%sMAG_PSF' % band][0]
                if aper_mag_desc is not None:
                    grb_mag_aper = mag_ecsvcleanSources[idx_GRBcleanpsf]['%sMAG_APER' % band][0]
                    grb_magerr_aper = mag_ecsvcleanSources[idx_GRBcleanpsf]['e_%sMAG_APER' % band][0]
                if auto_mag_desc is not None:
                    grb_mag_auto = mag_ecsvcleanSources[idx_GRBcleanpsf]['%sMAG_AUTO' % band][0]
                    grb_magerr_auto = mag_ecsvcleanSources[idx_GRBcleanpsf]['e_%sMAG_AUTO' % band][0]

                grb_ra = mag_ecsvcleanSources[idx_GRBcleanpsf]['ALPHA_J2000'][0]
                grb_dec = mag_ecsvcleanSources[idx_GRBcleanpsf]['DELTA_J2000'][0]
                grb_rad = mag_ecsvcleanSources[idx_GRBcleanpsf]['FLUX_RADIUS'][0]
                grb_snr = mag_ecsvcleanSources[idx_GRBcleanpsf]['SNR_WIN'][0]
                grb_elon = mag_ecsvcleanSources[idx_GRBcleanpsf]['ELONGATION'][0]
                grb_dist = d2d[0].to(u.arcsec)
                grb_dist = grb_dist / u.arcsec

                all_magtypes = set(MAGTYPES.keys())
                for mag in sorted(all_magtypes):
                    try:
                        print(f' %s {mag} magnitude of GRB is %.2f +/- %.2f' % (band,
                                                                                mag_ecsvcleanSources[idx_GRBcleanpsf][
                                                                                    f'{band}MAG_{mag}'][0],
                                                                                mag_ecsvcleanSources[idx_GRBcleanpsf][
                                                                                    f'e_{band}MAG_{mag}'][0]))
                    except KeyError:
                        pass

                # survey crsmtch check
                idx_both, idx_bothcleanpsf, d2d_crs, d3d_crs = massCatCoords.search_around_sky(
                    mag_ecsvsourceCatCoords[idx_GRBcleanpsf],
                    grb_rad * u.arcsec)
                prime_crs_cat = mag_ecsvcleanSources[idx_GRBcleanpsf]
                if len(idx_bothcleanpsf) > 0:
                    mag_diff_crs_ar, sep, survey_flg_ar = mag_diff_calc(survey_cat=good_cat_stars, prime_cat=prime_crs_cat,
                                                                  survey_idx=idx_bothcleanpsf,
                                                                  prime_idx=[idx_GRBcleanpsf][idx_both],
                                                                  d2d=d2d_crs, band=band)
                else:
                    mag_diff_crs_ar = [99] * mag_col_num
                    survey_flg_ar = non_cm_flag_calc(prime_source=prime_crs_cat, band=band, survey_lim_mag=survey_lim_mag)
                    sep = -1

                grbdata['RA'] = Column(np.round(np.array([grb_ra]), 5) * u.deg, description=ra_desc)
                grbdata['DEC'] = Column(np.round(np.array([grb_dec]), decimals=5) * u.deg, description=dec_desc)

                if psf_mag_desc is not None:
                    grbdata['%spsfMag' % band] = Column(np.round(np.array([grb_mag]), 3) * u.ABmag,
                        description=psf_mag_desc
                    )
                    grbdata['%spsfMag_Err' % band] = Column(np.round(np.array([grb_magerr]), 3) * u.ABmag,
                        description=psf_mag_err_desc
                    )
                    grbdata['%spsfME_CM' % band] = Column(
                        np.round(np.array([mag_and_flag_grabber(mag_diff_crs_ar, survey_flg_ar, 'mag', 'PSF')]), 3) * u.ABmag,
                        description=mag_diff_desc
                    )
                    grbdata['psfS_CM'] = Column(mag_and_flag_grabber(mag_diff_crs_ar, survey_flg_ar, 'flag', 'PSF'),
                                                description=flag_desc)

                if aper_mag_desc is not None:
                    grbdata['%saperMag' % band] = Column(np.round(np.array([grb_mag_aper]), 3) * u.ABmag,
                        description=aper_mag_desc
                    )
                    grbdata['%saperMag_Err' % band] = Column(np.round(np.array([grb_magerr_aper]), 3) * u.ABmag,
                        description=aper_mag_err_desc
                    )
                    grbdata['%saperME_CM' % band] = Column(
                        np.round(np.array([mag_and_flag_grabber(mag_diff_crs_ar, survey_flg_ar, 'mag', 'APER')]), 3) * u.ABmag,
                        description=mag_diff_desc
                    )
                    grbdata['aperS_CM'] = Column(mag_and_flag_grabber(mag_diff_crs_ar, survey_flg_ar, 'flag', 'APER'),
                                                 description=flag_desc)

                if auto_mag_desc is not None:
                    grbdata['%sautoMag' % band] = Column(np.round(np.array([grb_mag_auto]), 3) * u.ABmag,
                        description=auto_mag_desc
                    )
                    grbdata['%sautoMag_Err' % band] = Column(np.round(np.array([grb_magerr_auto]), 3) * u.ABmag,
                        description=auto_mag_err_desc
                    )
                    grbdata['%sautoME_CM' % band] = Column(
                        np.round(np.array([mag_and_flag_grabber(mag_diff_crs_ar, survey_flg_ar, 'mag', 'AUTO')]), 3) * u.ABmag,
                        description=mag_diff_desc
                    )
                    grbdata['autoS_CM'] = Column(mag_and_flag_grabber(mag_diff_crs_ar, survey_flg_ar, 'flag', 'AUTO'),
                                                 description=flag_desc)

                grbdata['Radius'] = Column(np.round(np.array([grb_rad]), decimals=2) * u.arcsec, description=rad_desc)
                grbdata['SNR'] = Column(np.round(np.array([grb_snr]), decimals=2), description=snr_desc)
                grbdata['Elongation'] = Column(np.round(np.array([grb_elon]), decimals=3), description=elon_desc)
                grbdata['Distance'] = Column(np.round(np.array([grb_dist]), decimals=5) * u.arcsec, description=dist_desc)
                grbdata['Separation'] = Column(np.round(np.array([sep]), decimals=5) * u.arcsec, description=sep_desc)

                grbdata.write('%s_%s_C%i_Data_%s_%s_loc_%d.ecsv' % (grbname, band, chip, survey, num, key),
                              overwrite=True)

                source_reg_gen(grb_ra, grb_dec, rad=grb_rad, append=True)
                print(' Generated GRB data table & source DS9 regions!')
            elif len(idx_GRBcleanpsf) > 1:
                print(' Multiple sources detected in search radius (ra = %s, dec = %s, rad = %s arcsec)'
                      ', refer to .ecsv file for source info!' % (ra, dec, photoDistThresh))
                mag_ar = []
                mag_err_ar = []
                apmag_ar = []
                apmag_err_ar = []
                automag_ar = []
                automag_err_ar = []
                ra_ar = []
                dec_ar = []
                rad_ar = []
                snr_ar = []
                elon_ar = []
                dist_ar = []
                sep_ar = []
                diff_ar = []
                crsmtch_ar = []
                idx_GRBcleanpsflist = idx_GRBcleanpsf.tolist()
                for i in idx_GRBcleanpsflist:
                    if psf_mag_desc is not None:
                        grb_mag = mag_ecsvcleanSources[idx_GRBcleanpsf]['%sMAG_PSF' % band][idx_GRBcleanpsflist.index(i)]
                        mag_ar.append(grb_mag)
                        grb_magerr = mag_ecsvcleanSources[idx_GRBcleanpsf]['e_%sMAG_PSF' % band][
                            idx_GRBcleanpsflist.index(i)]
                        mag_err_ar.append(grb_magerr)
                    if aper_mag_desc is not None:
                        grb_mag_aper = mag_ecsvcleanSources[idx_GRBcleanpsf]['%sMAG_APER' % band][idx_GRBcleanpsflist.index(i)]
                        apmag_ar.append(grb_mag_aper)
                        grb_magerr_aper = mag_ecsvcleanSources[idx_GRBcleanpsf]['e_%sMAG_APER' % band][
                            idx_GRBcleanpsflist.index(i)]
                        apmag_err_ar.append(grb_magerr_aper)
                    if auto_mag_desc is not None:
                        grb_mag_auto = mag_ecsvcleanSources[idx_GRBcleanpsf]['%sMAG_AUTO' % band][
                            idx_GRBcleanpsflist.index(i)]
                        automag_ar.append(grb_mag_auto)
                        grb_magerr_auto = mag_ecsvcleanSources[idx_GRBcleanpsf]['e_%sMAG_AUTO' % band][
                            idx_GRBcleanpsflist.index(i)]
                        automag_err_ar.append(grb_magerr_auto)
                    grb_ra = mag_ecsvcleanSources[idx_GRBcleanpsf]['ALPHA_J2000'][idx_GRBcleanpsflist.index(i)]
                    ra_ar.append(grb_ra)
                    grb_dec = mag_ecsvcleanSources[idx_GRBcleanpsf]['DELTA_J2000'][idx_GRBcleanpsflist.index(i)]
                    dec_ar.append(grb_dec)
                    grb_rad = mag_ecsvcleanSources[idx_GRBcleanpsf]['FLUX_RADIUS'][idx_GRBcleanpsflist.index(i)]
                    rad_ar.append(grb_rad)
                    grb_snr = mag_ecsvcleanSources[idx_GRBcleanpsf]['SNR_WIN'][idx_GRBcleanpsflist.index(i)]
                    snr_ar.append(grb_snr)
                    grb_elon = mag_ecsvcleanSources[idx_GRBcleanpsf]['ELONGATION'][idx_GRBcleanpsflist.index(i)]
                    elon_ar.append(grb_elon)
                    dist = (d2d[idx_GRBcleanpsflist.index(i)]).to(u.arcsec)
                    dist = dist / u.arcsec
                    dist_ar.append(dist)

                    # survey crsmtch check
                    checkcoords = SkyCoord(ra=[mag_ecsvsourceCatCoords[i].ra.deg],
                                           dec=[mag_ecsvsourceCatCoords[i].dec.deg],
                                           frame='icrs', unit='degree')
                    idx_both, idx_bothcleanpsf, d2d_crs, d3d_crs = massCatCoords.search_around_sky(
                        checkcoords, grb_rad * u.arcsec)
                    prime_crs_cat = mag_ecsvcleanSources[i]
                    if len(idx_bothcleanpsf) > 0:
                        mag_diff_crs_ar, sep, survey_flg_ar = mag_diff_calc(survey_cat=good_cat_stars,
                                                                      prime_cat=prime_crs_cat,
                                                                      survey_idx=idx_bothcleanpsf,
                                                                      prime_idx=[i][idx_both],
                                                                      d2d=d2d_crs, band=band)
                    else:
                        mag_diff_crs_ar = [99] * mag_col_num
                        survey_flg_ar = non_cm_flag_calc(prime_source=prime_crs_cat, band=band,
                                                         survey_lim_mag=survey_lim_mag)
                        sep = -1

                    diff_ar.append(mag_diff_crs_ar)
                    crsmtch_ar.append(survey_flg_ar)
                    sep_ar.append(sep)

                    source_reg_gen(grb_ra, grb_dec, rad=grb_rad, append=True)

                grbdata['RA'] = Column(np.round(np.array(ra_ar), 5) * u.deg, description=ra_desc)
                grbdata['DEC'] = Column(np.round(np.array(dec_ar), decimals=5) * u.deg, description=dec_desc)

                if psf_mag_desc is not None:
                    psf_diff_ar = [mag_and_flag_grabber(md, sf, 'mag', 'PSF') for md, sf in zip(diff_ar, crsmtch_ar)]
                    psf_flg_ar = [mag_and_flag_grabber(md, sf, 'flag', 'PSF') for md, sf in zip(diff_ar, crsmtch_ar)]

                    grbdata['%spsfMag' % band] = Column(np.round(np.array(mag_ar), 3) * u.ABmag,
                        description=psf_mag_desc
                    )
                    grbdata['%spsfMag_Err' % band] = Column(np.round(np.array(mag_err_ar), 3) * u.ABmag,
                        description=psf_mag_err_desc
                    )
                    grbdata['%spsfME_CM' % band] = Column(np.round(np.array(psf_diff_ar), 3) * u.ABmag,
                        description=mag_diff_desc
                    )
                    grbdata['psfS_CM'] = Column(np.array(psf_flg_ar), description=flag_desc)

                if aper_mag_desc is not None:
                    aper_diff_ar = [mag_and_flag_grabber(md, sf, 'mag', 'APER') for md, sf in zip(diff_ar, crsmtch_ar)]
                    aper_flg_ar = [mag_and_flag_grabber(md, sf, 'flag', 'APER') for md, sf in zip(diff_ar, crsmtch_ar)]

                    grbdata['%saperMag' % band] = Column(np.round(np.array(apmag_ar), 3) * u.ABmag,
                        description=aper_mag_desc
                    )
                    grbdata['%saperMag_Err' % band] = Column(np.round(np.array(apmag_err_ar), 3) * u.ABmag,
                        description=aper_mag_err_desc
                    )
                    grbdata['%saperME_CM' % band] = Column(np.round(np.array(aper_diff_ar), 3) * u.ABmag,
                        description=mag_diff_desc
                    )
                    grbdata['aperS_CM'] = Column(np.array(aper_flg_ar), description=flag_desc)

                if auto_mag_desc is not None:
                    auto_diff_ar = [mag_and_flag_grabber(md, sf, 'mag', 'AUTO') for md, sf in zip(diff_ar, crsmtch_ar)]
                    auto_flg_ar = [mag_and_flag_grabber(md, sf, 'flag', 'AUTO') for md, sf in zip(diff_ar, crsmtch_ar)]

                    grbdata['%sautoMag' % band] = Column(np.round(np.array(automag_ar), 3) * u.ABmag,
                        description=auto_mag_desc
                    )
                    grbdata['%sautoMag_Err' % band] = Column(np.round(np.array(automag_err_ar), 3) * u.ABmag,
                        description=auto_mag_err_desc
                    )
                    grbdata['%sautoME_CM' % band] = Column(np.round(np.array(auto_diff_ar), 3) * u.ABmag,
                        description=mag_diff_desc
                    )
                    grbdata['autoS_CM'] = Column(np.array(auto_flg_ar), description=flag_desc)

                grbdata['Radius'] = Column(np.round(np.array(rad_ar), decimals=2) * u.arcsec, description=rad_desc)
                grbdata['SNR'] = Column(np.round(np.array(snr_ar), decimals=2), description=snr_desc)
                grbdata['Elongation'] = Column(np.round(np.array(elon_ar), decimals=3), description=elon_desc)
                grbdata['Distance'] = Column(np.round(np.array(dist_ar), decimals=5) * u.arcsec, description=dist_desc)
                grbdata['Separation'] = Column(np.round(np.array(sep_ar), decimals=5) * u.arcsec, description=sep_desc)

                grbdata.write('%s_Multisource_%s_C%i_Data_%s_%s_loc_%d.ecsv' % (grbname, band, chip, survey, num, key),
                              overwrite=True)
                print(' Generated GRB data table & source DS9 regions!')
            else:
                print(
                    ' GRB source at inputted coords %s and %s not found in PRIME catalog, perhaps increase photoDistThresh or '
                    'alter sextractor params?' % (ra, dec))
    else:
        print(' idx size = %d' % len(idx_GRBcleanpsf))
        if len(idx_GRBcleanpsf) == 1:
            print(' GRB source at inputted coords %s and %s, rad = %s arcsec found!' % (ra, dec, photoDistThresh))

            if psf_mag_desc is not None:
                grb_mag = mag_ecsvcleanSources[idx_GRBcleanpsf]['%sMAG_PSF' % band][0]
                grb_magerr = mag_ecsvcleanSources[idx_GRBcleanpsf]['e_%sMAG_PSF' % band][0]
            if aper_mag_desc is not None:
                grb_mag_aper = mag_ecsvcleanSources[idx_GRBcleanpsf]['%sMAG_APER' % band][0]
                grb_magerr_aper = mag_ecsvcleanSources[idx_GRBcleanpsf]['e_%sMAG_APER' % band][0]
            if auto_mag_desc is not None:
                grb_mag_auto = mag_ecsvcleanSources[idx_GRBcleanpsf]['%sMAG_AUTO' % band][0]
                grb_magerr_auto = mag_ecsvcleanSources[idx_GRBcleanpsf]['e_%sMAG_AUTO' % band][0]

            grb_ra = mag_ecsvcleanSources[idx_GRBcleanpsf]['ALPHA_J2000'][0]
            grb_dec = mag_ecsvcleanSources[idx_GRBcleanpsf]['DELTA_J2000'][0]
            grb_rad = mag_ecsvcleanSources[idx_GRBcleanpsf]['FLUX_RADIUS'][0]
            grb_snr = mag_ecsvcleanSources[idx_GRBcleanpsf]['SNR_WIN'][0]
            grb_elon = mag_ecsvcleanSources[idx_GRBcleanpsf]['ELONGATION'][0]
            grb_dist = d2d[0].to(u.arcsec)
            grb_dist = grb_dist / u.arcsec

            print(
                ' Detected GRB ra = %.6f, dec = %.6f, with 50 percent flux radius (HWHM) = %.3f arcsec and SNR = %.3f' % (
                    grb_ra, grb_dec, grb_rad, grb_snr))

            all_magtypes = set(MAGTYPES.keys())
            for mag in sorted(all_magtypes):
                try:
                    print(f' %s {mag} magnitude of GRB is %.2f +/- %.2f' % (band,
                                                                     mag_ecsvcleanSources[idx_GRBcleanpsf][f'{band}MAG_{mag}'][0],
                                                                     mag_ecsvcleanSources[idx_GRBcleanpsf][f'e_{band}MAG_{mag}'][0]))
                except KeyError:
                    pass

            # survey crsmtch check
            idx_both, idx_bothcleanpsf, d2d_crs, d3d_crs = (
                massCatCoords.search_around_sky(mag_ecsvsourceCatCoords[idx_GRBcleanpsf], grb_rad * u.arcsec))
            prime_crs_cat = mag_ecsvcleanSources[idx_GRBcleanpsf]
            if len(idx_bothcleanpsf) > 0:
                # print(' Detected source crossmatched to existing %s source!' % survey)
                mag_diff_crs_ar, sep, survey_flg_ar = mag_diff_calc(survey_cat=good_cat_stars, prime_cat=prime_crs_cat,
                                                              survey_idx=idx_bothcleanpsf, prime_idx=idx_both,
                                                              d2d=d2d_crs, band=band)
            else:
                mag_diff_crs_ar = [99] * mag_col_num
                survey_flg_ar = non_cm_flag_calc(prime_source=prime_crs_cat, band=band,
                                                 survey_lim_mag=survey_lim_mag)
                sep = -1

            grbdata['RA'] = Column(np.round(np.array([grb_ra]), 5) * u.deg, description=ra_desc)
            grbdata['DEC'] = Column(np.round(np.array([grb_dec]), decimals=5) * u.deg, description=dec_desc)

            if psf_mag_desc is not None:
                grbdata['%spsfMag' % band] = Column(np.round(np.array([grb_mag]), 3) * u.ABmag,
                    description=psf_mag_desc
                )
                grbdata['%spsfMag_Err' % band] = Column(np.round(np.array([grb_magerr]), 3) * u.ABmag,
                    description=psf_mag_err_desc
                )
                grbdata['%spsfME_CM' % band] = Column(
                    np.round(np.array([mag_and_flag_grabber(mag_diff_crs_ar, survey_flg_ar, 'mag', 'PSF')]), 3) * u.ABmag,
                    description=mag_diff_desc
                )
                grbdata['psfS_CM'] = Column(mag_and_flag_grabber(mag_diff_crs_ar, survey_flg_ar, 'flag', 'PSF'),
                                            description=flag_desc)

            if aper_mag_desc is not None:
                grbdata['%saperMag' % band] = Column(np.round(np.array([grb_mag_aper]), 3) * u.ABmag,
                    description=aper_mag_desc
                )
                grbdata['%saperMag_Err' % band] = Column(np.round(np.array([grb_magerr_aper]), 3) * u.ABmag,
                    description=aper_mag_err_desc
                )
                grbdata['%saperME_CM' % band] = Column(
                    np.round(np.array([mag_and_flag_grabber(mag_diff_crs_ar, survey_flg_ar, 'mag', 'APER')]), 3) * u.ABmag,
                    description=mag_diff_desc
                )
                grbdata['aperS_CM'] = Column(mag_and_flag_grabber(mag_diff_crs_ar, survey_flg_ar, 'flag', 'APER'),
                                             description=flag_desc)

            if auto_mag_desc is not None:
                grbdata['%sautoMag' % band] = Column(np.round(np.array([grb_mag_auto]), 3) * u.ABmag,
                    description=auto_mag_desc
                )
                grbdata['%sautoMag_Err' % band] = Column(np.round(np.array([grb_magerr_auto]), 3) * u.ABmag,
                    description=auto_mag_err_desc
                )
                grbdata['%sautoME_CM' % band] = Column(
                    np.round(np.array([mag_and_flag_grabber(mag_diff_crs_ar, survey_flg_ar, 'mag', 'AUTO')]), 3) * u.ABmag,
                    description=mag_diff_desc
                )
                grbdata['autoS_CM'] = Column(mag_and_flag_grabber(mag_diff_crs_ar, survey_flg_ar, 'flag', 'AUTO'),
                                             description=flag_desc)

            grbdata['Radius'] = Column(np.round(np.array([grb_rad]), decimals=2) * u.arcsec, description=rad_desc)
            grbdata['SNR'] = Column(np.round(np.array([grb_snr]), decimals=2), description=snr_desc)
            grbdata['Elongation'] = Column(np.round(np.array([grb_elon]), decimals=3), description=elon_desc)
            grbdata['Distance'] = Column(np.round(np.array([grb_dist]), decimals=5) * u.arcsec, description=dist_desc)
            grbdata['Separation'] = Column(np.round(np.array([sep]), decimals=5) * u.arcsec, description=sep_desc)

            grbdata.write('%s_%s_C%i_Data_%s_%s.ecsv' % (grbname, band, chip, survey, num), overwrite=True)

            regprimename = source_reg_gen(grb_ra, grb_dec, rad=grb_rad)

            savename, threshname = grb_cutout(imageName, GRBcoords, photoDistThresh,
                                              regprimename=regprimename, regsurvname=regsurvname)

            html_gen(grbdata, directory, savename, threshname, band, survey, ra, dec, thresh, regsurvname, regprimename)

            print(' Generated GRB data table & source DS9 regions!')
        elif len(idx_GRBcleanpsf) > 1:
            print(' Multiple sources detected in search radius (ra = %s, dec = %s, rad = %s arcsec)'
                  ', refer to .ecsv file for source info!' % (ra, dec, photoDistThresh))
            mag_ar = []
            mag_err_ar = []
            apmag_ar = []
            apmag_err_ar = []
            automag_ar = []
            automag_err_ar = []
            ra_ar = []
            dec_ar = []
            rad_ar = []
            snr_ar = []
            elon_ar = []
            dist_ar = []
            sep_ar = []
            diff_ar = []
            crsmtch_ar = []
            idx_GRBcleanpsflist = idx_GRBcleanpsf.tolist()
            source_reg_gen()
            for i in idx_GRBcleanpsflist:
                if psf_mag_desc is not None:
                    grb_mag = mag_ecsvcleanSources[idx_GRBcleanpsf]['%sMAG_PSF' % band][idx_GRBcleanpsflist.index(i)]
                    mag_ar.append(grb_mag)
                    grb_magerr = mag_ecsvcleanSources[idx_GRBcleanpsf]['e_%sMAG_PSF' % band][idx_GRBcleanpsflist.index(i)]
                    mag_err_ar.append(grb_magerr)
                if aper_mag_desc is not None:
                    grb_mag_aper = mag_ecsvcleanSources[idx_GRBcleanpsf]['%sMAG_APER' % band][idx_GRBcleanpsflist.index(i)]
                    apmag_ar.append(grb_mag_aper)
                    grb_magerr_aper = mag_ecsvcleanSources[idx_GRBcleanpsf]['e_%sMAG_APER' % band][idx_GRBcleanpsflist.index(i)]
                    apmag_err_ar.append(grb_magerr_aper)
                    # grb_mag_aper2 = mag_ecsvcleanSources[idx_GRBcleanpsf]['%sMAG_APER' % band][:, 0][
                    #     idx_GRBcleanpsflist.index(i)]
                if auto_mag_desc is not None:
                    grb_mag_auto = mag_ecsvcleanSources[idx_GRBcleanpsf]['%sMAG_AUTO' % band][idx_GRBcleanpsflist.index(i)]
                    automag_ar.append(grb_mag_auto)
                    grb_magerr_auto = mag_ecsvcleanSources[idx_GRBcleanpsf]['e_%sMAG_AUTO' % band][
                        idx_GRBcleanpsflist.index(i)]
                    automag_err_ar.append(grb_magerr_auto)
                grb_ra = mag_ecsvcleanSources[idx_GRBcleanpsf]['ALPHA_J2000'][idx_GRBcleanpsflist.index(i)]
                ra_ar.append(grb_ra)
                grb_dec = mag_ecsvcleanSources[idx_GRBcleanpsf]['DELTA_J2000'][idx_GRBcleanpsflist.index(i)]
                dec_ar.append(grb_dec)
                grb_rad = mag_ecsvcleanSources[idx_GRBcleanpsf]['FLUX_RADIUS'][idx_GRBcleanpsflist.index(i)]
                rad_ar.append(grb_rad)
                grb_snr = mag_ecsvcleanSources[idx_GRBcleanpsf]['SNR_WIN'][idx_GRBcleanpsflist.index(i)]
                snr_ar.append(grb_snr)
                grb_elon = mag_ecsvcleanSources[idx_GRBcleanpsf]['ELONGATION'][idx_GRBcleanpsflist.index(i)]
                elon_ar.append(grb_elon)
                dist = (d2d[idx_GRBcleanpsflist.index(i)]).to(u.arcsec)
                dist = dist / u.arcsec
                dist_ar.append(dist)

                # survey crsmtch check
                checkcoords = SkyCoord(ra=[mag_ecsvsourceCatCoords[i].ra.deg], dec=[mag_ecsvsourceCatCoords[i].dec.deg],
                                       frame='icrs', unit='degree')
                idx_both, idx_bothcleanpsf, d2d_crs, d3d_crs = massCatCoords.search_around_sky(
                    checkcoords, grb_rad * u.arcsec)
                prime_crs_cat = mag_ecsvcleanSources[i]
                if len(idx_bothcleanpsf) > 0:
                    mag_diff_crs_ar, sep, survey_flg_ar = mag_diff_calc(survey_cat=good_cat_stars, prime_cat=prime_crs_cat,
                                                                  survey_idx=idx_bothcleanpsf, prime_idx=idx_both,
                                                                  d2d=d2d_crs, band=band)
                else:
                    mag_diff_crs_ar = [99] * mag_col_num
                    survey_flg_ar = non_cm_flag_calc(prime_source=prime_crs_cat, band=band,
                                                     survey_lim_mag=survey_lim_mag)
                    sep = -1

                diff_ar.append(mag_diff_crs_ar)
                crsmtch_ar.append(survey_flg_ar)
                sep_ar.append(sep)

                regprimename = source_reg_gen(grb_ra, grb_dec, rad=grb_rad, append=True)

            grbdata['RA'] = Column(np.round(np.array(ra_ar), 5) * u.deg, description=ra_desc)
            grbdata['DEC'] = Column(np.round(np.array(dec_ar), decimals=5) * u.deg, description=dec_desc)

            if psf_mag_desc is not None:
                psf_diff_ar = [mag_and_flag_grabber(md, sf, 'mag', 'PSF') for md, sf in zip(diff_ar, crsmtch_ar)]
                psf_flg_ar = [mag_and_flag_grabber(md, sf, 'flag', 'PSF') for md, sf in zip(diff_ar, crsmtch_ar)]

                grbdata['%spsfMag' % band] = Column(np.round(np.array(mag_ar), 3) * u.ABmag,
                    description=psf_mag_desc
                )
                grbdata['%spsfMag_Err' % band] = Column(np.round(np.array(mag_err_ar), 3) * u.ABmag,
                    description=psf_mag_err_desc
                )
                grbdata['%spsfME_CM' % band] = Column(np.round(np.array(psf_diff_ar), 3) * u.ABmag,
                    description=mag_diff_desc
                )
                grbdata['psfS_CM'] = Column(np.array(psf_flg_ar), description=flag_desc)

            if aper_mag_desc is not None:
                aper_diff_ar = [mag_and_flag_grabber(md, sf, 'mag', 'APER') for md, sf in zip(diff_ar, crsmtch_ar)]
                aper_flg_ar = [mag_and_flag_grabber(md, sf, 'flag', 'APER') for md, sf in zip(diff_ar, crsmtch_ar)]

                grbdata['%saperMag' % band] = Column(np.round(np.array(apmag_ar), 3) * u.ABmag,
                    description=aper_mag_desc
                )
                grbdata['%saperMag_Err' % band] = Column(np.round(np.array(apmag_err_ar), 3) * u.ABmag,
                    description=aper_mag_err_desc
                )
                grbdata['%saperME_CM' % band] = Column(np.round(np.array(aper_diff_ar), 3) * u.ABmag,
                    description=mag_diff_desc
                )
                grbdata['aperS_CM'] = Column(np.array(aper_flg_ar), description=flag_desc)

            if auto_mag_desc is not None:
                auto_diff_ar = [mag_and_flag_grabber(md, sf, 'mag', 'AUTO') for md, sf in zip(diff_ar, crsmtch_ar)]
                auto_flg_ar = [mag_and_flag_grabber(md, sf, 'flag', 'AUTO') for md, sf in zip(diff_ar, crsmtch_ar)]

                grbdata['%sautoMag' % band] = Column(np.round(np.array(automag_ar), 3) * u.ABmag,
                    description=auto_mag_desc
                )
                grbdata['%sautoMag_Err' % band] = Column(np.round(np.array(automag_err_ar), 3) * u.ABmag,
                    description=auto_mag_err_desc
                )
                grbdata['%sautoME_CM' % band] = Column(np.round(np.array(auto_diff_ar), 3) * u.ABmag,
                    description=mag_diff_desc
                )
                grbdata['autoS_CM'] = Column(np.array(auto_flg_ar), description=flag_desc)

            grbdata['Radius'] = Column(np.round(np.array(rad_ar), decimals=2) * u.arcsec, description=rad_desc)
            grbdata['SNR'] = Column(np.round(np.array(snr_ar), decimals=2), description=snr_desc)
            grbdata['Elongation'] = Column(np.round(np.array(elon_ar), decimals=3), description=elon_desc)
            grbdata['Distance'] = Column(np.round(np.array(dist_ar), decimals=5) * u.arcsec, description=dist_desc)
            grbdata['Separation'] = Column(np.round(np.array(sep_ar), decimals=5) * u.arcsec, description=sep_desc)

            grbdata.write('%s_Multisource_%s_C%i_Data_%s_%s.ecsv' % (grbname, band, chip, survey, num), overwrite=True)

            savename, threshname = grb_cutout(imageName, GRBcoords, photoDistThresh,
                                              regprimename=regprimename, regsurvname=regsurvname)

            html_gen(grbdata, directory, savename, threshname, band, survey, ra, dec, thresh, regsurvname, regprimename)

            print(' Generated GRB data table & source DS9 regions!')
        else:
            print(
                ' GRB source at inputted coords %s and %s not found in PRIME catalog, perhaps increase photoDistThresh or '
                'alter sextractor params?' % (ra, dec))


# New Source Search
def newsourcesearch(ra, dec, imageName, survey, band, thresh, massCatCoords, good_cat_stars, directory, chip, grbname=defaults['grb_name'],
                    mag_low_lim=PHOTOMETRY_MAG_LOWER_LIMIT, magtype=defaults['magtype'], **kwargs):
    """
    GRB functionality for large error radii.  When used in photometry.py, if the input grb threshold in > 60", GRB
    functionality defaults to this.  This does a series of pruning methods: 1. removes all PRIME sources crossmatched to
    input survey catalog.  2. prunes all sources outside the input error threshold. 3. prunes out all sources w/ mags
    brighter than the lower mag limit & mags dimmer than the survey limiting mag.  Outputs an .ecsv of sources left
    as candidates for the target event.


    Parameters
    ----------
    ra: float
        RA coordinate for search
    dec: float
        Dec coordinate for search
    imageName: str
        Filename of input image
    band: str
        Filter observation was taken in
    thresh: float
        Diameter of threshold to search for sources
    massCatCoords: SkyCoord
        Coordinate data of queried survey sources
    good_cat_stars: Table
        Astropy Table of AB-corrected survey sources (incl. ra, dec, mag, & mag error columns)
    directory: str
        Directory the input image is stored in
    chip: int
        Detector number of input image
    grbname: str
        Optional, custom name for all GRB data products (default = 'GRB')
    mag_low_lim: float
        Optional, specify a bright mag cutoff limit (default = 12.5)
    magtype: str
        Optional, specify the magtype for pruning (default = AUTO)
    """

    source_ra = ra
    source_dec = dec
    ab_cat_stars = good_cat_stars
    header = fits.getheader(imageName)
    w = WCS(header)

    if grbname != defaults['grb_name']:
        grb_name = grbname
        newsrcname = grbname
    else:
        grb_name = defaults['grb_name']
        newsrcname = 'New_PRIME'

    if not imageName.startswith(f'coadd.Open-{band}') and not imageName.endswith('.fits'):
        num = 'img'
    else:
        num = imageName[-16:-8]

    # large region postage stamp cutout fctn (png & fits)
    def grb_cutout(imageName, GRBcoords, photoDistThresh, regprimename=None, regsurvname=None):
        imgdata = fits.getdata(imageName)
        img = fits.open(imageName)
        head = img[0].header
        w = WCS(head)

        savename = '%s_%s_Cutout_%s_%s' % (grb_name, band, survey, num)
        threshname = '%s_query_thresh.reg' % grb_name

        size = 4 * photoDistThresh * u.arcsec
        try:
            cutout = Cutout2D(imgdata, GRBcoords, size, wcs=w, copy=True)
            region = CircleSkyRegion(center=GRBcoords[0], radius=Angle(thresh, unit='arcsec'))
            pix_region = region.to_pixel(cutout.wcs)
        except astropy.nddata.utils.NoOverlapError:
            print(' Area of GRB threshold not found within image, cannot generate cutout!')
            return savename, threshname

        mean, median, sigma_cut = sigma_clipped_stats(cutout.data)
        plt.figure(10, figsize=(8, 8))
        plt.imshow(cutout.data, vmin=median - 3 * sigma_cut, vmax=median + 3 * sigma_cut, origin='lower', cmap='viridis')
        pix_region.plot(color='cyan', ls='--', label='Input GRB threshold')

        if regprimename:
            primeregs = open('%s_srcs.reg' % newsrcname, 'r')
            plt_primeregs = []
            primeallregs = [reg for reg in primeregs if reg != 'fk5\n']
            for reg in primeallregs:
                nums = re.findall(r'[-+]?\d*\.?\d+', reg)
                srcra = float(nums[0])
                srcdec = float(nums[1])
                srcrad = float(nums[2])
                srccoords = SkyCoord(ra=[srcra], dec=[srcdec], frame='icrs', unit='degree')
                srcreg = CircleSkyRegion(center=srccoords[0], radius=Angle(srcrad, unit='arcsec'))
                plt_primeregs.append(srcreg)

            for reg in plt_primeregs:
                pix_reg = reg.to_pixel(cutout.wcs)
                pix_reg.plot(color='red', ls='-', label='PRIME Source')

        if regsurvname:
            survregs = open('%s_%s_srcs.reg' % (grb_name, survey), 'r')
            plt_survregs = []
            survallregs = [reg for reg in survregs if reg != 'fk5\n']
            for reg in survallregs:
                nums = re.findall(r'[-+]?\d*\.?\d+', reg)
                srcra = float(nums[0])
                srcdec = float(nums[1])
                srcrad = float(nums[2])
                srccoords = SkyCoord(ra=[srcra], dec=[srcdec], frame='icrs', unit='degree')
                srcreg = CircleSkyRegion(center=srccoords[0], radius=Angle(srcrad, unit='arcsec'))
                plt_survregs.append(srcreg)

            for reg in plt_survregs:
                pix_reg = reg.to_pixel(cutout.wcs)
                pix_reg.plot(color='magenta', ls='-', label='%s Source' % survey)

        handles, labels = plt.gca().get_legend_handles_labels()
        seen = set()
        filtered_handles = []
        filtered_labels = []
        for h, l in zip(handles, labels):
            if l not in seen:
                filtered_handles.append(h)
                filtered_labels.append(l)
                seen.add(l)
        plt.legend(filtered_handles, filtered_labels, loc='best')
        plt.savefig(savename + '.png', dpi=300)
        plt.close()

        fits.writeto(savename + '.fits', cutout.data, cutout.wcs.to_header(), overwrite=True)
        return savename, threshname

    # new source reg gen
    def source_reg_gen(src_ra=0, src_dec=0, rad=2, src_survey=None, append=False):
        if not src_survey:
            name = '%s_srcs.reg' % newsrcname
            color = 'red'
        else:
            name = '%s_%s_srcs.reg' % (grb_name, survey)
            color = 'yellow'
        if not append:
            newtext = open(name, 'w+')
            newtext.write('fk5')
            if src_ra != 0 and src_dec != 0:
                newtext.write(f'\ncircle({src_ra}, {src_dec}, {rad}") # color={color}')
        else:
            newtext = open(name, 'a')
            newtext.write(f'\ncircle({src_ra}, {src_dec}, {rad}") # color={color}')

        return name

    # html gen for new sources
    def html_gen(data, directory, savename, threshname, band, survey, ra, dec, thresh, primename=None):
        df = data.to_pandas()
        tbl_html = df.to_html(index=False, classes="my-table")

        with open(savename + '.png', "rb") as img_file:
            encoded = base64.b64encode(img_file.read()).decode("utf-8")

        fits_items = [savename + '.fits', threshname, primename]
        fits_items = [os.path.join(directory, f) for f in fits_items if f is not None]

        base = fits_items[0]
        regions = " ".join(f"-regions {f}" for f in fits_items[1:])

        ds9_command = f"ds9 {base} {regions} &"

        final_html = f"""
        <div style="text-align: center; font-family: Arial, sans-serif;">

            <h2>GRB Information</h2>
            <p style="margin-top: 0; margin-bottom: 10px; font-size: 16px; color: #555;">
                RA = {ra}, Dec = {dec}, threshold = {thresh}"
            </p>

            <img src="data:image/png;base64,{encoded}" alt="GRB Cutout Region" width="600" style="margin-bottom: 10px;">

            <p style="margin-top: 0px; margin-bottom: 10px; font-size: 16px; color: #555;">
                To see source regions, open GRB stamp in DS9 through terminal: 
            </p>
            <p style="margin-top: 0px; margin-bottom: 20px; font-size: 14px; color: #030303;">
                {ds9_command}
            </p>

            <div style="display: inline-block; text-align: left;">
                {tbl_html}
            </div>

        </div>
        """

        with open('%s_Source_%s_Data_%s_%s.html' % (newsrcname, band, survey, num), "w") as f:
            f.write(final_html)

    # crossmatch for all detected sources
    print('Large grb radius inputted, using new source search to find'
          ' detected sources brighter than survey lim mag w/ no crossmatch')

    # sexigesimal conversion
    try:
        float(source_ra)
        source_ra = source_ra
        source_dec = source_dec
        deci_sky_coords = SkyCoord(ra=[source_ra], dec=[source_dec], frame='icrs', unit='degree')
    except ValueError:
        coords = source_ra + ' ' + source_dec
        print('Sexagesimal RA = %s & Dec = %s' % (source_ra, source_dec))

        deci_sky_coords = SkyCoord(coords, frame='icrs', unit=(u.hourangle, u.deg))
        deci_coords = deci_sky_coords.to_string()
        deci_coords = deci_coords.split(' ')

        source_ra = deci_coords[0]
        source_dec = deci_coords[1]

    if len(ab_cat_stars.colnames) > 6:
        RA = 'ALPHA_J2000'
        DEC = 'DELTA_J2000'
    else:
        RA = ab_cat_stars.colnames[0]
        DEC = ab_cat_stars.colnames[1]

    massCatCoords = SkyCoord(ra=ab_cat_stars[RA], dec=ab_cat_stars[DEC], frame='icrs', unit='degree')
    print('Catalog cropped #:', len(massCatCoords))

    # initial crossmatch
    mag_ecsvname = '%s.%s.ecsv' % (imageName, survey)
    mag_ecsvtable = ascii.read(mag_ecsvname)
    mag_ecsvSources = mag_ecsvtable[(mag_ecsvtable['FLAGS'] == 0) & (mag_ecsvtable['FLAGS_MODEL'] == 0)]
    # print(len(mag_ecsvSources))
    if RA == 'ALPHA_J2000':
        mag_ecsvsourceCatCoords = SkyCoord(ra=mag_ecsvSources['ALPHA_J2000'], dec=mag_ecsvSources['DELTA_J2000'],
                                           frame='icrs',
                                           unit='degree')
    else:
        mag_ecsvsourceCatCoords = utils.pixel_to_skycoord(mag_ecsvSources['X_IMAGE'], mag_ecsvSources['Y_IMAGE'], w, origin=1)

    photoDistThresh = 1.0   # set higher to combat offset wcs in certain sources, maybe change for denser fields?
    idx_psfimage_noclean, idx_psfmass_noclean, d2d, d3d = massCatCoords.search_around_sky(mag_ecsvsourceCatCoords,
                                                                          photoDistThresh * u.arcsec)

    mask = np.ones(len(mag_ecsvSources), dtype=bool)
    mask[idx_psfimage_noclean] = False
    PSFsources_nomatch = mag_ecsvSources[mask]  # removing previous crossmatched sources
    print('# of sources found after removing survey crossmatches: %i' % len(PSFsources_nomatch))

    sourcecoords = SkyCoord(ra=[source_ra], dec=[source_dec], frame='icrs', unit='degree')
    PSFsources_nomatchCatCoords = SkyCoord(ra=PSFsources_nomatch['ALPHA_J2000'], dec=PSFsources_nomatch['DELTA_J2000'],
                                       frame='icrs',
                                       unit='degree')
    idx_inputcoords, idx_PSFsources_nomatch, d2dd, d3dd = PSFsources_nomatchCatCoords.search_around_sky(sourcecoords,
                                                                                                      thresh * u.arcsec)
    PSFsources_nomatch = PSFsources_nomatch[idx_PSFsources_nomatch]     # implementing error radius
    print('# of sources found after implementing err radius: %i' % len(PSFsources_nomatch))

    # lim mag pruning
    for f in PHOTOMETRY_LIM_MAGS.keys():
        if survey == f:
            lim_mag = PHOTOMETRY_LIM_MAGS[f]
            break
    else:
        # limiting mag est.
        # lim_mag = ab_cat_stars['%sMAG_PSF' % band][(ab_cat_stars['SNR_WIN'] < 6) & (ab_cat_stars['SNR_WIN'] > 5)]
        lim_mag = 19.7
        print(f'Error in finding lim mag for survey catalogs, no matching catalog found? '
              f'Using generous PRIME lim: {lim_mag}')

    print(f' {survey} limiting mag: {lim_mag}')
    PSFsources_new = PSFsources_nomatch[(PSFsources_nomatch[f'{band}MAG_{MAGTYPES[magtype]}'] < lim_mag) &
                                        (PSFsources_nomatch[f'{band}MAG_{MAGTYPES[magtype]}'] > mag_low_lim)]
    print('# of sources found after removing sources dimmer than %.2f & brighter than %.2f in %s mag: %i'
          % (mag_low_lim, lim_mag, MAGTYPES[magtype], len(PSFsources_new)))

    # Table gen

    if len(PSFsources_new) > 0:
        PSFsources_new.write('%s_Sources.%s.%s.%s.ecsv' % (newsrcname, imageName, survey, num), overwrite=True)
        print('New source full catalog written!')

        all_magtypes = set(MAGTYPES.keys())
        prime_mag_cols = sorted([col for col in PSFsources_new.colnames if any(mag in col for mag in all_magtypes)
                                 and band in col and f'{band}MAG' in col and 'e_' not in col])

        mag_auto_ar = []
        mag_auto_err_ar = []
        mag_psf_ar = []
        mag_psf_err_ar = []
        mag_aper_ar = []
        mag_aper_err_ar = []
        ra_ar = []
        dec_ar = []
        rad_ar = []
        snr_ar = []
        source_reg_gen()

        for i in PSFsources_new:
            if f'{band}MAG_AUTO' in prime_mag_cols:
                grb_mag_auto = i[f'{band}MAG_AUTO']
                mag_auto_ar.append(grb_mag_auto)
                grb_mag_auto_err = i[f'e_{band}MAG_AUTO']
                mag_auto_err_ar.append(grb_mag_auto_err)
            if f'{band}MAG_PSF' in prime_mag_cols:
                grb_mag_psf = i[f'{band}MAG_PSF']
                mag_psf_ar.append(grb_mag_psf)
                grb_mag_psf_err = i[f'e_{band}MAG_PSF']
                mag_psf_err_ar.append(grb_mag_psf_err)
            if f'{band}MAG_APER' in prime_mag_cols:
                grb_mag_aper = i[f'{band}MAG_APER']
                mag_aper_ar.append(grb_mag_aper)
                grb_mag_aper_err = i[f'e_{band}MAG_APER']
                mag_aper_err_ar.append(grb_mag_aper_err)
            grb_ra = i['ALPHA_J2000']
            ra_ar.append(grb_ra)
            grb_dec = i['DELTA_J2000']
            dec_ar.append(grb_dec)
            grb_rad = i['FLUX_RADIUS']
            rad_ar.append(grb_rad)
            grb_snr = i['SNR_WIN']
            snr_ar.append(grb_snr)

            regprimename = source_reg_gen(src_ra=grb_ra, src_dec=grb_dec, rad=grb_rad, append=True)

        grbdata = Table()
        grbdata['RA'] = np.round(np.array(ra_ar), decimals=5) * u.deg
        grbdata['DEC'] = np.round(np.array(dec_ar), decimals=5) * u.deg
        if f'{band}MAG_AUTO' in prime_mag_cols:
            grbdata[f'{band}autoMag'] = np.round(np.array(mag_auto_ar), decimals=3) * u.ABmag
            grbdata[f'{band}autoMag_Err'] = np.round(np.array(mag_auto_err_ar), decimals=3) * u.ABmag
        if f'{band}MAG_PSF' in prime_mag_cols:
            grbdata[f'{band}psfMag'] = np.round(np.array(mag_psf_ar), decimals=3) * u.ABmag
            grbdata[f'{band}psfMag_Err'] = np.round(np.array(mag_psf_err_ar), decimals=3) * u.ABmag
        if f'{band}MAG_APER' in prime_mag_cols:
            grbdata[f'{band}aperMag'] = np.round(np.array(mag_aper_ar), decimals=3) * u.ABmag
            grbdata[f'{band}aperMag_Err'] = np.round(np.array(mag_aper_err_ar), decimals=3) * u.ABmag
        grbdata['Radius'] = np.round(np.array(rad_ar), decimals=2) * u.arcsec
        grbdata['SNR'] = np.round(np.array(snr_ar), decimals=2)

        print('New source catalog writen!')
        grbdata.write('%s_Source_%s_Data_%s_%s.ecsv' % (newsrcname, band, survey, num), overwrite=True)

        # generate stamp & html
        print('Generating location cutout and html files!')
        savename, threshname = grb_cutout(imageName=imageName, GRBcoords=deci_sky_coords, photoDistThresh=thresh,
                                          regprimename=regprimename)
        html_gen(grbdata, directory, savename, threshname, band, survey, ra=source_ra,
                 dec=source_dec, thresh=thresh, primename=regprimename)

    else:
        print('No sources remaining after pruning!  Cannot write new source catalog!')