"""
Check-plot functions for photometry
"""

import numpy as np
import numpy.ma as ma
from astropy.table import Column
from astropy.table import Table
from astropy.io import fits

import matplotlib.pyplot as plt
import statsmodels.api as sm
from scipy.stats import skew
import re

from photometrus.settings import MAGTYPES

from photometrus.photometry.photom_utils import end_name_gen, open_fits_robust
from photometrus.photometry.photom_fits import single_fit_calc, photometric_fit_calc

from photometrus.utils.defaults import PROCESSING_DEFAULTS as defaults


# single photometry plots for int calibration
def single_plots(cleanPSFsources, PSFsources, data, imageName, survey, band, good_cat_stars, idx_mass,
                     idx_image, sigma, weights_noclip, clipped, crop=defaults['crop'], magtype=defaults['magtype']):
    """
    For use in multiprocessing photometry calibration, generates plots given any single magtype.  Plots currently
    include: flux density vs. ab mag, histogram of sigma clipped crossmatched residuals, WLS photometric fit line w/
    crossmatched source histogram (both sigma clipped and not), and half max histogram limiting mag plot.


    Parameters
    ----------
    cleanPSFsources: Table
        PRIME clean source table (pruning done to avoid bad sources for crossmatching and zero pt calculation)
    PSFsources: Table
        PRIME source table (unpruned)
    data: np.ndarray
        Input image data
    imageName: str
        Filename of input image
    survey: str
        Utilized survey for crossmatching and calibration
    band: str
        Filter observation was taken in
    good_cat_stars: Table
        Astropy Table of AB-corrected survey sources (incl. ra, dec, mag, & mag error columns)
    idx_mass: list
        Indices into queried survey table for sources crossmatched to PRIME
    idx_image: list
        Indices into PRIME source table for sources crossmatched to queried survey
    sigma: float
        Sigma clip value
    weights_noclip: Column
        Weights for PRIME sources before sigma clipping
    clipped: MaskedArray
        PRIME sources remaining after sigma clip
    crop: int
        Optional, Amount of corner cropping applied to clean source catalog to avoid bad corner noise (default = 300 px)
    magtype: str
        Optional, specify the magtype for pruning (default = AUTO)
    """

    # appropriate mag column
    colnames = good_cat_stars.colnames
    if len(colnames) > 6:
        magcol = f'{band}MAG_{magtype}'
        magerrcol = f'{band}MAG_{magtype}'
    else:
        magcol = colnames[2]
        magerrcol = colnames[3]

    try:
        chip = fits.getval(imageName,'CHIP')
    except KeyError:
        chip = 0
    if not imageName.startswith(f'coadd.Open-{band}') and not imageName.endswith('.fits'):
        num = 'img'
    else:
        num = imageName[-16:-8]

    # PSF-specific naming for parallel running
    end_name = end_name_gen('png', magtype=MAGTYPES[magtype])

    def predict_y_for(x, m, b):
        return m * x + b

    # PHOTOMETRIC FIT LINE MODELS

    model, model2 = single_fit_calc(cleanPSFsources, band, good_cat_stars, idx_mass, idx_image,
                         weights_noclip, clipped, sigma, with_plots=True, magtype=magtype)

    m = model.params[1]
    m_err = model.bse[1]
    b = model.params[0]
    b_err = model.bse[0]
    rsquare = model.rsquared
    rss = model.ssr
    model_resid = model.resid

    m2 = model2.params[1]
    m2err = model2.bse[1]
    b2 = model2.params[0]
    b2err = model2.bse[0]
    avg2 = np.average(model2.resid, weights=weights_noclip)

    # BINNED RESIDUAL STATISTICS

    # residual fit 3 sig clip - bin errors and stats
    mags = range(10, 22)
    x_arr = np.arange(10 + 0.5, 22 + 0.5, 1)

    res_errs = []
    res_means = []
    res_ranges = []
    res_skews = []
    res_nums = []
    for i in mags:
        mask = ((cleanPSFsources[f'{band}MAG_{MAGTYPES[magtype]}'][idx_image][~clipped.mask] > i) &
                     (cleanPSFsources[f'{band}MAG_{MAGTYPES[magtype]}'][idx_image][~clipped.mask] < i + 1))
        mask = model_resid[mask]
        res_err = np.std(mask)  # spread
        res_range = '%s - %s' % (i, i + 1)
        res_mean = np.nanmean(mask)  # mean
        res_skew = skew(mask, bias=False)  # skew
        res_num = len(mask)
        if ma.is_masked(res_err):
            res_err = 0
        res_errs.append(res_err)
        res_means.append(res_mean)
        res_ranges.append(res_range)
        res_skews.append(res_skew)
        res_nums.append(res_num)
    res_errs = np.nan_to_num(np.array(res_errs))
    res_means = np.nan_to_num(np.array(res_means))
    res_skews = np.nan_to_num(np.array(res_skews))
    res_errs_min = np.min((res_errs[res_errs != 0]))
    res_errs_max = np.max(res_errs)

    bintable = Table()
    bintable['Bin Range (mag)'] = res_ranges
    bintable['Mean (mag)'] = res_means
    bintable['Spread (mag)'] = res_errs
    bintable['Skew'] = res_skews
    bintable['Source Number'] = res_nums
    bintable.write('Resid_%s-sig_Data_%s_C%s_%s%s' % (sigma, band, chip, survey, end_name_gen(magtype=MAGTYPES[magtype])), overwrite=True)

    # ALL PLOTS

    plt.close('all')

    # PRIME flux vs catalog AB mag for crossmatches
    x_flx_lin = cleanPSFsources[f'{MAGTYPES[magtype]}_FLUX_DENSITY'][idx_image][~clipped.mask]
    x_flx = np.log(x_flx_lin)
    y_flx = good_cat_stars['%s' % magcol][idx_mass][~clipped.mask]
    x_const_flx = sm.add_constant(x_flx)
    model_flx = sm.WLS(y_flx, x_const_flx, weights=weights_noclip[~clipped.mask]).fit()
    # print(model_flx.params)
    m_flx = model_flx.params[1]
    m_flxerr = model_flx.bse[1]
    b_flx = model_flx.params[0]
    b_flxerr = model_flx.bse[0]

    x_grid = np.linspace(x_flx_lin.min(), x_flx_lin.max(), 100)
    log_x_grid = np.log(x_grid)
    log_x_grid_const = sm.add_constant(log_x_grid)
    y_pred = model_flx.predict(log_x_grid_const)

    plt.figure(9, figsize=(8, 8))
    plt.plot(cleanPSFsources[f'{MAGTYPES[magtype]}_FLUX_DENSITY'][idx_image][~clipped.mask],
             good_cat_stars['%s' % magcol][idx_mass][~clipped.mask],
             'r.', markersize=14, markeredgecolor='black')
    plt.plot(x_grid, y_pred, c='b')

    flx_txt = ('slope = %.4f' % m_flx + '\nslope err = %.4f' % m_flxerr +
               '\nint = %.4f' % b_flx + '\nint err = %.4f' % b_flxerr)

    plt.xlim(10, 50000)
    plt.ylim(12, 20.5)
    plt.title('PRIME Flux Density vs %s AB mag - %s Sigma Clip' % (survey, sigma))
    plt.xlabel(r'PRIME Flux Density ($\mu$Jy)', fontsize=15)
    plt.ylabel('%s %s AB Mags' % (survey, band), fontsize=15)
    plt.grid()
    plt.xscale('log')
    flx_box = dict(facecolor='white')
    plt.text(1000, 19, flx_txt, fontsize=12, bbox=flx_box)
    plt.savefig('%s_C%s_flux_mag_plot_sig_%s%s' % (survey, chip, num, end_name))
    plt.clf()
    print('Saved flux v. mag plot to dir!')

    # res plot y int, histogram
    if len(idx_image) >= 75000:
        bin_num_int = round(len(idx_image) / 500)
    elif 50000 <= len(idx_image) <= 75000:
        bin_num_int = round(len(idx_image) / 400)
    elif 25000 <= len(idx_image) <= 50000:
        bin_num_int = round(len(idx_image) / 300)
    elif 5000 <= len(idx_image) <= 25000:
        bin_num_int = round(len(idx_image) / 75)
    elif 1000 <= len(idx_image) <= 5000:
        bin_num_int = round(len(idx_image) / 50)
    elif len(idx_image) <= 1000:
        bin_num_int = 50

    # magtype specific coloring, etc.
    alpha = 0.7 if model else None
    cmap = 'gist_heat_r'
    if magtype == MAGTYPES['PSF']:
        cmap = 'gist_heat_r'
        lm_colors = ['red', 'blue']
    elif magtype == MAGTYPES['AUTO']:
        cmap = 'gist_earth_r'
        lm_colors = ['blue', 'black']
    elif magtype == MAGTYPES['APER']:
        cmap = 'pink_r'
        lm_colors = ['green', 'black']


    fig, ax2 = plt.subplots(1, 1, figsize=(9, 8))
    fig.suptitle(f'{survey} Residuals - {sigma} Sigma Clip - Density Histogram')

    hist = ax2.hist2d(
        x=good_cat_stars['%s' % magcol][idx_mass][~clipped.mask], y=model_resid,
        bins=[bin_num_int, bin_num_int], range=[[10, 21], [-1, 1]],
        cmap=cmap
    )

    cbar2 = fig.colorbar(hist[3], ax=ax2, pad=0.03)
    cbar2.set_label(f'{MAGTYPES[magtype]} Density')

    sig = ax2.scatter(x_arr, res_errs, marker='_', s=1625, c='black',
                          label=fr'{MAGTYPES[magtype]} ap. photom 1 $\sigma$ range = [%.3f - %.3f]' % (
                              res_errs_min, res_errs_max))
    ax2.scatter(x_arr, -res_errs, marker='_', s=1625, c='black')

    ax2.axhline(y=0, color='black', linestyle='--', linewidth=1)
    ax2.set_xlim(10, 21)
    ax2.set_ylim(-1, 1)
    ax2.yaxis.set_tick_params(labelleft=True)
    ax2.set_title(f"{MAGTYPES[magtype]} Fit Residuals")
    ax2.set_xlabel(f"{survey} {band} Mags")
    ax2.set_ylabel("Residuals")

    info = (
            f'{MAGTYPES[magtype]} fit'
            + '\nslope = %.4f +/- %.4f' % (m, m_err)
            + '\nintercept = %.3f +/- %.3f' % (b, b_err)
            + '\nR$^{2}$ = %.3f' % rsquare
            + '\nRSS = %d' % rss
            + '\nn_sources = %i' % len(cleanPSFsources[f'{MAGTYPES[magtype]}_FLUX_DENSITY'][idx_image][~clipped.mask])
    )

    ax2.text(10.5, 0.5, info, fontsize=9,
             bbox=dict(facecolor='white', edgecolor='black', pad=5.0))

    ax2.legend([sig], [sig.get_label()], loc='lower left', markerscale=0.5)

    plt.savefig('%s_C%s_residual_plot_int_hist_%s%s' % (survey, chip, num, end_name), dpi=300)
    plt.close(fig)

    print('Saved y-int residual plots to dir!')

    # WLS fit line over data plot

    txt = ('slope = %.4f' % m2 + '\nslope err = %.4f' % m2err + '\nint = %.4f' % b2 + '\nint err = %.4f' % b2err +
           '\nn_sources = %i' % len(cleanPSFsources[idx_image]))

    # WLS hist density plot
    if len(idx_image) >= 75000:
        bin_num = round(len(idx_image) / 500)
    elif 50000 <= len(idx_image) <= 75000:
        bin_num = round(len(idx_image) / 350)
    elif 25000 <= len(idx_image) <= 50000:
        bin_num = round(len(idx_image) / 150)
    elif 5000 <= len(idx_image) <= 25000:
        bin_num = round(len(idx_image) / 50)
    elif 1000 <= len(idx_image) <= 5000:
        bin_num = round(len(idx_image) / 20)
    elif len(idx_image) <= 1000:
        bin_num = 100

    plt.figure(5, figsize=(10, 8))
    plt.clf()
    plt.xlim(10, 22)
    plt.ylim(10, 22)
    plt.title(f'{survey} vs PRIME {MAGTYPES[magtype]} Mag w/ Weighted Fit - Density Histogram')
    plt.grid()
    plt.ylabel('PRIME %s Mags' % band, fontsize=15)
    plt.xlabel('%s %s Mags' % (survey, band), fontsize=15)
    plt.hist2d(x=good_cat_stars['%s' % magcol][idx_mass], y=cleanPSFsources[f'{band}MAG_{MAGTYPES[magtype]}'][idx_image],
               bins=[bin_num, bin_num], range=[[10, 22], [10, 22]], cmap=cmap)
    plt.plot(good_cat_stars['%s' % magcol][idx_mass],
             predict_y_for(good_cat_stars['%s' % magcol][idx_mass], m2, b2), c='b')
    plt.colorbar(label='Density')
    box = dict(facecolor='white')
    plt.text(11, 18, txt, fontsize=12, bbox=box)
    plt.savefig('%s_C%s_WLS_fit_hist_plot_%s%s' % (survey, chip, num, end_name), dpi=300)
    plt.clf()

    print('Saved WLS fit plots to dir!')

    # WLS 3 sig hist density plot

    fig, ax2 = plt.subplots(1, 1, figsize=(9, 8))
    fig.suptitle(f'{survey} vs PRIME w/ Weighted Fit - {sigma} Sigma Clip - Density Histogram')

    hist = ax2.hist2d(
        x=good_cat_stars['%s' % magcol][idx_mass][~clipped.mask],
        y=cleanPSFsources[f'{band}MAG_{MAGTYPES[magtype]}'][idx_image][~clipped.mask],
        bins=[bin_num, bin_num], range=[[10, 22], [10, 22]],
        cmap=cmap
    )

    line = ax2.plot(good_cat_stars['%s' % magcol][idx_mass][~clipped.mask],
                        predict_y_for(good_cat_stars['%s' % magcol][idx_mass][~clipped.mask], m, b),
                        c='r')

    cbar2 = fig.colorbar(hist[3], ax=ax2, pad=0.03)
    cbar2.set_label(f'{MAGTYPES[magtype]} Density')

    ax2.grid()
    ax2.set_xlim(10, 22)
    ax2.set_ylim(10, 22)
    ax2.yaxis.set_tick_params(labelleft=True)
    ax2.set_title(f"{survey} vs PRIME {MAGTYPES[magtype]} Plot w/ WLS fit line")
    ax2.set_xlabel(f"{survey} {band} Mags")
    ax2.set_ylabel(f"PRIME {band} Mags")

    ax2.text(11, 18, info, fontsize=12,
             bbox=dict(facecolor='white', edgecolor='black'))

    plt.savefig('%s_C%s_WLS_fit_3sig_hist_plot_%s%s' % (survey, chip, num, end_name), dpi=300)
    plt.close(fig)

    print('Saved WLS 3 sig fit plots to dir!')

    # Limiting Mag Plot

    def lim_mag_calc(bin_vals, all_mags):
        idxs = np.digitize(all_mags, bins=bin_vals)

        indices = {i: [] for i in range(len(bin_vals))}
        for idx, value in enumerate(idxs):
            indices[value].append(idx)
        sorted_indices_lists = list(indices.values())

        all_sources = []
        for i in sorted_indices_lists:
            number = len(i)
            all_sources.append(number)

        # plotting
        idx_arr = np.where(np.isclose(bin_vals, 12.5))
        min_x_idx = idx_arr[0][0]  # avoid saturated <12.5 mag sources from affecting maximum

        filtered_sources = all_sources[min_x_idx:]
        idxmax = filtered_sources.index(max(filtered_sources)) + min_x_idx
        split_sources = all_sources[idxmax:]

        halfmax = max(filtered_sources) / 2
        halfmaxpt = list(max(enumerate(split_sources), key=lambda x: -abs(halfmax - x[1])))
        halfmaxpt = [halfmaxpt[0] + idxmax, halfmaxpt[1]]

        limmag = round(bin_vals[halfmaxpt[0]], 1)

        return all_sources, halfmax, limmag

    max_x = data.shape[0]
    max_y = data.shape[1]
    all_mags_all = PSFsources[(PSFsources[f'{band}MAG_{MAGTYPES[magtype]}'] < 25) &
                              (PSFsources['XWIN_IMAGE'] < (max_x - crop)) & (PSFsources['XWIN_IMAGE'] > crop) &
                              (PSFsources['YWIN_IMAGE'] < (max_y) - crop) & (PSFsources['YWIN_IMAGE'] > crop)
                              ]
    all_mags = all_mags_all[f'{band}MAG_{MAGTYPES[magtype]}']

    bin_vals = np.array(np.arange(12, 25.5, 0.1))

    all_sources, halfmax, limmag = lim_mag_calc(bin_vals=bin_vals, all_mags=all_mags)

    print(f'{MAGTYPES[magtype]} Lim Mag = ', limmag)

    plt.figure(8, figsize=(24, 8))

    lmdata = plt.bar(bin_vals, height=all_sources, width=0.1, align='edge', color=lm_colors[0], edgecolor='black',
                       alpha=1, label=f'{MAGTYPES[magtype]} Binned Sources')
    lmhalf = plt.axhline(halfmax, linestyle='--', color=lm_colors[1],
                           label=f'{MAGTYPES[magtype]} Half Max = %s' % round(halfmax, 1))
    limmag_line = plt.axvline(limmag, color=lm_colors[1], linewidth=2,
                             label=f'{MAGTYPES[magtype]} Limiting Mag = %s' % round(limmag, 1))

    xticks = np.arange(12, 25.5, 0.5)
    plt.xticks(xticks, fontsize=10)
    plt.grid()
    plt.yscale('log')
    plt.title('PRIME Limiting Mag Plot')
    plt.ylabel('Number of Sources')
    plt.xlabel('%s Magnitude' % band)

    handles = [lmdata, lmhalf, limmag_line]
    labels = [h.get_label() for h in handles]
    plt.legend(handles, labels, fontsize=15, loc='upper right')
    plt.savefig('%s_C%s_lim_mag_plot_%s%s' % (survey, chip, num, end_name), dpi=300)
    print('Saved lim mag plot to dir!')
    plt.clf()

    plt.close('all')

    # source location regions
    newtext = open('PRIME_%s_C%s_crsmtched_srcs_%s.reg' % (survey, chip, num), 'w+')
    newtext.write('fk5')
    for a, d, rad in zip(cleanPSFsources['ALPHA_J2000'][idx_image], cleanPSFsources['DELTA_J2000'][idx_image],
                         cleanPSFsources['FLUX_RADIUS'][idx_image]):
        newtext.write(f'\ncircle({a}, {d}, {rad}") # color=red')

    newtext_all = open('PRIME_%s_C%s_all_srcs_%s.reg' % (survey, chip, num), 'w+')
    newtext_all.write('fk5')
    for a, d, rad in zip(PSFsources['ALPHA_J2000'], PSFsources['DELTA_J2000'],
                         PSFsources['FLUX_RADIUS']):
        newtext_all.write(f'\ncircle({a}, {d}, {rad}") # color=green')

    print('Source location reg files saved!')

    print('Writing relevant plot info to image header...')
    with open_fits_robust(imageName) as hdul:
        hdr = hdul[0].header
        hdr.set(f'{MAGTYPES[magtype]}_M', m, 'WLS %s sig fit slope' % sigma, after='SURVEY')
        hdr.set(f'E_{MAGTYPES[magtype]}_M', m_err, 'Error in WLS %s sig fit slope' % sigma, after=f'{MAGTYPES[magtype]}_M')
        hdr.set(f'{MAGTYPES[magtype]}_B', b, 'WLS %s sig fit intercept' % sigma, after=f'E_{MAGTYPES[magtype]}_M')
        hdr.set(f'E_{MAGTYPES[magtype]}_B', b_err, 'Error in WLS %s sig fit intercept' % sigma, after=f'{MAGTYPES[magtype]}_B')
        hdr.set(f'LM_{MAGTYPES[magtype]}', limmag, 'Source Histogram FWHM Limiting Mag', after=f'E_{MAGTYPES[magtype]}_B')

    return m, b, round(3 * b_err, 4)


#%% old plots function, MOSTLY DEPRECIATED (only used w/ -no_int_cal flag)

def photometry_plots(cleanPSFsources, PSFsources, data, imageName, survey, band, good_cat_stars, idx_psfmass, idx_psfimage,
                     psfweights_noclip, psf_clipped, sigma, aperweights_noclip, aper_clipped_all, autoweights_noclip, auto_clipped,
                     magtype=defaults['magtype']):

    # TODO temporary disabling of in progress aper mag plotting
    aperweights_noclip = []
    aper_clipped_all = []

    # appropriate mag column
    colnames = good_cat_stars.colnames
    if len(colnames) > 6:
        magcol = f'{band}MAG_{magtype}'
        magerrcol = f'{band}MAG_{magtype}'
    else:
        magcol = colnames[2]
        magerrcol = colnames[3]

    # aperture sizes
    if len(aperweights_noclip) > 0:
        apers_str = fits.getheader(imageName)['APERS']
        aper_arr = apers_str.split(',')

    try:
        chip = fits.getval(imageName,'CHIP')
    except KeyError:
        chip = 0
    if not imageName.startswith(f'coadd.Open-{band}') and not imageName.endswith('.fits'):
        num = 'img'
    else:
        num = imageName[-16:-8]

    # PSF-specific naming for parallel running
    end_name = end_name_gen('png')

    def predict_y_for(x, m, b):
        return m * x + b

    # PHOTOMETRIC FIT LINE MODELS

    model2, model_sig, model_sig_resid, model_auto, model_auto_resid, aper_model_sigs = (
        photometric_fit_calc(cleanPSFsources, band, good_cat_stars, idx_psfmass, idx_psfimage,
                         psfweights_noclip, psf_clipped, sigma, aperweights_noclip, aper_clipped_all,
                         autoweights_noclip, auto_clipped, with_plots=True))

    if model_auto:
        m_auto = model_auto.params[1]
        m_autoerr = model_auto.bse[1]
        b_auto = model_auto.params[0]
        b_autoerr = model_auto.bse[0]
        rsquare_auto = model_auto.rsquared
        rss_auto = model_auto.ssr
        model_auto_resid = model_auto.resid

        m2 = model2.params[1]
        m2err = model2.bse[1]
        b2 = model2.params[0]
        b2err = model2.bse[0]
        rsquare2 = model2.rsquared
        rss2 = model2.ssr
        avg2 = np.average(model2.resid, weights=autoweights_noclip)
        var2 = np.average((model2.resid - avg2) ** 2, weights=autoweights_noclip)

    if model_sig:
        m_sig = model_sig.params[1]
        m_sigerr = model_sig.bse[1]
        b_sig = model_sig.params[0]
        b_sigerr = model_sig.bse[0]
        rsquare_sig = model_sig.rsquared
        rss_sig = model_sig.ssr
        model_sig_resid = model_sig.resid

        m2 = model2.params[1]
        m2err = model2.bse[1]
        b2 = model2.params[0]
        b2err = model2.bse[0]
        rsquare2 = model2.rsquared
        rss2 = model2.ssr
        avg2 = np.average(model2.resid, weights=psfweights_noclip)
        var2 = np.average((model2.resid - avg2) ** 2, weights=psfweights_noclip)

    # BINNED RESIDUAL STATISTICS

    # auto aperture residual fit 3 sig clip - bin errors and stats
    mags = range(10, 22)
    x_arr = np.arange(10 + 0.5, 22 + 0.5, 1)

    auto_res_errs = []
    auto_res_means = []
    auto_res_ranges = []
    auto_res_skews = []
    auto_res_nums = []
    for i in mags:
        auto_mask = ((cleanPSFsources['%sMAG_AUTO' % band][idx_psfimage][~auto_clipped.mask] > i) &
                (cleanPSFsources['%sMAG_AUTO' % band][idx_psfimage][~auto_clipped.mask] < i+1))
        auto_mask = model_auto_resid[auto_mask]
        resauto_err = np.std(auto_mask)   # spread
        resauto_range = '%s - %s' % (i, i + 1)
        resauto_mean = np.nanmean(auto_mask)  # mean
        resauto_skew = skew(auto_mask, bias=False)    # skew
        resauto_num = len(auto_mask)
        if ma.is_masked(resauto_err):
            resauto_err = 0
        auto_res_errs.append(resauto_err)
        auto_res_means.append(resauto_mean)
        auto_res_ranges.append(resauto_range)
        auto_res_skews.append(resauto_skew)
        auto_res_nums.append(resauto_num)
    auto_res_errs = np.nan_to_num(np.array(auto_res_errs))
    auto_res_means = np.nan_to_num(np.array(auto_res_means))
    auto_res_skews = np.nan_to_num(np.array(auto_res_skews))
    auto_res_errs_min = np.min((auto_res_errs[auto_res_errs != 0]))
    auto_res_errs_max = np.max(auto_res_errs)

    # residual fit 3 sig clip - bin errors and stats
    if model_sig:
        res_errs = []
        res_means = []
        res_ranges = []
        res_skews = []
        res_nums = []
        for i in mags:
            mask = ((cleanPSFsources['%sMAG_PSF' % band][idx_psfimage][~psf_clipped.mask] > i) &
                    (cleanPSFsources['%sMAG_PSF' % band][idx_psfimage][~psf_clipped.mask] < i+1))
            res_mask = model_sig_resid[mask]
            # avgsig = np.average(res_mask, weights=psfweights_noclip[~psf_clipped.mask][mask])
            # varsig = np.average((res_mask - avgsig)**2, weights=psfweights_noclip[~psf_clipped.mask][mask])
            # ressig_err = np.sqrt(varsig)    # spread
            ressig_err = np.std(res_mask)   # spread
            ressig_range = '%s - %s' % (i, i + 1)
            ressig_mean = np.nanmean(res_mask)  # mean
            ressig_skew = skew(res_mask, bias=False)    # skew
            ressig_num = len(res_mask)
            if ma.is_masked(ressig_err):
                ressig_err = 0
            res_errs.append(ressig_err)
            res_means.append(ressig_mean)
            res_ranges.append(ressig_range)
            res_skews.append(ressig_skew)
            res_nums.append(ressig_num)
        res_errs = np.nan_to_num(np.array(res_errs))
        res_means = np.nan_to_num(np.array(res_means))
        res_skews = np.nan_to_num(np.array(res_skews))
        res_errs_min = np.min((res_errs[res_errs != 0]))
        res_errs_max = np.max(res_errs)

        # bin errors / stats table
        bintable = Table()
        bintable['Bin Range (mag)'] = [res_ranges, auto_res_ranges]
        bintable['Mean (mag)'] = [res_means, auto_res_means]
        bintable['Spread (mag)'] = [res_errs, auto_res_errs]
        bintable['Skew'] = [res_skews, auto_res_skews]
        bintable['Source Number'] = [res_nums, auto_res_nums]

        bintable.meta['Description'] = ('Table of data on residuals, binned by mag. Each column is a vector w/ the 1st idx '
                                        'being the PSF fit data and the 2nd idx being auto aperture photom fit data')

        bintable.write('Resid_%s-sig_Data_%s_C%s_%s%s' % (sigma, band, chip, survey, end_name_gen()), overwrite=True)
        print('Residual data table for PSF and auto photometry written!')

    else:
        bintable = Table()
        bintable['Bin Range (mag)'] = auto_res_ranges
        bintable['Mean (mag)'] = auto_res_means
        bintable['Spread (mag)'] = auto_res_errs
        bintable['Skew'] = auto_res_skews
        bintable['Source Number'] = auto_res_nums

        bintable.meta['Description'] = 'Table of data on residuals for PSF fit data, binned by mag'
        bintable.write('Resid_%s-sig_Data_%s_C%s_%s.ecsv' % (sigma, band, chip, survey), overwrite=True)
        print('Residual data table for %s sigma clip written!' % sigma)

    # aper photom residual fit sig clip - bin errors and stats
    if len(aperweights_noclip) > 0:
        aper_res_errs_all = []
        aper_res_means_all = []
        aper_res_ranges_all = []
        aper_res_skews_all = []
        aper_res_nums_all = []
        aper_res_max_min_all = []

        for idx, (aperclipped, apermodel) in enumerate(zip(aper_clipped_all, aper_model_sigs)):
            ap_res_errs = []
            ap_res_means = []
            ap_res_ranges = []
            ap_res_skews = []
            ap_res_nums = []
            for mag in mags:
                ap_mask = ((cleanPSFsources['%sMAG_APER' % band][idx_psfimage][~aperclipped.mask] > mag) &
                        (cleanPSFsources['%sMAG_PSF' % band][idx_psfimage][~aperclipped.mask] < mag+1))
                ap_res_mask = apermodel.resid[ap_mask]
                ap_ressig_err = np.std(ap_res_mask)   # spread
                ap_ressig_range = '%s - %s' % (mag, mag + 1)
                ap_ressig_mean = np.nanmean(ap_res_mask)  # mean
                ap_ressig_skew = skew(ap_res_mask, bias=False)    # skew
                ap_ressig_num = len(ap_res_mask)
                if ma.is_masked(ap_ressig_err):
                    ap_ressig_err = 0
                ap_res_errs.append(ap_ressig_err)
                ap_res_means.append(ap_ressig_mean)
                ap_res_ranges.append(ap_ressig_range)
                ap_res_skews.append(ap_ressig_skew)
                ap_res_nums.append(ap_ressig_num)
            ap_res_errs = np.nan_to_num(np.array(ap_res_errs))
            ap_res_means = np.nan_to_num(np.array(ap_res_means))
            ap_res_skews = np.nan_to_num(np.array(ap_res_skews))
            ap_res_errs_min = np.min((ap_res_errs[ap_res_errs != 0]))
            ap_res_errs_max = np.max(ap_res_errs)
            ap_res_max_min = (ap_res_errs_max, ap_res_errs_min)

            aper_res_errs_all.append(ap_res_errs)
            aper_res_means_all.append(ap_res_means)
            aper_res_ranges_all.append(ap_res_ranges)
            aper_res_skews_all.append(ap_res_skews)
            aper_res_nums_all.append(ap_res_nums)
            aper_res_max_min_all.append(ap_res_max_min)

    # ALL PLOTS

    plt.close('all')
    # mag comparison plot
    plt.figure(1, figsize=(8, 8))
    plt.plot(cleanPSFsources[f'{band}MAG_AUTO'][idx_psfimage], good_cat_stars['%s' % magcol][idx_psfmass],
             'r.', markersize=14, markeredgecolor='black')
    plt.xlim(10, 22)
    plt.ylim(10, 22)
    plt.title('PRIME Mags vs %s Mags' % survey)
    plt.xlabel('PRIME %s Mags' % band, fontsize=15)
    plt.ylabel('%s %s Mags' % (survey, band), fontsize=15)
    plt.grid()
    # plt.savefig('%s_C%s_mag_comp_plot_%s.png' % (survey, chip, num))
    plt.clf()
    # print('Saved mag comparison plot to dir!')

    # PRIME flux vs catalog AB mag for crossmatches
    x_flx_lin = cleanPSFsources[f'AUTO_FLUX_DENSITY'][idx_psfimage][~auto_clipped.mask]
    x_flx = np.log(x_flx_lin)
    y_flx = good_cat_stars['%s' % magcol][idx_psfmass][~auto_clipped.mask]
    x_const_flx = sm.add_constant(x_flx)
    model_flx = sm.WLS(y_flx, x_const_flx, weights=autoweights_noclip[~auto_clipped.mask]).fit()
    # print(model_flx.params)
    m_flx = model_flx.params[1]
    m_flxerr = model_flx.bse[1]
    b_flx = model_flx.params[0]
    b_flxerr = model_flx.bse[0]

    x_grid = np.linspace(x_flx_lin.min(), x_flx_lin.max(), 100)
    log_x_grid = np.log(x_grid)
    log_x_grid_const = sm.add_constant(log_x_grid)
    y_pred = model_flx.predict(log_x_grid_const)

    plt.figure(9, figsize=(8, 8))
    plt.plot(cleanPSFsources['AUTO_FLUX_DENSITY'][idx_psfimage][~auto_clipped.mask],
             good_cat_stars['%s' % magcol][idx_psfmass][~auto_clipped.mask],
             'r.', markersize=14, markeredgecolor='black')
    plt.plot(x_grid, y_pred, c='b')

    flx_txt = ('slope = %.4f' % m_flx + '\nslope err = %.4f' % m_flxerr +
           '\nint = %.4f' % b_flx + '\nint err = %.4f' % b_flxerr)

    plt.xlim(10, 50000)
    plt.ylim(12, 20.5)
    plt.title('PRIME Flux Density vs %s AB mag - %s Sigma Clip' % (survey, sigma))
    plt.xlabel(r'PRIME Flux Density ($\mu$Jy)', fontsize=15)
    plt.ylabel('%s %s AB Mags' % (survey, band), fontsize=15)
    plt.grid()
    plt.xscale('log')
    flx_box = dict(facecolor='white')
    plt.text(1000, 19, flx_txt, fontsize=12, bbox=flx_box)
    plt.savefig('%s_C%s_flux_mag_plot_sig_%s%s' % (survey, chip, num, end_name))
    plt.clf()
    print('Saved flux v. mag plot to dir!')


    # residual plot - y int forced to zero
    """
    plt.figure(2, figsize=(8, 6))
    plt.scatter(cleanPSFsources['%sMAG_PSF' % band][idx_psfimage], model.resid, color='red')
    plt.ylim(-1.5, 1.5)
    plt.xlim(10,21)
    plt.title('PRIME vs %s Residuals' % survey)
    plt.ylabel('Residuals')
    plt.xlabel('Mags')
    plt.axhline(y=res_err, color='tab:orange', linestyle='--', linewidth=1)
    plt.axhline(y=res_err * 2, color='green', linestyle='--', linewidth=1)
    plt.axhline(y=-res_err, color='tab:orange', linestyle='--', linewidth=1)
    plt.axhline(y=-res_err * 2, color='green', linestyle='--', linewidth=1)
    plt.axhline(y=0, color='black', linestyle='--', linewidth=1)
    plt.legend(['Residuals', r'1 $\sigma$ = %.3f' % res_err, r'2 $\sigma$ = %.3f' % (2*res_err)], loc='lower left')
    info = ('eqn: y = mx'+'\nslope = %.5f +/- %.5f' % (m,merr))+('\nR$^{2}$ = %.3f' % rsquare)+('\nRSS = %d' % rss)
    #+('\nintercept = %.3f +/- %.3f' % (b,berr))
    plt.text(15,-1.25,info,bbox=dict(facecolor='white',edgecolor='black',alpha=1,pad=5.0))
    plt.savefig('%s_C%s_residual_plot_%s.png' % (survey,chip,num),dpi=300)
    print('Saved residual plot to dir!')
    """

    # res plot, y int include, PSF and aperture fits

    if len(aperweights_noclip) > 0:
        plt.figure(2, figsize=(8, 6))
        psf_sc = plt.scatter(good_cat_stars['%s' % magcol][idx_psfmass][~psf_clipped.mask], model_sig_resid, color='red', alpha=0.5,
                    label='PSF Residuals')
        plt.ylim(-1, 1)
        plt.xlim(10, 21)
        plt.title('%s Residuals - %s Sigma Clip' % (survey, sigma))
        plt.ylabel('Residuals')
        plt.xlabel('%s %s Mags' % (survey, band))
        # plt.axhline(y=ressig_err, color='blue', linestyle='--', linewidth=1)
        # plt.axhline(y=ressig_err * 2, color='green', linestyle='--', linewidth=1)
        # plt.axhline(y=-ressig_err, color='blue', linestyle='--', linewidth=1)
        # plt.axhline(y=-ressig_err * 2, color='green', linestyle='--', linewidth=1)
        plt.scatter(x_arr, res_errs, marker='_', s=1625, c='black')
        plt.scatter(x_arr, -res_errs, marker='_', s=1625, c='black')
        plt.axhline(y=0, color='black', linestyle='--', linewidth=1)
        info2 = ('PSF eqn: y = mx+b' + '\nslope = %.4f +/- %.4f' % (m_sig, m_sigerr)) + (
                    '\nintercept = %.3f +/- %.3f' % (b_sig, b_sigerr)) + ('\nR$^{2}$ = %.3f' % rsquare_sig) + ('\nRSS = %d' % rss_sig)
        plt.text(15, -0.9, info2, fontsize=9, bbox=dict(facecolor='white', edgecolor='black', pad=5.0))

        scatters = []
        color_arr = ['blue', 'green', 'magenta', 'yellow', 'tab:orange']
        for aperclipped, apermodel, aperminmaxes, apererrs, apersize, color in (
                zip(aper_clipped_all, aper_model_sigs, aper_res_max_min_all, aper_res_errs_all, aper_arr, color_arr)):
            sc = plt.scatter(good_cat_stars['%s' % magcol][idx_psfmass][~aperclipped.mask], apermodel.resid, color=color,
                        alpha=0.4, label=f'{apersize * 0.498}" Aperture Residuals')
            scatters.append(sc)
            plt.scatter(x_arr, apererrs, marker='_', s=1625, c='black', alpha=0.4)
            plt.scatter(x_arr, -apererrs, marker='_', s=1625, c='black', alpha=0.4)

        handles = scatters + [psf_sc]
        labels = [h.get_label() for h in handles]
        plt.legend(handles, labels, loc='lower left',
                   markerscale=0.5)

        plt.savefig('%s_C%s_residual_plot_all_%s%s' % (survey, chip, num, end_name), dpi=300)
        plt.clf()

    # res plot y int, histogram
    if len(idx_psfimage) >= 75000:
        bin_num_int = round(len(idx_psfimage) / 500)
    elif 50000 <= len(idx_psfimage) <= 75000:
        bin_num_int = round(len(idx_psfimage) / 400)
    elif 25000 <= len(idx_psfimage) <= 50000:
        bin_num_int = round(len(idx_psfimage) / 300)
    elif 5000 <= len(idx_psfimage) <= 25000:
        bin_num_int = round(len(idx_psfimage) / 75)
    elif 1000 <= len(idx_psfimage) <= 5000:
        bin_num_int = round(len(idx_psfimage) / 50)
    elif len(idx_psfimage) <= 1000:
        bin_num_int = 50

    alpha = 0.7 if model_auto else None

    # PSF PANEL
    if model_sig:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 8), sharey=True)
        fig.subplots_adjust(wspace=0.12)
        fig.suptitle(f'{survey} Residuals - {sigma} Sigma Clip - Density Histogram')

        psfhist = ax1.hist2d(x=good_cat_stars['%s' % magcol][idx_psfmass][~psf_clipped.mask],y=model_sig_resid,
            bins=[bin_num_int, bin_num_int], range=[[10, 21], [-1, 1]],
            cmap='gist_heat_r'
        )
        cbar1 = fig.colorbar(psfhist[3], ax=ax1, pad=0.03)
        cbar1.set_label('PSF Density')

        psfsig = ax1.scatter(x_arr, res_errs, marker='_', s=1625, c='blue',
                             label=r'PSF photom 1 $\sigma$ range = [%.3f - %.3f]' % (res_errs_min, res_errs_max))
        ax1.scatter(x_arr, -res_errs, marker='_', s=1625, c='blue')

        ax1.axhline(y=0, color='black', linestyle='--', linewidth=1)
        ax1.set_xlim(10, 21)
        ax1.set_ylim(-1, 1)
        ax1.set_title(f"PSF Fit Residuals")
        ax1.set_xlabel(f"{survey} {band} Mags")
        ax1.set_ylabel("Residuals")

        infohist = (
                'PSF fit'
                + '\nslope = %.4f +/- %.4f' % (m_sig, m_sigerr)
                + '\nintercept = %.3f +/- %.3f' % (b_sig, b_sigerr)
                + '\nR$^{2}$ = %.3f' % rsquare_sig
                + '\nRSS = %d' % rss_sig
                + '\nn_sources = %i' % len(cleanPSFsources['PSF_FLUX_DENSITY'][idx_psfimage][~psf_clipped.mask])
        )

        ax1.text(10.5, 0.5, infohist, fontsize=9,
                 bbox=dict(facecolor='white', edgecolor='black', pad=5.0))

        ax1.legend([psfsig], [psfsig.get_label()], loc='lower left', markerscale=0.5)

    else:
        fig, ax2 = plt.subplots(1, 1, figsize=(9, 8))
        fig.suptitle(f'{survey} Residuals - {sigma} Sigma Clip - Density Histogram')

    # AUTO PANEL
    autohist = ax2.hist2d(
        x=good_cat_stars['%s' % magcol][idx_psfmass][~auto_clipped.mask],y=model_auto_resid,
        bins=[bin_num_int, bin_num_int], range=[[10, 21], [-1, 1]],
        cmap='gist_earth_r'
    )

    cbar2 = fig.colorbar(autohist[3], ax=ax2, pad=0.03)
    cbar2.set_label('Auto Ap. Density')

    autosig = ax2.scatter(x_arr, auto_res_errs, marker='_', s=1625, c='black',
                          label=r'Auto ap. photom 1 $\sigma$ range = [%.3f - %.3f]' % (
                          auto_res_errs_min, auto_res_errs_max))
    ax2.scatter(x_arr, -auto_res_errs, marker='_', s=1625, c='black')

    ax2.axhline(y=0, color='black', linestyle='--', linewidth=1)
    ax2.set_xlim(10, 21)
    ax2.set_ylim(-1, 1)
    ax2.yaxis.set_tick_params(labelleft=True)
    ax2.set_title(f"Auto Aperture Fit Residuals")
    ax2.set_xlabel(f"{survey} {band} Mags")
    ax2.set_ylabel("Residuals")

    infoauto = (
            'Auto Ap. fit'
            + '\nslope = %.4f +/- %.4f' % (m_auto, m_autoerr)
            + '\nintercept = %.3f +/- %.3f' % (b_auto, b_autoerr)
            + '\nR$^{2}$ = %.3f' % rsquare_auto
            + '\nRSS = %d' % rss_auto
            + '\nn_sources = %i' % len(cleanPSFsources['AUTO_FLUX_DENSITY'][idx_psfimage][~auto_clipped.mask])
    )

    ax2.text(10.5, 0.5, infoauto, fontsize=9,
             bbox=dict(facecolor='white', edgecolor='black', pad=5.0))

    ax2.legend([autosig], [autosig.get_label()], loc='lower left', markerscale=0.5)

    plt.savefig('%s_C%s_residual_plot_int_hist_%s%s' % (survey, chip, num, end_name), dpi=300)
    plt.close(fig)

    print('Saved y-int residual plots to dir!')

    # WLS fit line over data plot

    txt = ('slope = %.4f' % m2 + '\nslope err = %.4f' % m2err + '\nint = %.4f' % b2 + '\nint err = %.4f' % b2err +
           '\nn_sources = %i' % len(cleanPSFsources[idx_psfimage]))
    #
    # plt.figure(4, figsize=(8, 8))
    # plt.xlim(10, 22)
    # plt.ylim(10, 22)
    # plt.title('%s vs PRIME w/ Weighted Fit' % survey)
    # plt.grid()
    # plt.ylabel('PRIME %s Mags' % band, fontsize=15)
    # plt.xlabel('%s %s Mags' % (survey, band), fontsize=15)
    # plt.scatter(good_cat_stars['%s' % magcol][idx_psfmass], cleanPSFsources['%sMAG_AUTO' % band][idx_psfimage])
    # plt.plot(good_cat_stars['%s' % magcol][idx_psfmass],
    #          predict_y_for(good_cat_stars['%s' % magcol][idx_psfmass], m2, b2), c='r')
    # box = dict(facecolor='white')
    # plt.text(11, 18, txt, fontsize=12, bbox=box)
    # plt.savefig('%s_C%s_WLS_fit_plot_%s.png' % (survey, chip, num), dpi=300)

    # WLS hist density plot
    if len(idx_psfimage) >= 75000:
        bin_num = round(len(idx_psfimage) / 500)
    elif 50000 <= len(idx_psfimage) <= 75000:
        bin_num = round(len(idx_psfimage) / 350)
    elif 25000 <= len(idx_psfimage) <= 50000:
        bin_num = round(len(idx_psfimage) / 150)
    elif 5000 <= len(idx_psfimage) <= 25000:
        bin_num = round(len(idx_psfimage) / 50)
    elif 1000 <= len(idx_psfimage) <= 5000:
        bin_num = round(len(idx_psfimage) / 20)
    elif len(idx_psfimage) <= 1000:
        bin_num = 100

    plt.figure(5, figsize=(10, 8))
    plt.clf()
    plt.xlim(10, 22)
    plt.ylim(10, 22)
    plt.title('%s vs PRIME w/ Weighted Fit - Density Histogram' % survey)
    plt.grid()
    plt.ylabel('PRIME %s Mags' % band, fontsize=15)
    plt.xlabel('%s %s Mags' % (survey, band), fontsize=15)
    plt.hist2d(x=good_cat_stars['%s' % magcol][idx_psfmass], y=cleanPSFsources['%sMAG_AUTO' % band][idx_psfimage],
               bins=[bin_num, bin_num], range=[[10, 22],[10, 22]], cmap='gist_heat_r')
    plt.plot(good_cat_stars['%s' % magcol][idx_psfmass],
             predict_y_for(good_cat_stars['%s' % magcol][idx_psfmass], m2, b2), c='b')
    plt.colorbar(label='Density')
    box = dict(facecolor='white')
    plt.text(11, 18, txt, fontsize=12, bbox=box)
    plt.savefig('%s_C%s_WLS_fit_hist_plot_%s%s' % (survey, chip, num, end_name), dpi=300)
    plt.clf()

    print('Saved WLS fit plots to dir!')

    # WLS fit for 3 sig clip of data
    # sigtxt = ('slope = %.4f' % m_sig + '\nslope err = %.4f' % m_sigerr + '\nint = %.4f' % b_sig +
    #           '\nint err = %.4f' % b_sigerr + '\nn_sources = %i' % len(cleanPSFsources[idx_psfimage][~psf_clipped.mask]))
    #
    # plt.figure(6, figsize=(8, 8))
    # plt.clf()
    # plt.xlim(10, 22)
    # plt.ylim(10, 22)
    # plt.title('%s vs PRIME w/ Weighted Fit - %s Sigma Clip' % (survey, sigma))
    # plt.grid()
    # plt.ylabel('PRIME %s Mags' % band, fontsize=15)
    # plt.xlabel('%s %s Mags' % (survey, band), fontsize=15)
    # plt.scatter(good_cat_stars['%s' % magcol][idx_psfmass][~psf_clipped.mask],
    #             cleanPSFsources['%sMAG_AUTO' % band][idx_psfimage][~psf_clipped.mask])
    # plt.plot(good_cat_stars['%s' % magcol][idx_psfmass][~psf_clipped.mask],
    #          predict_y_for(good_cat_stars['%s' % magcol][idx_psfmass][~psf_clipped.mask], m_sig, b_sig), c='r')
    # box = dict(facecolor='white')
    # plt.text(11, 18, sigtxt, fontsize=12, bbox=box)
    # plt.savefig('%s_C%s_WLS_fit_3sig_plot_%s.png' % (survey, chip, num), dpi=300)

    # WLS 3 sig hist density plot

    # PSF PANEL
    if model_sig:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 8), sharey=True)
        fig.subplots_adjust(wspace=0.12)
        fig.suptitle(f'{survey} vs PRIME w/ Weighted Fit - {sigma} Sigma Clip - Density Histogram')

        psfhist = ax1.hist2d(x=good_cat_stars['%s' % magcol][idx_psfmass][~psf_clipped.mask],
                   y=cleanPSFsources['%sMAG_PSF' % band][idx_psfimage][~psf_clipped.mask],
                   bins=[bin_num, bin_num], range=[[10, 22],[10, 22]], cmap='gist_heat_r')

        psfline = ax1.plot(good_cat_stars['%s' % magcol][idx_psfmass][~psf_clipped.mask],
                 predict_y_for(good_cat_stars['%s' % magcol][idx_psfmass][~psf_clipped.mask], m_sig, b_sig), c='b')

        cbar1 = fig.colorbar(psfhist[3], ax=ax1, pad=0.03)
        cbar1.set_label('PSF Density')

        ax1.grid()
        ax1.set_xlim(10, 22)
        ax1.set_ylim(10, 22)
        ax1.set_title(f"{survey} vs PRIME PSF Mag Plot w/ WLS fit line")
        ax1.set_xlabel(f"{survey} {band} Mags")
        ax1.set_ylabel(f"PRIME {band} Mags")

        ax1.text(11, 18, infohist, fontsize=12,
                 bbox=dict(facecolor='white', edgecolor='black'))

    else:
        fig, ax2 = plt.subplots(1, 1, figsize=(9, 8))
        fig.suptitle(f'{survey} vs PRIME w/ Weighted Fit - {sigma} Sigma Clip - Density Histogram')

    # AUTO PANEL
    autohist = ax2.hist2d(
        x=good_cat_stars['%s' % magcol][idx_psfmass][~auto_clipped.mask],
        y=cleanPSFsources['%sMAG_AUTO' % band][idx_psfimage][~auto_clipped.mask],
        bins=[bin_num, bin_num], range=[[10, 22],[10, 22]],
        cmap='gist_earth_r'
    )

    autoline = ax2.plot(good_cat_stars['%s' % magcol][idx_psfmass][~auto_clipped.mask],
                       predict_y_for(good_cat_stars['%s' % magcol][idx_psfmass][~auto_clipped.mask], m_auto, b_auto),
                       c='r')

    cbar2 = fig.colorbar(autohist[3], ax=ax2, pad=0.03)
    cbar2.set_label('Auto Ap. Density')

    ax2.grid()
    ax2.set_xlim(10, 22)
    ax2.set_ylim(10, 22)
    ax2.yaxis.set_tick_params(labelleft=True)
    ax2.set_title(f"{survey} vs PRIME Auto Ap. Plot w/ WLS fit line")
    ax2.set_xlabel(f"{survey} {band} Mags")
    ax2.set_ylabel(f"PRIME {band} Mags")

    ax2.text(11, 18, infoauto, fontsize=12,
             bbox=dict(facecolor='white', edgecolor='black'))

    plt.savefig('%s_C%s_WLS_fit_3sig_hist_plot_%s%s' % (survey, chip, num, end_name), dpi=300)
    plt.close(fig)

    print('Saved WLS 3 sig fit plots to dir!')

    # flux vs mag plot - histogram vers.
    plt.figure(10, figsize=(8, 8))
    plt.hist2d(x=cleanPSFsources['AUTO_FLUX_DENSITY'][idx_psfimage], y=good_cat_stars['%s' % magcol][idx_psfmass],
              bins=[bin_num, bin_num], range=[[100, 50000],[12, 20.5]], cmap='gist_heat_r')
    plt.colorbar(label='Density')
    plt.xlim(10, 50000)
    plt.ylim(12, 20.5)
    plt.title('PRIME Flux Density vs %s AB mag - Histogram' % survey)
    plt.xlabel(r'PRIME Flux Density ($\mu$Jy)', fontsize=15)
    plt.ylabel('%s %s AB Mags' % (survey, band), fontsize=15)
    plt.grid()
    plt.xscale('log')
    # plt.savefig('%s_C%s_flux_mag_hist_plot_%s.png' % (survey, chip, num))
    plt.clf()
    # print('Saved flux v. mag hist plot to dir!')

    # Limiting Mag Plot

    def lim_mag_calc(bin_vals, all_mags):
        idxs = np.digitize(all_mags, bins=bin_vals)

        indices = {i: [] for i in range(len(bin_vals))}
        for idx, value in enumerate(idxs):
            indices[value].append(idx)
        sorted_indices_lists = list(indices.values())

        all_sources = []
        for i in sorted_indices_lists:
            number = len(i)
            all_sources.append(number)

        # plotting
        idx_arr = np.where(np.isclose(bin_vals, 12.5))
        min_x_idx = idx_arr[0][0]   # avoid saturated <12.5 mag sources from affecting maximum

        filtered_sources = all_sources[min_x_idx:]
        idxmax = filtered_sources.index(max(filtered_sources)) + min_x_idx
        split_sources = all_sources[idxmax:]

        halfmax = max(filtered_sources) / 2
        halfmaxpt = list(max(enumerate(split_sources), key=lambda x: -abs(halfmax - x[1])))
        halfmaxpt = [halfmaxpt[0] + idxmax, halfmaxpt[1]]

        limmag = round(bin_vals[halfmaxpt[0]], 1)

        return all_sources, halfmax, limmag

    all_mags_all = PSFsources[PSFsources['%sMAG_AUTO' % band] < 25]
    all_mags = all_mags_all['%sMAG_AUTO' % band]

    bin_vals = np.array(np.arange(12, 25.5, 0.1))

    all_sources, halfmax, limmag = lim_mag_calc(bin_vals=bin_vals, all_mags=all_mags)

    print('Auto Ap. Lim Mag = ', limmag)

    plt.figure(8, figsize=(24, 8))

    autodata = plt.bar(bin_vals, height=all_sources, width=0.1, align='edge', color='blue', edgecolor='black',
                       alpha=1, label='Auto Ap. Binned Sources')
    autohalf = plt.axhline(halfmax, linestyle='--', color='black',
                           label='Auto ap. Half Max = %s' % round(halfmax, 1))
    autolimmag = plt.axvline(limmag, color='black', linewidth=2,
                             label='Auto ap. Limiting Mag = %s' % round(limmag, 1))

    psfhalf = psf_limmag = []
    if model_sig:
        psf_mags = PSFsources[PSFsources['%sMAG_PSF' % band] < 25]
        all_psf_mags = psf_mags['%sMAG_PSF' % band]

        all_psf_sources, psf_halfmax, psf_limmag = lim_mag_calc(bin_vals=bin_vals, all_mags=all_psf_mags)

        psfdata = plt.bar(bin_vals, height=all_psf_sources, width=0.1, align='edge', color='red', edgecolor='black',
                          alpha=alpha,
                          label='PSF Binned Sources')
        psfhalf = plt.axhline(psf_halfmax, linestyle='--', label='PSF Half Max = %s' % round(psf_halfmax, 1))
        psflimmag = plt.axvline(psf_limmag, color='b', linewidth=2, label='PSF Limiting Mag = %s' % round(psf_limmag, 1))
        print('PSF Lim Mag = ', psf_limmag)

    xticks = np.arange(12, 25.5, 0.5)
    plt.xticks(xticks, fontsize=10)
    plt.grid()
    plt.yscale('log')
    plt.title('PRIME Limiting Mag Plot')
    plt.ylabel('Number of Sources')
    plt.xlabel('%s Magnitude' % band)

    if psfhalf:
        handles = [psfdata, psfhalf, psflimmag, autodata, autohalf, autolimmag]
    else:
        handles = [autohalf, autolimmag]
    labels = [h.get_label() for h in handles]
    plt.legend(handles, labels, fontsize=15, loc='upper right')
    plt.savefig('%s_C%s_lim_mag_plot_%s%s' % (survey, chip, num, end_name), dpi=300)
    print('Saved lim mag plot to dir!')
    plt.clf()

    # Crossmatch location check plot
    # if not os.path.isfile('%s_C%s_source_check_plot_%s.png' % (survey, chip, num)):
    # mean, median, sigma_plot = sigma_clipped_stats(data)
    #
    # fig = plt.figure(figsize=(10, 10))
    # ax = fig.gca()
    #
    # im = ax.imshow(
    #     data,
    #     vmin=median - 1.5 * sigma_plot,
    #     vmax=median + 1.5 * sigma_plot,
    #     origin='lower'
    # )
    #
    # # Draw circles
    # circles = [
    #     plt.Circle(
    #         (cleanPSFsources['X_IMAGE'][idx_psfimage][i],
    #          cleanPSFsources['Y_IMAGE'][idx_psfimage][i]),
    #         radius=5,
    #         edgecolor='r',
    #         facecolor='None'
    #     ) for i in range(len(cleanPSFsources['X_IMAGE'][idx_psfimage]))
    # ]
    # for c in circles:
    #     ax.add_artist(c)
    #
    # cbar = fig.colorbar(im, ax=ax)
    # cbar.set_label("Pixel Value")
    #
    # plt.savefig('%s_C%s_source_check_plot_%s%s' % (survey, chip, num, end_name), dpi=150)
    # print('Saved source location check plot to dir!')
    # plt.clf()

    plt.close('all')

    # source location regions
    newtext = open('%s_C%s_crsmtched_srcs_%s%s' % (survey, chip, num, end_name_gen('reg')), 'w+')
    newtext.write('fk5')
    for a, d, rad in zip(cleanPSFsources['ALPHA_J2000'][idx_psfimage], cleanPSFsources['DELTA_J2000'][idx_psfimage],
                    cleanPSFsources['FLUX_RADIUS'][idx_psfimage]):
        newtext.write(f'\ncircle({a}, {d}, {rad}") # color=red')

    newtext = open('%s_C%s_all_srcs_%s.reg' % (survey, chip, num), 'w+')
    newtext.write('fk5')
    for a, d, rad in zip(PSFsources['ALPHA_J2000'], PSFsources['DELTA_J2000'],
                    PSFsources['FLUX_RADIUS']):
        newtext.write(f'\ncircle({a}, {d}, {rad}") # color=green')

    print('Source location reg files saved!')

    print('Writing relevant plot info to image header...')
    with open_fits_robust(imageName) as hdul:
        hdr = hdul[0].header
        hdr.set('AUTO_M', m_auto, 'WLS %s sig fit slope' % sigma, after='SURVEY')
        hdr.set('E_AUTO_M', m_autoerr, 'Error in WLS %s sig fit slope' % sigma, after='AUTO_M')
        hdr.set('AUTO_B', b_auto, 'WLS %s sig fit intercept' % sigma, after='E_AUTO_M')
        hdr.set('E_AUTO_B', b_autoerr, 'Error in WLS %s sig fit intercept' % sigma, after='AUTO_B')
        hdr.set('LM_AUTO', limmag, 'Source Histogram FWHM Limiting Mag', after='E_AUTO_B')
        if model_sig:
            try:
                hdr.set('PSF_M', m_sig, 'WLS %s sig fit slope' % sigma, after='auto_fit_m')
            except KeyError:
                hdr.set('PSF_M', m_sig, 'WLS %s sig fit slope' % sigma, after='SURVEY')
            hdr.set('E_PSF_M', m_sigerr, 'Error in WLS %s sig fit slope' % sigma, after='PSF_M')
            hdr.set('PSF_B', b_sig, 'WLS %s sig fit intercept' % sigma, after='E_PSF_M')
            hdr.set('E_PSF_B', b_sigerr, 'Error in WLS %s sig fit intercept' % sigma, after='PSF_B')
            hdr.set('LM_PSF', psf_limmag, 'Source Histogram FWHM Limiting Mag', after='E_PSF_B')

    # return m_sig, b_sig, round(3 * m_sigerr, 4)
    # TODO to run calibration based on auto aperture photom, uncomment line below, comment above
    return m_auto, b_auto, round(3 * b_autoerr, 4)
