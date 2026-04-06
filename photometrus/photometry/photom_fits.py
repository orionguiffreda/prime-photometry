"""
Photometric fit line calculation functions
"""

import numpy as np
import statsmodels.api as sm

from photometrus.settings import MAGTYPES

from photometrus.utils.defaults import PROCESSING_DEFAULTS as defaults


def single_fit_calc(cleanPSFsources, band, good_cat_stars, idx_psfmass, idx_psfimage,
                         weights_noclip, clipped, sigma, with_plots=False, magtype=defaults['magtype']):
    """Calculates WLS fit line used for determining the goodness of the photometry, for 1 magtype"""

    # appropriate mag column
    colnames = good_cat_stars.colnames
    if len(colnames) > 6:
        magcol = f'{band}MAG_{magtype}'
        magerrcol = f'{band}MAG_{magtype}'
    else:
        magcol = colnames[2]
        magerrcol = colnames[3]

    # sigma residual fit
    x = good_cat_stars['%s' % magcol][idx_psfmass][~clipped.mask]
    y = cleanPSFsources[f'{band}MAG_{MAGTYPES[magtype]}'][idx_psfimage][~clipped.mask]

    x_const = sm.add_constant(x)
    model = sm.WLS(y, x_const, weights=weights_noclip[~clipped.mask]).fit()

    x_nc = good_cat_stars['%s' % magcol][idx_psfmass]
    y_nc = cleanPSFsources[f'{band}MAG_{MAGTYPES[magtype]}'][idx_psfimage]
    x_const_nc = sm.add_constant(x_nc)
    model2 = sm.WLS(y_nc, x_const_nc, weights=weights_noclip).fit()

    m = model.params[1]
    m_err = model.bse[1]
    b = model.params[0]
    b_err = model.bse[0]
    model_resid = model.resid

    print('Num of crossmatched sources used in %s %s sig fit: %i'
          % (MAGTYPES[magtype], sigma, len(cleanPSFsources[f'{band}MAG_{MAGTYPES[magtype]}'][idx_psfimage][~clipped.mask])))
    print(' %s sig fit: slope = %.4f +/- %.4f' % (sigma, m, m_err))
    print(' %s sig fit: y-int = %.4f +/- %.4f' % (sigma, b, b_err))

    if with_plots:
        return model, model2
    else:
        full_b_err = round(3 * b_err, 4)
        if full_b_err >= 0.25:
            print(f'3 sig int. error larger than expected.. {full_b_err}')
            print('Reducing error (3 sig -> 0.5 sig) threshold significantly, should investigate image / results')
            return m, b, round(0.5 * b_err, 4)
        return m, b, full_b_err


#%% old fit calculation function, MOSTLY DEPRECIATED (only used w/ -no_int_cal flag)

def photometric_fit_calc(cleanPSFsources, band, good_cat_stars, idx_psfmass, idx_psfimage,
                         psfweights_noclip, psf_clipped, sigma, aperweights_noclip, aper_clipped_all,
                         autoweights_noclip, auto_clipped, with_plots=False, magtype=defaults['magtype']):
    """Calculates WLS fit lines used for determining the goodness of the photometry"""

    # appropriate mag column
    colnames = good_cat_stars.colnames
    if len(colnames) > 6:
        magcol = f'{band}MAG_{magtype}'
        magerrcol = f'{band}MAG_{magtype}'
    else:
        magcol = colnames[2]
        magerrcol = colnames[3]

    # sigma residual fit for auto aperture photometry
    x_auto = good_cat_stars['%s' % magcol][idx_psfmass][~auto_clipped.mask]
    y_auto = cleanPSFsources['%sMAG_AUTO' % band][idx_psfimage][~auto_clipped.mask]
    x_const_auto = sm.add_constant(x_auto)
    model_auto = sm.WLS(y_auto, x_const_auto, weights=autoweights_noclip[~auto_clipped.mask]).fit()

    x_auto_nc = good_cat_stars['%s' % magcol][idx_psfmass]
    y_auto_nc = cleanPSFsources['%sMAG_AUTO' % band][idx_psfimage]
    x_const_auto_nc = sm.add_constant(x_auto_nc)
    model2 = sm.WLS(y_auto_nc, x_const_auto_nc, weights=autoweights_noclip).fit()

    m_auto = model_auto.params[1]
    m_autoerr = model_auto.bse[1]
    b_auto = model_auto.params[0]
    b_autoerr = model_auto.bse[0]
    model_auto_resid = model_auto.resid

    print('Num of crossmatched sources used in auto aperture %s sig fit: %i'
          % (sigma, len(cleanPSFsources['%sMAG_AUTO' % band][idx_psfimage][~auto_clipped.mask])))
    print(' %s sig fit: slope = %.4f +/- %.4f' % (sigma, m_auto, m_autoerr))
    print(' %s sig fit: y-int = %.4f +/- %.4f' % (sigma, b_auto, b_autoerr))

    aper_model_sigs = model_sig = model_sig_resid = 0
    aperweights_noclip = []

    if len(psfweights_noclip) > 0:
        model2 = 0
        x = cleanPSFsources['%sMAG_PSF' % band][idx_psfimage]
        y = good_cat_stars['%s' % magcol][idx_psfmass]
        x_const = sm.add_constant(x)
        # model = sm.WLS(y, x, weights=psfweights).fit()
        model2 = sm.WLS(y, x_const, weights=psfweights_noclip).fit()

        avg2 = np.average(model2.resid, weights=psfweights_noclip)
        var2 = np.average((model2.resid - avg2) ** 2, weights=psfweights_noclip)

        # residual fit - 3 sigma clip
        x_sig = good_cat_stars['%s' % magcol][idx_psfmass][~psf_clipped.mask]
        y_sig = cleanPSFsources['%sMAG_PSF' % band][idx_psfimage][~psf_clipped.mask]
        x_const_sig = sm.add_constant(x_sig)
        model_sig = sm.WLS(y_sig, x_const_sig, weights=psfweights_noclip[~psf_clipped.mask]).fit()
        m_sig = model_sig.params[1]
        m_sigerr = model_sig.bse[1]
        b_sig = model_sig.params[0]
        b_sigerr = model_sig.bse[0]
        model_sig_resid = model_sig.resid

        print('Num of crossmatched sources used in PSF %s sig fit: %i'
              % (sigma, len(cleanPSFsources['%sMAG_PSF' % band][idx_psfimage][~psf_clipped.mask])))
        print(' %s sig fit: slope = %.4f +/- %.4f' % (sigma, m_sig, m_sigerr))
        print(' %s sig fit: y-int = %.4f +/- %.4f' % (sigma, b_sig, b_sigerr))

    # sigma clipped resids for aper mags, disabled right now
    if len(aperweights_noclip) > 0:
        aperweights_noclip = aperweights_noclip.tolist()
        aper_model_sigs= []
        for idx, (aperweight_nc, aperclipped) in enumerate(zip(aperweights_noclip, aper_clipped_all)):
            x_sig_ap = good_cat_stars['%s' % magcol][idx_psfmass][~aperclipped.mask]
            y_sig_ap = cleanPSFsources['%sMAG_APER' % band][:,idx][idx_psfimage][~aperclipped.mask]
            x_const_sig_ap = sm.add_constant(x_sig_ap)
            model_sig_ap = sm.WLS(y_sig_ap, x_const_sig_ap, weights=aperweight_nc[~aperclipped.mask]).fit()
            aper_model_sigs.append(model_sig_ap)

    if with_plots:
        return model2, model_sig, model_sig_resid, model_auto, model_auto_resid, aper_model_sigs
    else:
        if magtype == 'AUTO':
            return m_auto, b_auto, round(3 * b_autoerr, 4)
        elif magtype == 'PSF':
            b_sig_err = round(3 * b_sigerr, 4)
            if b_sig_err > 0.1:
                print(' PSF intercept err > 0.1, reducing acc. err value from 3 sigma to 1.5')
                return m_sig, b_sig, round(1.5 * b_sigerr, 4)
            else:
                return m_sig, b_sig, round(3 * b_sigerr, 4)