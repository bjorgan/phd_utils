"""
Simple ad-hoc inverse model, based on https://github.com/ntnu-bioopt/gpudm.

See fit_reflectance(...) for entry method for doing the full chain, estimating both
muam694 and skin parameters in dermis.
"""

from photon_transport import forward_model
from . import derivatives
import scipy.optimize
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def estimate_melanin(wavelengths, reflectance):
    """
    Estimate melanin content (as melanin absorption at 694 nm, muam694) using
    the method in https://dx.doi.org/10.1117/1.JBO.19.6.066003.

    General philosophy is to set a low melanin content value, try to fit
    chromophores including melanin in dermis, remove melanin from the dermis
    and then try to fit the melanin in epidermis given the dermal properties,
    in two iterations.

    Parameters
    ----------
    wavelengths: ndarray
        Wavelength array
    reflectance: ndarray
        Reflectance
    Returns
    -------
    melanin: float
        Melanin content
    """
    start_wlen = 730
    end_wlen = 820
    num_iterations = 2

    subset = (wavelengths >= start_wlen) & (wavelengths <= end_wlen)
    wavelengths = wavelengths[subset]
    reflectance = reflectance[subset]

    scattering = forward_model.get_scattering_spectra(wavelengths)
    absorption = forward_model.get_absorption_spectra(wavelengths)

    muam694 = 100
    oxy = 0.5
    bvf = 0.01

    for i in range(0, num_iterations):
        skin_model = forward_model.normal_skin(absorption, scattering, muam694=muam694, oxy=oxy, bvf=bvf)

        #get dermal absorption coefficient given current muam694
        muad = estimate_muad(skin_model, reflectance)

        #fit properties
        result, fitted_muad, used_spectra = unmix_absorption(absorption[['wavelength', 'melanin', 'muabo', 'muabd']], muad)

        #remove melanin from dermis
        fitted_muad.absorption -= result.melanin * used_spectra.melanin.values
        skin_model[1]['mua'] = fitted_muad.absorption

        #get epidermal absorption given current dermal absorption coefficient minus melanin
        muae = estimate_muae(skin_model, reflectance)
        result, fitted_muae, used_spectra = unmix_absorption(absorption[['wavelength', 'melanin']], muae)
        muam694 = result.melanin
    return muam694

def reconstruct_reflectance(muam694, fitted_dermal_absorption):
    """
    Reconstruct diffuse reflectance from fitted absorption coefficient as
    obtained from unmix_absorption(...).

    Parameters
    ----------
    muam694: float
        Melanin content in epidermis.
    fitted_dermal_absorption: ndarray
        Fitted dermal absorption coefficient.
    Returns
    -------
    reflectance: ndarray
        Reflectance
    """
    skin_model = forward_model.normal_skin(forward_model.get_absorption_spectra(fitted_dermal_absorption.wavelength),
                                           forward_model.get_scattering_spectra(fitted_dermal_absorption.wavelength),
                                           muam694=muam694,
                                           oxy=0.6,
                                           bvf=0.01)
    skin_model[1]['mua'] = fitted_dermal_absorption.absorption
    return forward_model.get_reflectance(skin_model)

def params_to_oxy_bvf(params):
    """
    Convert parameters as obtained from unmix_absorption(...) to oxy, bvf

    Parameters
    ----------
    params: pd.DataFrame
        Pandas dataframe/series containing columns 'muabo' and 'muabd'
    Returns
    -------
    bvf_oxy: pd.Series
        Pandas series containing 'bvf', 'oxy'
    """
    bvf = params['muabo'] + params['muabd']
    oxy = params['muabo']/bvf
    return pd.Series({'oxy': oxy, 'bvf': bvf})

def fit_results_to_physical_parameters(fit_results):
    """
    Convert fitted skin parameters as obtained from fit_reflectance(...) to oxy, bvf.

    Parameters
    ----------
    fit_result:
        Fitted results as series or dataframe
    Returns
    -------
    res:
        Series or Dataframe containing the physical properties determined from the fitting parameters
    """
    return pd.concat([params_to_oxy_bvf(fit_results.params_lower), params_to_oxy_bvf(fit_results.params_upper)], keys=['params_lower', 'params_upper'], axis=1)

def fit_reflectance(wavelengths, reflectance, muam694=None):
    """
    Entry method for fitting melanin content and skin parameters in dermis in
    two wavelength intervals, using the general method from
    https://dx.doi.org/10.1117/1.JBO.19.6.066003.

    General philosophy is to fit melanin using estimate_melanin(...), and then
    fit the necessary absorption coefficient in the lower layer to match
    measured reflectance and the melanin content. This absorption coefficient
    is then fit using various absorption spectra for shorter wavelengths (low
    penetration depth) and for longer wavelengths (high penetration depth)
    using unmix_absorption(...).

    Parameters
    ----------
    wavelength:
        Wavelengths
    reflectance:
        Measured reflectance
    muam694:
        Melanin content, estimated if unspecified
    Returns
    -------
    skin_param:
        Skin parameters and fitted results as a pandas series.
    """
    if muam694 is None:
        #estimate melanin content
        muam694 = estimate_melanin(wavelengths, reflectance)

    #estimate dermal absorption coefficient
    scattering = forward_model.get_scattering_spectra(wavelengths)
    absorption = forward_model.get_absorption_spectra(wavelengths)
    skin_model = forward_model.normal_skin(absorption, scattering, muam694=muam694, oxy=0.6, bvf=0.01)
    muad = estimate_muad(skin_model, reflectance)
    muad_fit = pd.Series({'wavelength': wavelengths, 'absorption': muad})

    #fit dermal absorption coefficient in two wavelength intervals
    params_500, fitted_500, absorption_500 = unmix_absorption(absorption[['wavelength', 'constant', 'melanin', 'muabd', 'muabo']], muad, (510, 590))
    params_700, fitted_700, absorption_700 = unmix_absorption(absorption[['wavelength', 'constant', 'melanin', 'muabd', 'muabo', 'segelstein81']], muad, (690, 820))

    #reconstruct reflectance from fitted absorptions
    fitted_500['reflectance'] = reconstruct_reflectance(muam694, fitted_500)
    fitted_700['reflectance'] = reconstruct_reflectance(muam694, fitted_700)

    #prepare return values
    ret_series = pd.concat([params_500,
                      params_700,
                      fitted_500,
                      fitted_700,
                      muad_fit],
                    keys=['params_lower',
                          'params_upper',
                          'fitted_lower',
                          'fitted_upper',
                          'muad'])
    ret_series[('params_both', 'muam694')] = muam694
    return ret_series

"""Number of iterations used in estimation of muad."""
NUM_NEWTON_ITERATIONS = 15

def replace_nan(array):
    if hasattr(array, '__len__'):
        array[np.isnan(array)] = 0
    elif np.isnan(array):
        array = 0
    return array

def _root_function(mua, reflectance, skin_model):
    return forward_model.ReflIsoL2(skin_model[0]['mua'], mua, skin_model[0]['musr'], skin_model[1]['musr'], skin_model[0]['d']) - reflectance

def _jacobian(mua, skin_model):
    _, derivative = derivatives.ReflIsoL2DerivMuad(skin_model[0]['d'], skin_model[0]['mua'], mua, skin_model[0]['musr'], skin_model[1]['musr'])
    jac = np.diag(derivative)
    return jac

def estimate_muad(skin_model, reflectance, num_newton_iterations=None):
    """
    Estimate dermal absorption coefficient from human skin reflectance, given a two-layered skin model and the melanin content in epidermis.

    Parameters
    ----------
    skin_model: dict
        Skin model as obtained from forward_model.normal_skin(...).
    reflectance: ndarray, float
        Reflectance
    num_newton_iterations: optional
        Number of iterations of Newton's method. If set to None, iterates until
        convergence. NOTE: This was previously set to 15, as determined from
        the master thesis, but turns out that the number of iterations is
        dependent on the scattering parameters... See
        improve_inverse.../parameter_and_crossvalidation_investigation/check_newton_fitting_errors.py.
        If set to None, will estimate using scipy.optimize.root instead.

    Returns
    -------
    mua2: ndarray, float
        Absorption coefficient in dermis for each reflectance value
    """

    mua2 = skin_model[1]['mua']

    if num_newton_iterations is not None:
        for i in range(0, num_newton_iterations):
            simulated_reflectance, derivative = derivatives.ReflIsoL2DerivMuad(skin_model[0]['d'], skin_model[0]['mua'], mua2, skin_model[0]['musr'], skin_model[1]['musr'])
            mua2 = mua2 - (simulated_reflectance - reflectance)/derivative

    else:
        mua2 = scipy.optimize.root(lambda mua: _root_function(mua, reflectance, skin_model), mua2, jac=lambda mua: _jacobian(mua, skin_model)).x

    return mua2

def estimate_muae(skin_model, reflectance):
    """
    Estimate epidermal absorption coefficient from human skin reflectance,
    given a two-layered skin model and fixed properties in dermis.

    Parameters
    ----------
    skin_model:
        Skin model as obtained from forward_model.normal_skin(...).
    Returns
    -------
    reflectance:
        Reflectance
    """
    mua1 = skin_model[0]['mua']
    for i in range(0, NUM_NEWTON_ITERATIONS):
        simulated_reflectance, derivative = derivatives.ReflIsoL2DerivMuae(skin_model[0]['d'], mua1, skin_model[1]['mua'], skin_model[0]['musr'], skin_model[1]['musr'])
        mua1 = mua1 - (simulated_reflectance - reflectance)/derivative
    return mua1

def unmix_absorption(absorption_spectra, absorption, wlen_interval = (None, None)):
    """
    Unmix an absorption spectrum using absorption specta from known
    constitutents.  Applies standard non-negative least squares.

    Parameters
    ----------
    absorption_spectra:
        Absorption spectrum matrix, as obtained from forward_model.get_absorption_spectra(...) (and possibly reduced to a smaller subset of chromophores)
    absorption:
        Absorption coefficients to unmix
    start_wlen:
        Start wavelength, set to the first wavelength if unspecified
    end_wlen:
        End wavelength, set to last wavelength if unspecified
    Returns
    -------
    res:
        Tuple containing fitted parameters, reconstructed absorption coefficient and the absorption matrix reduced to the specified wavelength range
    """
    #replace nan values by neighbouring values
    nan_values = np.isnan(absorption)
    absorption[nan_values] = np.interp(np.flatnonzero(nan_values), np.flatnonzero(~nan_values), absorption[~nan_values])

    start_wlen = wlen_interval[0]
    end_wlen = wlen_interval[1]

    #prepare data matrices using start and end wavelengths
    if start_wlen is None:
        start_wlen = absorption_spectra.wavelength.iloc[0]
    if end_wlen is None:
        end_wlen = absorption_spectra.wavelength.iloc[-1]

    subset = (absorption_spectra.wavelength >= start_wlen) & (absorption_spectra.wavelength <= end_wlen)
    wavelength = absorption_spectra.wavelength[subset]
    absorption_spectra = absorption_spectra[subset].drop('wavelength', axis=1)
    absorption = absorption[subset]

    #fit using non-negative least squares
    result, residual = scipy.optimize.nnls(absorption_spectra, absorption)
    fitted_absorption = pd.Series({'wavelength': wavelength.values, 'absorption': np.dot(absorption_spectra, result)})
    result = pd.Series(data=result, index=absorption_spectra.columns)
    return result, fitted_absorption, absorption_spectra
