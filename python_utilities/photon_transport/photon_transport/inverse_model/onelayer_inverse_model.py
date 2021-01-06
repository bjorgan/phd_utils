"""
Inverse model based on only one layer. Assumes a scattering spectrum, and
estimates the absorption spectrum.
"""

from photon_transport import forward_model
import scipy.optimize
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def fit_reflectance(wavelengths, reflectance, musrr500=forward_model.musrr500_bashkatov, musmr500=forward_model.musmr500_bashkatov, return_musr=False):
    """
    Entry method. Fits ReflIsoL1 to the input reflectance by assuming a
    scattering spectrum (default scattering from the forward model module).

    Parameters
    ----------
    wavelength:
        Wavelengths
    reflectance
        Measured reflectance
    Returns
    -------
    muad:
        Absorption spectrum
    musr: optional
        Scattering spectrum, if return_musr is set to True
    """

    #construct scattering spectrum
    scattering_params = forward_model.get_scattering_spectra(wavelengths)
    musr = (100*(scattering_params.rayleigh_reduced * musrr500 + scattering_params.mie_reduced*musmr500)).as_matrix()

    #estimate absorption coefficient
    muad = estimate_mua(musr, reflectance)

    if return_musr:
        return muad, musr
    else:
        return muad

def estimate_mua(musr, reflectance):
    """
    Estimate absorption coefficients from the input reflectance assuming an
    input scattering spectrum.
    """
    mua_ret = np.ones(len(musr))

    i=0
    for mua in mua_ret:
        initial_value = mua
        mua_ret[i] = scipy.optimize.newton(
                     lambda mua_input:
                     (forward_model.ReflIsoL1(mua_input, musr[i]) - reflectance[i]),
                     initial_value)
        i += 1
    return mua_ret
