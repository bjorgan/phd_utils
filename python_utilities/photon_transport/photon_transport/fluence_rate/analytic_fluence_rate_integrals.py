"""
Functions for calculating the analytic integral over fluence rates and
source terms.
"""

from . import fluence_rate_threelayers
from . import fluence_rate
import numpy as np
import photon_transport
import matplotlib.pyplot as plt
import photon_transport.inverse_model.simple_inverse_model as inverse

def analytic_integral_twolayer(model):
    """
    Calculate integral over the fluence rate in a two-layered model.

    Parameters
    ----------
    model: skin model structure (see forward_model.normal_skin())
        Skin model, assumed to be two-layered.
    Returns
    -------
    layer_1:
        Integral over layer 1
    layer_2:
        Integral over layer 2
    """
    layer_1, layer_2 = analytic_depthresolved_integral_twolayer(model)
    return layer_1, layer_2

def analytic_integral_threelayer(model):
    """
    Calculate integral over the fluence rate in a three-layered model.

    Parameters
    ----------
    model: skin model structure (see forward_model.normal_skin())
        Skin model, assumed to be three-layered.
    Returns
    -------
    layer_1:
        Integral over layer 1.
    layer_2:
        Integral over layer 2.
    layer_3:
        Integral over layer 3.
    """

    mua_1, mua_2, mua_3, musr_1, musr_2, musr_3, d_1, d_2 = photon_transport.forward_model.skin_model_to_optical_parameters(model)

    mutr_1 = musr_1 + mua_1
    mutr_2 = musr_2 + mua_2
    mutr_3 = musr_3 + mua_3

    D_1 = 1.0/(3*mutr_1)
    D_2 = 1.0/(3*mutr_2)
    D_3 = 1.0/(3*mutr_3)

    del_1 = np.sqrt(D_1/mua_1)
    del_2 = np.sqrt(D_2/mua_2)
    del_3 = np.sqrt(D_3/mua_3)

    #source terms
    S_1 = (del_1**2)*musr_1/(D_1*(1-(mutr_1**2)*(del_1**2)))
    S_2 = (del_2**2)*musr_2/(D_2*(1-(mutr_2**2)*(del_2**2)))*np.exp(-mutr_1*d_1)*np.exp(mutr_2*d_1)
    S_3 = (del_3**2)*musr_3/(D_3*(1-(mutr_3**2)*(del_3**2)))*np.exp(-mutr_1*d_1)*np.exp(-mutr_2*(d_2-d_1))*np.exp(mutr_3*d_2)

    #coefficients
    A_1, A_2, A_3, A_4, A_5 = fluence_rate_threelayers.analytic_coefficients(mua_1, mua_2, mua_3, musr_1, musr_2, musr_3, d_1, d_2)

    #integral terms
    layer_1_contrib = lambda x: -S_1/mutr_1 * np.exp(-mutr_1*x) - A_1*del_1*np.exp(-x/del_1) + A_2*del_1*np.exp(x/del_1)
    layer_2_contrib = lambda x: -S_2/mutr_2 * np.exp(-mutr_2*x) - A_3*del_2*np.exp(-x/del_2) + A_4*del_2*np.exp(x/del_2)
    layer_3_contrib = lambda x: -S_3/mutr_3 * np.exp(-mutr_3*x) - A_5*del_3*np.exp(-x/del_3)

    layer_1_integral = layer_1_contrib(d_1) - layer_1_contrib(0)
    layer_2_integral = layer_2_contrib(d_2) - layer_2_contrib(d_1)
    layer_3_integral = -layer_3_contrib(d_2)

    return layer_1_integral, layer_2_integral, layer_3_integral

def analytic_depthresolved_integral_twolayer(model, depth = None):
    """
    Integral over the two-layered fluence rate, with possibility for integrating
    over arbitrary depths.

    Parameters
    ----------
    model: skin model
        Two-layered skin model.
    depth: float, optional
        If set to None, will integrate over the first and second layer and
        return them separatedly. If set to a specific depth, will integrate
        from the start depth of the second layer and to the specified depth.
    Returns
    -------
    (layer_1_full_integral, layer_2_full_integral): tuple
        Integrals over each separate layer, if depth was set to None.
    layer_2_variable_integral: float
        Integral to the specified depth, if depth is not none.
    """

    mua_1, mua_2, musr_1, musr_2, d_1 = photon_transport.forward_model.skin_model_to_optical_parameters(model)

    A_1, A_2, A_3 = fluence_rate.analytic_coefficients(mua_1, mua_2, musr_1, musr_2, d_1)

    mutr_1 = musr_1 + mua_1
    mutr_2 = musr_2 + mua_2

    D_1 = 1.0/(3*mutr_1)
    D_2 = 1.0/(3*mutr_2)

    del_1 = np.sqrt(D_1/mua_1)
    del_2 = np.sqrt(D_2/mua_2)

    S_1 = (del_1**2)*musr_1/(D_1*(1-(mutr_1**2)*(del_1**2)))
    S_2 = (del_2**2)*musr_2/(D_2*(1-(mutr_2**2)*(del_2**2)))*np.exp(-mutr_1*d_1)*np.exp(mutr_2*d_1)

    layer_1_contrib = lambda x: -S_1/mutr_1 * np.exp(-mutr_1*x) - A_1*del_1*np.exp(-x/del_1) + A_2*del_1*np.exp(x/del_1)
    layer_2_contrib = lambda x: -S_2/mutr_2 * np.exp(-mutr_2*x) - A_3*del_2*np.exp(-x/del_2)

    layer_1_full_integral = layer_1_contrib(d_1) - layer_1_contrib(0)

    layer_2_full_integral = -layer_2_contrib(d_1)

    if depth is None:
        return layer_1_full_integral, layer_2_full_integral
    else:
        layer_2_variable_integral = layer_2_contrib(depth)-layer_2_contrib(d_1)
        return layer_2_variable_integral

def analytic_dermis_integral_threelayer(model, return_separate=False):
    """
    (Legacy)
    Does the same as analytic_integral_threelayer, but kept for compatibility
    reasons as some scripts assume the existence of this function, pluss needs
    the option for return_separate=False.

    Parameters
    ----------
    model: skin model
        Three-layer model
    return_separate: optional, boolean
        Whether to return the separate layer terms as in
        analytic_integral_threelayer, or whether to multiply by absorption
        coefficients and return them along with the sum over the integrals
    Returns
    -------
    Fluence rate integrals, dependent on the return_separate flag.
    """

    mua_1, mua_2, mua_3, musr_1, musr_2, musr_3, d_1, d_2 = photon_transport.forward_model.skin_model_to_optical_parameters(model)
    layer_1, layer_2, layer_3 = analytic_integral_threelayer(model)

    if return_separate:
        return layer_1_integral, layer_2_integral, layer_3_integral
    else:
        return mua_1*layer_1_integral + mua_2*layer_2_integral + mua_3*layer_3_integral, layer_1_integral + layer_2_integral + layer_3_integral

def analytic_dermis_integral_twolayer(model, return_separate=False):
    """
    See analytic_dermis_integral_twolayer().
    """
    mua_1, mua_2, musr_1, musr_2, d_1 = photon_transport.forward_model.skin_model_to_optical_parameters(model)

    layer_1_integral, layer_2_integral = analytic_integral_twolayer(model)

    if return_separate:
        return layer_1_integral, layer_2_integral
    else:
        return layer_1_integral*mua_1 + layer_2_integral*mua_2, layer_1_integral + layer_2_integral

def source_term_integral_full(model):
    """
    Integral over the isotropic source term in the diffusion equation, from 0 to infinity.
    """
    dermal_contribution = source_term_integral(model)

    #estimate contribution from first layer
    mutr_1 = model[0]['mua'] + model[0]['musr']
    musr_1 = model[0]['musr']
    epidermal_contribution = musr_1/mutr_1*(1-np.exp(-mutr_1*model[0]['d']))
    return epidermal_contribution + dermal_contribution

def source_term_integral(model):
    """
    Integral over the isotropic source term in the diffusion equation, from d_1 to infinity.
    """

    mutr_1 = model[0]['mua'] + model[0]['musr']
    if len(model) == 2:
        mutr_2_2 = model[1]['mua'] + model[1]['musr']
        S_integral = model[1]['musr']/mutr_2_2*np.exp(-mutr_1*model[0]['d'])

    if len(model) == 3:
        musr_2_3 = model[1]['musr']
        musr_3_3 = model[2]['musr']
        mutr_2_3 = model[1]['mua'] + musr_2_3
        mutr_3_3 = model[2]['mua'] + musr_3_3
        d_1 = model[0]['d']
        d_2 = model[1]['d']
        S_integral = np.exp(-mutr_1*d_1)*(musr_2_3/mutr_2_3 + np.exp(-mutr_2_3*(d_2 - d_1))*(musr_3_3/mutr_3_3 - musr_2_3/mutr_2_3))

    return S_integral
