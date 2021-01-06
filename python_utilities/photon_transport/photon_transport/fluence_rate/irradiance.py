"""
Solutions for irradiance.

Uses the solution in fluence_rate to calculate fluence rate and photon flux,
and is thus subject to the same constraints and numerical error as those
models.

Irradiance is radiance integrated over the half-sphere, and is thus dependent
on the direction we are facing.
"""

import numpy as np
from . import fluence_rate
import photon_transport

def irradiance(phi, j, direction):
    """
    Calculate the irradiance given fluence rate phi
    and photon flux j.
    """
    return phi/4 -1*direction*j/2

def two_layers_isotropic(x, mua_1, mua_2, musr_1, musr_2, direction, d_1 = 100e-06):
    """
    Calculate the irradiance for a two-layered model with isotropic source
    functions.
    """
    phi, j = fluence_rate.two_layers_isotropic(x, mua_1, mua_2, musr_1, musr_2, d_1)
    return irradiance(phi, j, direction)

def two_layers_isotropic_from_skin_model(x, skin_model, direction, wavelength_index=None):
    """
    Calculate the irradiance for a two-layered skin model with isotropic source
    functions, with skin_model as input.

    Parameters
    ----------
    x: float or array
        Depth in tissue
    skin_model: array
        Skin model as obtained from forward_model.normal_skin()
    direction: integer
        Direction in tissue, 1 corresponding to forward direction (positive
        z-direction), -1 corresponding to irradiance travelling back
    wavelength_index: integer, optional
        Wavelength index, will calculate for all wavelength indices present in
        skin_model if set to None (default)
    Returns
    -------
    irradiance: float or array
        Irradiance, as float value or array over floats, depending on input type.
    """

    mua_1, mua_2, musr_1, musr_2, d_1 = photon_transport.forward_model.skin_model_to_optical_parameters(skin_model)

    if wavelength_index is not None:
        return two_layers_isotropic(x, mua_1[wavelength_index], mua_2[wavelength_index], musr_1[wavelength_index], musr_2[wavelength_index], direction, d_1)
    else:
        return two_layers_isotropic(x, mua_1, mua_2, musr_1, musr_2, direction, d_1)

def infinite_isotropic_medium(x, mua, musr):
    """
    Calculate the irradiance for an infinite medium with isotropic source
    functions.
    """
    phi, j = fluence_rate.infinite_isotropic_medium(x, mua, musr)
    return irradiance(phi, j)
