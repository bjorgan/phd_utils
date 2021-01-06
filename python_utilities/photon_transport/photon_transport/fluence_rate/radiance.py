"""
Solutions for radiance.

Uses the solution in fluence_rate to calculate fluence rate and photon flux, and is thus subject to the same constraints and numerical error as those models.

Radiance is dependent on the direction we are facing. We input only the z-component of the direction vector, since we model only the z-component of the direction vector j.
"""

import numpy as np
import fluence_rate
import photon_transport

def radiance(phi, j, direction_vector_z_component):
    """
    Calculate the radiance L from the irradiance phi and photon flux j.
    """
    return 1.0/(4.0*np.pi)*phi + 3.0/(4.0*np.pi)*j*direction_vector_z_component

def two_layers_isotropic(x, mua_1, mua_2, musr_1, musr_2, direction_vector_z_component, d_1 = 100e-06):
    """
    Calculate the radiance from a two-layered model with isotropic source function.
    """
    phi, j = fluence_rate.two_layers_isotropic(x, mua_1, mua_2, musr_1, musr_2, d_1)
    return radiance(phi, j, direction_vector_z_component)

def two_layers_isotropic_from_skin_model(x, skin_model, direction_vector_z_component, wavelength_index=0):
    """
    Calculate the radiance from a two-layered model with isotropic source
    functions, using a skin model structure (as it would be obtained from
    normal_skin() in forward_model).
    """
    mua_1, mua_2, musr_1, musr_2, d_1 = photon_transport.forward_model.skin_model_to_optical_parameters(skin_model)
    return two_layers_isotropic(x, mua_1[wavelength_index], mua_2[wavelength_index], musr_1[wavelength_index], musr_2[wavelength_index], direction_vector_z_component, d_1)
