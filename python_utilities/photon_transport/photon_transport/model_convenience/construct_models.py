"""
Simplified interfaces for constructing two- or three-layered models with given
parameters, with preloaded absorption and scattering properties.
"""

import photon_transport
import photon_transport.inverse_model.simple_inverse_model as inverse
import numpy as np
import pandas as pd

wlens = np.arange(400, 800)
absorption = photon_transport.forward_model.get_absorption_spectra(wlens)
scattering = photon_transport.forward_model.get_scattering_spectra(wlens)

muam694=300
oxys=[0.3, 0.99]
bvfs=[0.03, 0.02]
thicknesses=[100e-06, 300e-06]

def get_twolayer_model(threelayer_model):
    """
    Fit a two-layer model to a three-layer model, using same epidermal
    properties.
    """

    threelayer_refl = photon_transport.forward_model.get_reflectance(threelayer_model)
    two_layer_model = [dict(threelayer_model[0]), dict(threelayer_model[1])]
    two_layer_model[1]['mua'] = inverse.estimate_muad(two_layer_model, threelayer_refl).copy()
    two_layer_model[1]['d'] = 1.0
    return two_layer_model

def get_threelayer_model(muam694=muam694, oxys=oxys, bvfs=bvfs, thicknesses=thicknesses, return_full=False):
    """
    Get a three-layer model with the given parameters.
    """

    model_params = pd.Series({
        'muam694': muam694,
        'muabo_1': bvfs[0]*oxys[0],
        'muabd_1': bvfs[0]*(1-oxys[0]),
        'muabo_2': bvfs[1]*oxys[1],
        'muabd_2': bvfs[1]*(1-oxys[1]),
        'd_1': thicknesses[0],
        'd_2': thicknesses[0] + thicknesses[1]
        })

    model = photon_transport.forward_model.normal_skin_threelayers(absorption, scattering, muam694=muam694, oxys=oxys, bvfs=bvfs, thicknesses=thicknesses)
    reflectance = photon_transport.forward_model.get_reflectance(model)

    if return_full:
        return model_params, model, reflectance
    else:
        return model
