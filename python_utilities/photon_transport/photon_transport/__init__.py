"""
Routines and utilities related to diffusion-model based photon transport
modeling of human tissue.

List over submodules:
    * fluence_rate/: Routines for calculating fluence rate, irradiance and
    radiance in various model configurations, and the integrals over such
    properties.
    * inverse_model/: Routines related to applying the model in inverse, either
    convenience methods or full inverse methods.
    * optimization/: Also contains routines for inverse models, but with a more
    convenient interface, onelayer and twolayer. Also has routines for
    crossvalidation.
    * model_convenience/: Convenience methods for constructing forward models.
    * opt_params/: Data for optical properties of various materials.
    * forward_model.py: Methods for constructing a model of human skin and
    predicting the diffuse reflectance.
"""

from . import forward_model
from . import diffusion_model_solutions
