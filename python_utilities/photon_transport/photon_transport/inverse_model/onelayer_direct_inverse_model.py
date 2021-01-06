"""
Calculate musr/mua directly from a reflectance spectrum using
ReflIsoL1 in reverse, analytically.
"""

import numpy as np

from photon_transport.forward_model import A_default

def get_musr_divided_by_mua(refl, A=A_default):
    """
    Assuming a one-layered model with isotropic source functions (ReflIsoL1),
    calculate musr/mua from the given reflectance spectrum.

    Parameters
    ----------
    refl:
        Reflectance
    Returns
    -------
    (sol_1, sol_2): tuple
        Solutions to the equation ReflIsoL1(musr/mua) = refl
    """

    a = A**2*(1-refl)**2
    b = -2*A*(1-refl)*refl*(A+1) - refl**2*((1+3*A)/np.sqrt(3))**2
    c = refl**2*(A+1)**2 - refl**2*((1 + 3*A)/(np.sqrt(3)))**2

    discriminant = np.sqrt(b**2 - 4*a*c)

    solution_1 = (-b + discriminant)/(2*a)
    solution_2 = (-b - discriminant)/(2*a)
    return solution_1, solution_2

