"""
Solution for fluence rate in three-layered medium.

The solution will probably fail in horrible ways for some singularities due to
the source terms.
"""

import numpy as np
import photon_transport


def three_layers_isotropic_from_skin_model(x, skin_model, wavelength_index=None):
    """
    Calculate the fluence rate for a skin_model structure.
    """

    mua_1, mua_2, mua_3, musr_1, musr_2, musr_3, d_1, d_2 = photon_transport.forward_model.skin_model_to_optical_parameters(skin_model)
    if wavelength_index is None:
        return three_layers_isotropic(x, mua_1, mua_2, mua_3, musr_1, musr_2, musr_3, d_1, d_2)
    else:
        return three_layers_isotropic(x, mua_1[wavelength_index], mua_2[wavelength_index], mua_3[wavelength_index], musr_1[wavelength_index], musr_2[wavelength_index], musr_3[wavelength_index], d_1, d_2)

def three_layers_isotropic(x, mua_1, mua_2, mua_3, musr_1, musr_2, musr_3, d_1, d_2):
    """
    Fluence rate in a three-layered medium, assuming isotropic source
    functions.  From L. O. Svaasand et al., "Tissue parameters determining the
    visual appearance of normal skin and port-wine stains", Laser Med. Sci., 10
    (1995)
    """

    A = photon_transport.forward_model.A

    #transmitted light density after specular reflection
    P0 = 1

    mutr_1 = musr_1 + mua_1
    mutr_2 = musr_2 + mua_2
    mutr_3 = musr_3 + mua_3

    D_1 = 1.0/(3*mutr_1)
    D_2 = 1.0/(3*mutr_2)
    D_3 = 1.0/(3*mutr_3)

    del_1 = np.sqrt(D_1/mua_1)
    del_2 = np.sqrt(D_2/mua_2)
    del_3 = np.sqrt(D_3/mua_3)

    #solution terms in layer 1
    phi_source_constant_layer_1 = P0*(del_1**2)*musr_1/(D_1*(1-(mutr_1**2)*(del_1**2)))
    phi_source_layer_1 = lambda x: phi_source_constant_layer_1*np.exp(-mutr_1*x)
    phi_term_1_layer_1 = lambda x: np.exp(-x/del_1)
    phi_term_2_layer_1 = lambda x: np.exp(x/del_1)

    j_source_layer_1 = lambda x: -D_1*(-mutr_1*phi_source_layer_1(x))
    j_term_1_layer_1 = lambda x: -D_1*(-1.0/del_1*phi_term_1_layer_1(x))
    j_term_2_layer_1 = lambda x: -D_1*(1.0/del_1*phi_term_2_layer_1(x))

    #solution terms in layer 2
    phi_source_constant_layer_2 = P0*(del_2**2)*musr_2/(D_2*(1-(mutr_2**2)*(del_2**2)))
    phi_source_layer_2 = lambda x: phi_source_constant_layer_2*np.exp(-mutr_1*d_1) * np.exp(-mutr_2*(x-d_1))
    phi_term_1_layer_2 = lambda x: np.exp(-x/del_2)
    phi_term_2_layer_2 = lambda x: np.exp(x/del_2)

    j_source_layer_2 = lambda x: -D_2*(-mutr_2*phi_source_layer_2(x))
    j_term_1_layer_2 = lambda x: -D_2*(-1.0/del_2*phi_term_1_layer_2(x))
    j_term_2_layer_2 = lambda x: -D_2*(1.0/del_2*phi_term_2_layer_2(x))

    #solution terms in layer 3
    phi_source_constant_layer_3 = P0*(del_3**2)*musr_3/(D_3*(1-(mutr_3**2)*(del_3**2)))
    phi_source_layer_3 = lambda x: phi_source_constant_layer_3*np.exp(-mutr_1*d_1) * np.exp(-mutr_2*(d_2-d_1)) * np.exp(-mutr_3*(x - d_2))
    phi_term_layer_3 = lambda x: np.exp(-x/del_3)

    j_source_layer_3 = lambda x: -D_3*(-mutr_3*phi_source_layer_3(x))
    j_term_layer_3 = lambda x: -D_3*(-1.0/del_3*phi_term_layer_3(x))

    #get coefficients
    A_1, A_2, A_3, A_4, A_5 = analytic_coefficients(mua_1, mua_2, mua_3, musr_1, musr_2, musr_3, d_1, d_2)

    #full solution
    phi_1 = lambda x: phi_source_layer_1(x) + A_1*phi_term_1_layer_1(x) + A_2*phi_term_2_layer_1(x)
    phi_2 = lambda x: phi_source_layer_2(x) + A_3*phi_term_1_layer_2(x) + A_4*phi_term_2_layer_2(x)
    phi_3 = lambda x: phi_source_layer_3(x) + A_5*phi_term_layer_3(x)

    j_1 = lambda x: j_source_layer_1(x) + A_1*j_term_1_layer_1(x) + A_2*j_term_2_layer_1(x)
    j_2 = lambda x: j_source_layer_2(x) + A_3*j_term_1_layer_2(x) + A_4*j_term_2_layer_2(x)
    j_3 = lambda x: j_source_layer_3(x) + A_5*j_term_layer_3(x)

    if not hasattr(x, '__len__'):
        x = np.array([x])

    layer_1 = x < d_1
    layer_2 = (x >= d_1) & (x < d_2)
    layer_3 = x > d_2

    x_1 = x[layer_1]
    x_2 = x[layer_2]
    x_3 = x[layer_3]

    phi = np.zeros(len(x))
    phi[layer_1] = phi_1(x_1)
    phi[layer_2] = phi_2(x_2)
    phi[layer_3] = phi_3(x_3)

    j = np.zeros(len(x))
    j[layer_1] = j_1(x_1)
    j[layer_2] = j_2(x_2)
    j[layer_3] = j_3(x_3)

    return phi, j


def analytic_coefficients(mua_1, mua_2, mua_3, musr_1, musr_2, musr_3, d_1, d_2):
    """
    Analytic solutions for the cofficients A1, A2, A3, A4 and A5 in the
    solution for phi given in Svaasand1995.  Used for calculating the fluence
    rate solution above.

    See ReflIsoL3.mw in the penetration_depth folder for derivation.
    """
    A = photon_transport.forward_model.A

    mutr_1 = musr_1 + mua_1
    mutr_2 = musr_2 + mua_2
    mutr_3 = musr_3 + mua_3

    D_1 = 1.0/(3*mutr_1)
    D_2 = 1.0/(3*mutr_2)
    D_3 = 1.0/(3*mutr_3)

    del_1 = np.sqrt(D_1/mua_1)
    del_2 = np.sqrt(D_2/mua_2)
    del_3 = np.sqrt(D_3/mua_3)

    K_1 = del_1*del_1*musr_1/(D_1*(1 - mutr_1**2 * del_1**2))
    K_2 = del_2*del_2*musr_2/(D_2*(1 - mutr_2**2 * del_2**2))
    K_3 = del_3*del_3*musr_3/(D_3*(1 - mutr_3**2 * del_3**2))

    A_1 = ((-2 * np.exp(-mutr_1 * d_1) * D_2 * del_2 * (A * del_1 - D_1) * np.exp(d_1 / del_2) * (del_3 * mutr_2 * D_2 * K_2 - D_3 * (K_3 * del_3 * mutr_3 + K_2 - K_3)) * np.exp(mutr_2 * (-d_2 + d_1)) - np.exp(d_2 / del_2) * (D_2 * del_3 + D_3 * del_2) * (K_1 * (D_1 * del_2 + D_2 * del_1) * (D_1 * mutr_1 + A) * np.exp(d_1 / del_1) + np.exp(-mutr_1 * d_1) * (A * del_1 - D_1) * ((-K_2 * del_2 * mutr_2 - K_1 + K_2) * D_2 + del_2 * mutr_1 * D_1 * K_1))) * np.exp(-d_1 / del_2) + (-K_1 * (-D_1 * del_2 + D_2 * del_1) * (D_1 * mutr_1 + A) * np.exp(d_1 / del_1) + np.exp(-mutr_1 * d_1) * (A * del_1 - D_1) * ((-K_2 * del_2 * mutr_2 + K_1 - K_2) * D_2 + del_2 * mutr_1 * D_1 * K_1)) * np.exp(d_1 / del_2) * np.exp(-d_2 / del_2) * (-D_2 * del_3 + D_3 * del_2)) * del_1 / (np.exp(d_2 / del_2) * (D_2 * del_3 + D_3 * del_2) * (-(-D_1 * del_2 + D_2 * del_1) * (A * del_1 - D_1) * np.exp(-d_1 / del_1) + np.exp(d_1 / del_1) * (D_1 * del_2 + D_2 * del_1) * (A * del_1 + D_1)) * np.exp(-d_1 / del_2) + np.exp(d_1 / del_2) * (-(D_1 * del_2 + D_2 * del_1) * (A * del_1 - D_1) * np.exp(-d_1 / del_1) + np.exp(d_1 / del_1) * (-D_1 * del_2 + D_2 * del_1) * (A * del_1 + D_1)) * np.exp(-d_2 / del_2) * (-D_2 * del_3 + D_3 * del_2))

    A_2 = -((-np.exp(d_2 / del_2) * K_1 * (D_2 * del_3 + D_3 * del_2) * (-D_1 * del_2 + D_2 * del_1) * (D_1 * mutr_1 + A) * np.exp(-d_1 / del_1) - np.exp(-mutr_1 * d_1) * (2 * D_2 * del_2 * np.exp(d_1 / del_2) * (del_3 * mutr_2 * D_2 * K_2 - D_3 * (K_3 * del_3 * mutr_3 + K_2 - K_3)) * np.exp(mutr_2 * (-d_2 + d_1)) + np.exp(d_2 / del_2) * (D_2 * del_3 + D_3 * del_2) * ((-K_2 * del_2 * mutr_2 - K_1 + K_2) * D_2 + del_2 * mutr_1 * D_1 * K_1)) * (A * del_1 + D_1)) * np.exp(-d_1 / del_2) + (-K_1 * (D_1 * del_2 + D_2 * del_1) * (D_1 * mutr_1 + A) * np.exp(-d_1 / del_1) + np.exp(-mutr_1 * d_1) * ((-K_2 * del_2 * mutr_2 + K_1 - K_2) * D_2 + del_2 * mutr_1 * D_1 * K_1) * (A * del_1 + D_1)) * np.exp(d_1 / del_2) * np.exp(-d_2 / del_2) * (-D_2 * del_3 + D_3 * del_2)) * del_1 / (np.exp(d_2 / del_2) * (D_2 * del_3 + D_3 * del_2) * (-(-D_1 * del_2 + D_2 * del_1) * (A * del_1 - D_1) * np.exp(-d_1 / del_1) + np.exp(d_1 / del_1) * (D_1 * del_2 + D_2 * del_1) * (A * del_1 + D_1)) * np.exp(-d_1 / del_2) + np.exp(d_1 / del_2) * (-(D_1 * del_2 + D_2 * del_1) * (A * del_1 - D_1) * np.exp(-d_1 / del_1) + np.exp(d_1 / del_1) * (-D_1 * del_2 + D_2 * del_1) * (A * del_1 + D_1)) * np.exp(-d_2 / del_2) * (-D_2 * del_3 + D_3 * del_2))

    A_3 = del_2 * ((-np.exp(-mutr_1 * d_1) * (A * del_1 - D_1) * np.exp(d_1 / del_2) * (del_3 * mutr_2 * D_2 * K_2 - D_3 * (K_3 * del_3 * mutr_3 + K_2 - K_3)) * (D_1 * del_2 + D_2 * del_1) * np.exp(mutr_2 * (-d_2 + d_1)) - np.exp(d_2 / del_2) * (D_2 * del_3 + D_3 * del_2) * (2 * del_1 * D_1 * K_1 * (D_1 * mutr_1 + A) * np.exp(d_1 / del_1) + ((K_1 * del_1 * mutr_1 - K_1 + K_2) * D_1 - del_1 * mutr_2 * D_2 * K_2) * np.exp(-mutr_1 * d_1) * (A * del_1 - D_1))) * np.exp(-d_1 / del_1) + np.exp(d_1 / del_1) * np.exp(-mutr_1 * d_1) * ((-D_1 * del_2 + D_2 * del_1) * np.exp(d_1 / del_2) * (del_3 * mutr_2 * D_2 * K_2 - D_3 * (K_3 * del_3 * mutr_3 + K_2 - K_3)) * np.exp(mutr_2 * (-d_2 + d_1)) + np.exp(d_2 / del_2) * (D_2 * del_3 + D_3 * del_2) * ((K_1 * del_1 * mutr_1 + K_1 - K_2) * D_1 - del_1 * mutr_2 * D_2 * K_2)) * (A * del_1 + D_1)) / (-(A * del_1 - D_1) * (np.exp(d_2 / del_2) * (D_2 * del_3 + D_3 * del_2) * (-D_1 * del_2 + D_2 * del_1) * np.exp(-d_1 / del_2) + np.exp(d_1 / del_2) * np.exp(-d_2 / del_2) * (-D_2 * del_3 + D_3 * del_2) * (D_1 * del_2 + D_2 * del_1)) * np.exp(-d_1 / del_1) + np.exp(d_1 / del_1) * (np.exp(d_2 / del_2) * (D_2 * del_3 + D_3 * del_2) * (D_1 * del_2 + D_2 * del_1) * np.exp(-d_1 / del_2) + np.exp(d_1 / del_2) * np.exp(-d_2 / del_2) * (-D_2 * del_3 + D_3 * del_2) * (-D_1 * del_2 + D_2 * del_1)) * (A * del_1 + D_1))

    A_4 = -del_2 * ((-(-D_2 * del_3 + D_3 * del_2) * (2 * del_1 * D_1 * K_1 * (D_1 * mutr_1 + A) * np.exp(d_1 / del_1) + ((K_1 * del_1 * mutr_1 - K_1 + K_2) * D_1 - del_1 * mutr_2 * D_2 * K_2) * np.exp(-mutr_1 * d_1) * (A * del_1 - D_1)) * np.exp(-d_2 / del_2) + np.exp(mutr_2 * (-d_2 + d_1)) * (-D_1 * del_2 + D_2 * del_1) * np.exp(-mutr_1 * d_1) * np.exp(-d_1 / del_2) * (A * del_1 - D_1) * (del_3 * mutr_2 * D_2 * K_2 - D_3 * (K_3 * del_3 * mutr_3 + K_2 - K_3))) * np.exp(-d_1 / del_1) + np.exp(d_1 / del_1) * np.exp(-mutr_1 * d_1) * (A * del_1 + D_1) * ((-D_2 * del_3 + D_3 * del_2) * ((K_1 * del_1 * mutr_1 + K_1 - K_2) * D_1 - del_1 * mutr_2 * D_2 * K_2) * np.exp(-d_2 / del_2) - np.exp(mutr_2 * (-d_2 + d_1)) * np.exp(-d_1 / del_2) * (del_3 * mutr_2 * D_2 * K_2 - D_3 * (K_3 * del_3 * mutr_3 + K_2 - K_3)) * (D_1 * del_2 + D_2 * del_1))) / (-(A * del_1 - D_1) * (np.exp(d_2 / del_2) * (D_2 * del_3 + D_3 * del_2) * (-D_1 * del_2 + D_2 * del_1) * np.exp(-d_1 / del_2) + np.exp(d_1 / del_2) * np.exp(-d_2 / del_2) * (-D_2 * del_3 + D_3 * del_2) * (D_1 * del_2 + D_2 * del_1)) * np.exp(-d_1 / del_1) + np.exp(d_1 / del_1) * (np.exp(d_2 / del_2) * (D_2 * del_3 + D_3 * del_2) * (D_1 * del_2 + D_2 * del_1) * np.exp(-d_1 / del_2) + np.exp(d_1 / del_2) * np.exp(-d_2 / del_2) * (-D_2 * del_3 + D_3 * del_2) * (-D_1 * del_2 + D_2 * del_1)) * (A * del_1 + D_1))

    A_5 = 2 * (((-np.exp(-mutr_1 * d_1) * ((K_2 * del_2 * mutr_2 - K_2 + K_3) * D_2 - del_2 * mutr_3 * D_3 * K_3) * (A * del_1 - D_1) * np.exp(d_1 / del_2) * (D_1 * del_2 + D_2 * del_1) * np.exp(mutr_2 * (-d_2 + d_1)) / 2 - np.exp(d_2 / del_2) * D_2 * del_2 * (2 * del_1 * D_1 * K_1 * (D_1 * mutr_1 + A) * np.exp(d_1 / del_1) + ((K_1 * del_1 * mutr_1 - K_1 + K_2) * D_1 - del_1 * mutr_2 * D_2 * K_2) * np.exp(-mutr_1 * d_1) * (A * del_1 - D_1))) * np.exp(-d_1 / del_1) + np.exp(d_1 / del_1) * ((-D_1 * del_2 + D_2 * del_1) * ((K_2 * del_2 * mutr_2 - K_2 + K_3) * D_2 - del_2 * mutr_3 * D_3 * K_3) * np.exp(d_1 / del_2) * np.exp(mutr_2 * (-d_2 + d_1)) / 2 + np.exp(d_2 / del_2) * D_2 * del_2 * ((K_1 * del_1 * mutr_1 + K_1 - K_2) * D_1 - del_1 * mutr_2 * D_2 * K_2)) * np.exp(-mutr_1 * d_1) * (A * del_1 + D_1)) * np.exp(-d_2 / del_2) + np.exp(d_2 / del_2) * ((K_2 * del_2 * mutr_2 + K_2 - K_3) * D_2 - del_2 * mutr_3 * D_3 * K_3) * np.exp(-mutr_1 * d_1) * np.exp(-d_1 / del_2) * np.exp(mutr_2 * (-d_2 + d_1)) * (-(-D_1 * del_2 + D_2 * del_1) * (A * del_1 - D_1) * np.exp(-d_1 / del_1) + np.exp(d_1 / del_1) * (D_1 * del_2 + D_2 * del_1) * (A * del_1 + D_1)) / 2) * del_3 / (np.exp(d_2 / del_2) * (D_2 * del_3 + D_3 * del_2) * (-(-D_1 * del_2 + D_2 * del_1) * (A * del_1 - D_1) * np.exp(-d_1 / del_1) + np.exp(d_1 / del_1) * (D_1 * del_2 + D_2 * del_1) * (A * del_1 + D_1)) * np.exp(-d_1 / del_2) + np.exp(d_1 / del_2) * (-(D_1 * del_2 + D_2 * del_1) * (A * del_1 - D_1) * np.exp(-d_1 / del_1) + np.exp(d_1 / del_1) * (-D_1 * del_2 + D_2 * del_1) * (A * del_1 + D_1)) * np.exp(-d_2 / del_2) * (-D_2 * del_3 + D_3 * del_2)) / np.exp(-d_2 / del_3)


    return (A_1, A_2, A_3, A_4, A_5)
