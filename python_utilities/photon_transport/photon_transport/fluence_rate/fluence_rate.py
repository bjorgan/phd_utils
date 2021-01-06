"""
Solution for fluence rates.

The solutions in the two-layered case will fail if mua = musr/2, or close to
this.
"""

import numpy as np
import photon_transport

def two_layers_isotropic(x, mua_1, mua_2, musr_1, musr_2, d_1 = 100e-06, return_as_separate_terms=False):
    """
    Fluence rate in a two-layered medium, assuming isotropic source functions.
    From L. O. Svaasand et al., "Tissue parameters determining the visual
    appearance of normal skin and port-wine stains", Laser Med. Sci., 10 (1995)
    """
    return two_layers_isotropic_coefficient_based(x, mua_1, mua_2, musr_1, musr_2, d_1, use_analytic_coefficients=True, return_as_separate_terms=return_as_separate_terms)

def two_layers_isotropic_from_skin_model(x, skin_model, return_as_separate_terms=False, wavelength_index=None):
    mua_1, mua_2, musr_1, musr_2, d_1 = photon_transport.forward_model.skin_model_to_optical_parameters(skin_model)
    if wavelength_index is not None:
        return two_layers_isotropic(x, mua_1[wavelength_index], mua_2[wavelength_index], musr_1[wavelength_index], musr_2[wavelength_index], d_1, return_as_separate_terms)
    else:
        return two_layers_isotropic(x, mua_1, mua_2, musr_1, musr_2, d_1, return_as_separate_terms)

def two_layers_isotropic_analytic(x, mua_1, mua_2, musr_1, musr_2, d_1 = 100e-06):
    """
    Complete analytic solution to the fluence rate in two layers.

    Not necessarily more numerically stable than the other solutions close to singularities.

    Obtained from FluenceRate_isoL2.mw.
    """
    A = photon_transport.forward_model.A

    mutr_1 = musr_1 + mua_1
    mutr_2 = musr_2 + mua_2
    D_1 = 1.0/(3*mutr_1)
    D_2 = 1.0/(3*mutr_2)
    del_1 = np.sqrt(D_1/mua_1)
    del_2 = np.sqrt(D_2/mua_2)

    phi_1 = lambda x: del_1 * ((del_2 * mutr_2 + 1) * del_1 * (D_1 * del_2 - D_2 * del_1) * musr_1 * (del_1 * (D_1 * mutr_1 + A) * np.exp(x / del_1) - np.exp(-mutr_1 * x) * (A * del_1 - D_1)) * np.exp(-d_1 / del_1) + (del_1 ** 2 * musr_1 * (D_1 * del_2 + D_2 * del_1) * (del_2 * mutr_2 + 1) * (D_1 * mutr_1 + A) * np.exp(d_1 / del_1) - ((mutr_1 * ((-musr_1 * mutr_2 + musr_2 * mutr_1) * del_2 - musr_1) * del_2 * D_1 + D_2 * musr_1 * (del_2 * mutr_2 + 1)) * del_1 ** 2 - D_1 * del_2 ** 2 * musr_2) * (A * del_1 - D_1) * np.exp(-mutr_1 * d_1)) * np.exp(-x / del_1) + (-np.exp(d_1 / del_1) * np.exp(-mutr_1 * x) * del_1 * musr_1 * (D_1 * del_2 + D_2 * del_1) * (del_2 * mutr_2 + 1) + ((mutr_1 * ((-musr_1 * mutr_2 + musr_2 * mutr_1) * del_2 - musr_1) * del_2 * D_1 + D_2 * musr_1 * (del_2 * mutr_2 + 1)) * del_1 ** 2 - D_1 * del_2 ** 2 * musr_2) * np.exp(x / del_1) * np.exp(-mutr_1 * d_1)) * (A * del_1 + D_1)) / D_1 / (del_1 * mutr_1 - 1) / (del_2 * mutr_2 + 1) / ((D_1 * del_2 - D_2 * del_1) * (A * del_1 - D_1) * np.exp(-d_1 / del_1) + np.exp(d_1 / del_1) * (D_1 * del_2 + D_2 * del_1) * (A * del_1 + D_1)) / (del_1 * mutr_1 + 1)

    phi_2 = lambda x: 1 / (del_2 * mutr_2 + 1) * (-(del_1 * mutr_1 - 1) * np.exp(-mutr_1 * d_1) * ((D_1 * del_2 - D_2 * del_1) * (A * del_1 - D_1) * np.exp(-d_1 / del_1) + np.exp(d_1 / del_1) * (D_1 * del_2 + D_2 * del_1) * (A * del_1 + D_1)) * del_2 * np.exp(-d_1 / del_2) * musr_2 * (del_1 * mutr_1 + 1) * np.exp(-mutr_2 * (x - d_1)) + np.exp(-x / del_2) * ((2 * D_2 * del_1 ** 3 * musr_1 * (del_2 * mutr_2 - 1) * (del_2 * mutr_2 + 1) * (D_1 * mutr_1 + A) * np.exp(d_1 / del_1) - (del_1 * mutr_1 - 1) * (D_2 * (mutr_2 * (-musr_1 * mutr_2 + musr_2 * mutr_1) * del_2 ** 2 + musr_1) * del_1 ** 2 - del_2 ** 2 * musr_2 * (D_1 * mutr_1 - D_2 * mutr_2) * del_1 - D_1 * del_2 ** 2 * musr_2) * np.exp(-mutr_1 * d_1) * (A * del_1 - D_1)) * np.exp(-d_1 / del_1) + (D_2 * (mutr_2 * (-musr_1 * mutr_2 + musr_2 * mutr_1) * del_2 ** 2 + musr_1) * del_1 ** 2 + del_2 ** 2 * musr_2 * (D_1 * mutr_1 - D_2 * mutr_2) * del_1 - D_1 * del_2 ** 2 * musr_2) * np.exp(-mutr_1 * d_1) * (A * del_1 + D_1) * np.exp(d_1 / del_1) * (del_1 * mutr_1 + 1))) / (del_2 * mutr_2 - 1) * del_2 / (del_1 * mutr_1 - 1) / D_2 / ((D_1 * del_2 - D_2 * del_1) * (A * del_1 - D_1) * np.exp(-d_1 / del_1) + np.exp(d_1 / del_1) * (D_1 * del_2 + D_2 * del_1) * (A * del_1 + D_1)) / np.exp(-d_1 / del_2) / (del_1 * mutr_1 + 1)

    if hasattr(x, "__len__"):
        x_1 = x[x < d_1]
        x_2 = x[x >= d_1]

        phi = np.zeros(len(x))
        phi[x < d_1] = phi_1(x_1)
        phi[x >= d_1] = phi_2(x_2)
    else:
        if (x < d_1):
            phi = phi_1(x)
        else:
            phi = phi_2(x)

    return phi

def two_layers_isotropic_coefficient_based(x, mua_1, mua_2, musr_1, musr_2, d_1 = 100e-06, use_analytic_coefficients=False, return_as_separate_terms=False):
    """
    Finds the coefficients A_1, A_2 and A_3 and constructs the solution from
    these, either numerically or by given expressions.
    """
    A = photon_transport.forward_model.A

    #transmitted light density after specular reflection
    P0 = 1

    mutr_1 = musr_1 + mua_1
    mutr_2 = musr_2 + mua_2
    D_1 = 1.0/(3*mutr_1)
    D_2 = 1.0/(3*mutr_2)
    del_1 = np.sqrt(D_1/mua_1)
    del_2 = np.sqrt(D_2/mua_2)

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
    phi_term_layer_2 = lambda x: np.exp(-x/del_2)

    j_source_layer_2 = lambda x: -D_2*(-mutr_2*phi_source_layer_2(x))
    j_term_layer_2 = lambda x: -D_2*(-1.0/del_2*phi_term_layer_2(x))

    #find coefficients A_1, A_2, A_3 from boundary conditions
    #{...}_terms * [A_1, A_2, A_3].T = {...}_constant
    if not use_analytic_coefficients:
        #continuity of phi at x = d_1
        phi_continuity_terms = np.array([phi_term_1_layer_1(x=d_1),
                                        phi_term_2_layer_1(x=d_1),
                                        (-1)*phi_term_layer_2(x=d_1)])
        phi_continuity_constant = -phi_source_layer_1(x=d_1) + phi_source_layer_2(x=d_1)

        #continuity of j at x = d_1
        j_continuity_terms = np.array([j_term_1_layer_1(x=d_1),
                                        j_term_2_layer_1(x=d_1),
                                        (-1)*j_term_layer_2(x=d_1)])
        j_continuity_constant = -j_source_layer_1(x=d_1) + j_source_layer_2(x=d_1)

        #boundary condition at air-tissue surface: j = -A*phi (NB: Note negative sign, opposite in Svaasand's paper)
        air_tissue_terms = np.array([j_term_1_layer_1(x=0) + A*phi_term_1_layer_1(x=0),
                                        j_term_2_layer_1(x=0) + A*phi_term_2_layer_1(x=0),
                                        0])
        air_tissue_constant = -j_source_layer_1(x=0) - A*phi_source_layer_1(x=0)

        terms = np.stack((phi_continuity_terms,
                            j_continuity_terms,
                            air_tissue_terms), axis=0)
        constants = np.array([phi_continuity_constant,
                            j_continuity_constant,
                            air_tissue_constant])
        A_constants = np.linalg.solve(terms, constants)
        A_1 = A_constants[0]
        A_2 = A_constants[1]
        A_3 = A_constants[2]
    else:
        A_1, A_2, A_3 = analytic_coefficients(mua_1, mua_2, musr_1, musr_2, d_1)

    #full solution, but split into source term and solution terms
    phi_1 = lambda x: (phi_source_layer_1(x), A_1*phi_term_1_layer_1(x) + A_2*phi_term_2_layer_1(x))
    phi_2 = lambda x: (phi_source_layer_2(x), A_3*phi_term_layer_2(x))
    j_1 = lambda x: (j_source_layer_1(x), A_1*j_term_1_layer_1(x) + A_2*j_term_2_layer_1(x))
    j_2 = lambda x: (j_source_layer_2(x), A_3*j_term_layer_2(x))

    if not return_as_separate_terms:
        if hasattr(x, "__len__"):
            x_1 = x[x < d_1]
            x_2 = x[x >= d_1]

            phi = np.zeros(len(x))
            phi[x < d_1] = phi_1(x_1)[0] + phi_1(x_1)[1]
            phi[x >= d_1] = phi_2(x_2)[0] + phi_2(x_2)[1]

            j = np.zeros(len(x))
            j[x < d_1] = j_1(x_1)[0] + j_1(x_1)[1]
            j[x >= d_1] = j_2(x_2)[0] + j_2(x_2)[1]
        else:
            if (x < d_1):
                phi = phi_1(x)[0] + phi_1(x)[1]
                j = j_1(x)[0] + j_1(x)[1]
            else:
                phi = phi_2(x)[0] + phi_2(x)[1]
                j = j_2(x)[0] + j_2(x)[1]

        return phi, j

    else:
        if hasattr(x, "__len__"):
            x_1 = x[x < d_1]
            x_2 = x[x >= d_1]

            phi_source = np.zeros(len(x))
            phi_solution = np.zeros(len(x))

            phi_source[x < d_1] = phi_1(x_1)[0]
            phi_source[x >= d_1] = phi_2(x_2)[0]

            phi_solution[x < d_1] = phi_1(x_1)[1]
            phi_solution[x >= d_1] = phi_2(x_2)[1]

            j_source = np.zeros(len(x))
            j_solution = np.zeros(len(x))

            j_source[x < d_1] = j_1(x_1)[0]
            j_source[x >= d_1] = j_2(x_2)[0]

            j_solution[x < d_1] = j_1(x_1)[1]
            j_solution[x >= d_1] = j_2(x_2)[1]

            return phi_source, phi_solution, j_source, j_solution
        else:
            if (x < d_1):
                phi_source = phi_1(x)[0]
                phi_solution = phi_1(x)[1]

                j_source = j_1(x)[0]
                j_solution = j_1(x)[1]
            else:
                phi_source = phi_2(x)[0]
                phi_solution = phi_2(x)[1]

                j_source = j_2(x)[0]
                j_solution = j_2(x)[1]
            return phi_source, phi_solution, j_source, j_solution

def infinite_isotropic_medium(x, mua, musr):
    """
    Solution for infinite medium, with an infinitely broad source at x = 0

    From Wang et al., Biomedical Optics: Principles and imaging.
    """
    mutr = musr + mua
    D = 1.0/(3*mutr)
    pendepth = np.sqrt(D/mua)
    phi = 1/(2*pendepth*mua)*np.exp(-np.abs(x)/pendepth)
    j = -D*(-1/pendepth*phi)
    return phi, j


def analytic_coefficients(mua_1, mua_2, musr_1, musr_2, d_1=100e-06):
    """
    Analytic solutions for the cofficients A1, A2 and A3 in the solution for
    phi given in Svaasand1995.

    Used for comparison with numerically obtained coefficients in
    two_layers_isotropic_coefficient_based_solution() (they fare slightly
    better).

    See ReflIsoL2.mw in the penetration_depth folder for derivation.
    """
    A = photon_transport.forward_model.A

    mutr_1 = musr_1 + mua_1
    mutr_2 = musr_2 + mua_2
    D_1 = 1.0/(3*mutr_1)
    D_2 = 1.0/(3*mutr_2)
    del_1 = np.sqrt(D_1/mua_1)
    del_2 = np.sqrt(D_2/mua_2)

    K_1 = del_1*del_1*musr_1/(D_1*(1 - mutr_1**2 * del_1**2))
    K_2 = del_2*del_2*musr_2/(D_2*(1 - mutr_2**2 * del_2**2))

    A_1 = (-K_1 * (D_1 * del_2 + D_2 * del_1) * (D_1 * mutr_1 + A) * np.exp(d_1 / del_1) + np.exp(-mutr_1 * d_1) * (A * del_1 - D_1) * (-mutr_1 * D_1 * K_1 * del_2 + D_2 * (K_2 * del_2 * mutr_2 + K_1 - K_2))) * del_1 / ((D_1 * del_2 - D_2 * del_1) * (A * del_1 - D_1) * np.exp(-d_1 / del_1) + np.exp(d_1 / del_1) * (D_1 * del_2 + D_2 * del_1) * (A * del_1 + D_1))
    A_2 = -del_1 * (K_1 * (D_1 * del_2 - D_2 * del_1) * (D_1 * mutr_1 + A) * np.exp(-d_1 / del_1) + (A * del_1 + D_1) * np.exp(-mutr_1 * d_1) * (-mutr_1 * D_1 * K_1 * del_2 + D_2 * (K_2 * del_2 * mutr_2 + K_1 - K_2))) / ((D_1 * del_2 - D_2 * del_1) * (A * del_1 - D_1) * np.exp(-d_1 / del_1) + np.exp(d_1 / del_1) * (D_1 * del_2 + D_2 * del_1) * (A * del_1 + D_1))
    A_3 = -((2 * D_1 * K_1 * del_1 * (D_1 * mutr_1 + A) * np.exp(d_1 / del_1) - (A * del_1 - D_1) * np.exp(-mutr_1 * d_1) * ((-K_1 * del_1 * mutr_1 + K_1 - K_2) * D_1 + mutr_2 * D_2 * K_2 * del_1)) * np.exp(-d_1 / del_1) + (A * del_1 + D_1) * ((-K_1 * del_1 * mutr_1 - K_1 + K_2) * D_1 + mutr_2 * D_2 * K_2 * del_1) * np.exp(-mutr_1 * d_1) * np.exp(d_1 / del_1)) * del_2 / np.exp(-d_1 / del_2) / ((D_1 * del_2 - D_2 * del_1) * (A * del_1 - D_1) * np.exp(-d_1 / del_1) + np.exp(d_1 / del_1) * (D_1 * del_2 + D_2 * del_1) * (A * del_1 + D_1))


    return (A_1, A_2, A_3)

def one_layers_isotropic(x, mua, musr):
    """
    Fluence rate in a one-layered, semi-infinite medium.

    See results.tex for derivation.
    """
    A = photon_transport.forward_model.A
    mutr = musr + mua
    D = 1.0/(3*mutr)
    delta = np.sqrt(D/mua)
    source_term = (delta**2 * musr)/(D*(1-mutr**2 * delta**2))
    exp_term_1 = np.exp(-mutr*x)
    exp_term_2 = -np.exp(-x/delta)*(D*mutr + A)/(D/delta + A)

    phi = source_term*(exp_term_1 + exp_term_2)
    j = -D*source_term*(-mutr*exp_term_1 - 1.0/delta * exp_term_2)
    return phi, j
