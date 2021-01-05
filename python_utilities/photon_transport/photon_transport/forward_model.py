"""
Photon transport using the diffusion approximation, forward model-related
methods. normal_skin(...) and get_reflectance(...) are the most likely entry
methods to use.
"""

import numpy as np
import pandas as pd
import os.path
from os import listdir
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d

from .diffusion_model_solutions import ReflIsoL1, ReflIsoL2, ReflIsoL3, ReflE2L2, ReflE2L3, A_default
import enum

def get_absorption_spectra(wavelength, interpolation_kind='cubic'):
    """
    Obtain optical absorption coefficients at discrete wavelengths for all
    files located in opt_params/*.dat, and melanin.

    Parameters
    ----------
    wavelength: list
        Wavelength array
    Returns
    -------
    absorption: pd.DataFrame
        Pandas DataFrame with each column corresponding to the absorption
        spectrum for a specific material, in addition to one column for the
        wavelength values.
    """
    absorption_files_path = os.path.dirname(os.path.abspath(__file__)) + "/opt_params/"

    #T. Spott, L. O. Svaasand, R. E. Anderson, and P. F. Schmedling. Application of optical diusion
    #theory to transcutaneous bilirubinometry. In Proc. SPIE 3195, Laser-Tissue Interaction, Tissue
    #Optics, and Laser Welding III, 234 , volume 3195, pages 234245, 1998.
    melanin_absorption = (694.0/wavelength)**3.46

    #construct absorption coefficient matrix out of files located in opt_params/
    absorption = pd.DataFrame({'wavelength': wavelength, 'constant': np.repeat(1, len(wavelength)), 'melanin': melanin_absorption})
    files = [f for f in listdir(absorption_files_path) if os.path.isfile(os.path.join(absorption_files_path, f))]
    for f in files:
        value_matrix = np.loadtxt(os.path.join(absorption_files_path, f))
        absorption[f.replace('.dat', '')] = interp1d(value_matrix[:,0], value_matrix[:,1], kind=interpolation_kind, bounds_error=False)(wavelength)

    return absorption

def musrr_fun(wavelength):
    """
    Rayleigh scattering.

    Parameters
    ----------
    wavelength: float or ndarray
        Wavelength
    Returns
    -------
    musrr: ndarray
        Rayleigh scattering wavelength dependency
    """
    return (wavelength/500.0)**float(-4.0)

#Obtained from A. N. Bashkatov, E. A. Genina, V. V. Tuchin. Optical properties
#of skin, subcutaneous and muscle tissues: A review. J. Inoov. Opt. Health
#Sci., 4(1):9-38, 2011
b_mie_bashkatov = 0.22

def musmr_fun(wavelength, b_mie = b_mie_bashkatov):
    """
    Mie scattering.

    Parameters
    ----------
    wavelength: float or ndarray
        Wavelength
    b_mie: float
        Power coefficient
    Returns
    -------
    musmr: ndarray
        Mie scattering wavelength dependency
    """
    return (wavelength/500.0)**float(-b_mie)

def get_scattering_spectra(wavelength, b_mie=b_mie_bashkatov):
    """
    Obtain optical scattering coefficients (and related properties) at
    input wavelengths.

    Parameters
    ----------
    wavelength: ndarray
        Wavelengths
    Returns
    -------
    scattering: pd.DataFrame
        Pandas dataframe with each column corresponding to a
        scattering property
    """
    scattering = pd.DataFrame({'wavelength': wavelength})

    #Average cosine.  M. J. C. van Gemert, S. L. Jacques, H. J. C. M.
    #Sterenborg, and W. M. Star. Skin optics.  IEEE Trans Biomed Eng ,
    #36:11461154, 1989.
    scattering['g'] = 0.62 + wavelength*29e-5

    scattering['rayleigh_reduced'] = musrr_fun(wavelength)
    scattering['mie_reduced'] = musmr_fun(wavelength, b_mie)

    #bulk refraction index
    scattering['n'] = np.repeat(1, len(wavelength))

    return scattering

def num_layers(skin_model):
    """
    Get number of layers in a skin model.
    """

    return len(skin_model)

def skin_model_to_optical_parameters(skin_model):
    """
    Get the optical parameters iconstituting a skin model.
    """

    if num_layers(skin_model) == 1:
        return skin_model[0]['mua'], skin_model[0]['musr']
    if num_layers(skin_model) == 2:
        return skin_model[0]['mua'], skin_model[1]['mua'], skin_model[0]['musr'], skin_model[1]['musr'], skin_model[0]['d']
    if num_layers(skin_model) == 3:
        return skin_model[0]['mua'], skin_model[1]['mua'], skin_model[2]['mua'], skin_model[0]['musr'], skin_model[1]['musr'], skin_model[2]['musr'], skin_model[0]['d'], skin_model[1]['d']

class source_function(enum.Enum):
    """
    Used for selecting source function in get_reflectance().
    """
    ISOTROPIC = 1
    DELTA_EDDINGTON = 2

def get_reflectance(skin_properties, source=SourceFunction.ISOTROPIC, A=A_default):
    """
    Obtain reflectance given the skin model structure as obtained from
    normal_skin(...) or similar.

    Parameters
    ----------
    skin_properties: list of dicts containing 'mua', 'musr', ...
        Skin model
    Returns
    -------
    refl: ndarray
        Simulated reflectance spectrum
    """
    if source == SourceFunction.ISOTROPIC:
        #isotropic source functions
        if num_layers(skin_properties) == 1:
            mua, musr = skin_model_to_optical_parameters(skin_properties)
            return ReflIsoL1(mua, musr, A)
        if num_layers(skin_properties) == 2:
            mua1, mua2, musr1, musr2, d_1 = skin_model_to_optical_parameters(skin_properties)
            return ReflIsoL2(mua1, mua2, musr1, musr2, d_1, A)
        if num_layers(skin_properties) == 3:
            mua1, mua2, mua3, musr1, musr2, musr3, d_1, d_2 = skin_model_to_optical_parameters(skin_properties)
            return ReflIsoL3(mua1, mua2, mua3, musr1, musr2, musr3, d_1, d_2, A)
    elif source == SourceFunction.DELTA_EDDINGTON:
        #delta-eddington

        #optical properties first layer
        mua1 = skin_properties[0]['mua']
        mus1 = skin_properties[0]['mus']
        g1 = skin_properties[0]['g']
        d_1 = skin_properties[0]['d']

        #second layer
        if num_layers(skin_properties) > 1:
            mua2 = skin_properties[1]['mua']
            mus2 = skin_properties[1]['mus']
            g2 = skin_properties[1]['g']

        #third layer
        if num_layers(skin_properties) > 2:
            d_2 = skin_properties[1]['d']
            mua3 = skin_properties[2]['mua']
            mus3 = skin_properties[2]['mus']
            g3 = skin_properties[2]['g']

        if num_layers(skin_properties) == 2:
            return ReflE2L2(mua1, mua2, mus1, mus2, g1, g2, d_1, A)
        if num_layers(skin_properties) == 3:
            return ReflE3L3(mua1, mua2, mua3, mus1, mus2, mus3, g1, g2, g3, d_1, d_2, A)

def physical_scattering_parameters_to_coefficients(musr_500, f_rayleigh):
    """
    Convert scattering parameters in terms of scattering coeff at 500 nm and
    fraction of rayleigh scattering to coefficients as expected by normal_skin.
    """

    musrr_500 = musr_500*f_rayleigh
    musmr_500 = musr_500*(1-f_rayleigh)
    return musmr_500, musrr_500

#Rayleigh scattering coefficient at 500 nm in cm-1, recalculated from A. N.
#Bashkatov, E. A. Genina, V. V. Tuchin. Optical properties of skin,
#subcutaneous and muscle tissues: A review. J. Inoov. Opt. Health Sci.,
#4(1):9-38, 2011 to a (wlen/500)-formulation)
musrr500_bashkatov = 17.6

#Mie scattering coefficient at 500 nm in cm-1 (reference above)
musmr500_bashkatov = 18.78

#default mua_other constant absorption
mua_other_default=25

def normal_skin_scattering(scattering_spectra, musrr500=musrr500_bashkatov, musmr500=musmr500_bashkatov):
    """
    Calculate scattering properties for normal skin.

    Parameters
    ----------
    scattering_spectra: pd.DataFrame
        Obtained from get_scattering_spectra() above.
    musrr500: float, optional
        Rayleigh scattering factor.
    musmr500: float, optional
        Mie scattering factor.
    Returns
    -------
    musr: ndarray
        Reduced scattering coefficient
    musr: ndarray
        Scattering coefficient
    g: ndarray
        Anisotropy factor
    n: ndarray
        Refractive index
    """
    musr = 100*(scattering_spectra.rayleigh_reduced * musrr500 + scattering_spectra.mie_reduced*musmr500).values
    mus = musr/(1-scattering_spectra.g).values
    g = scattering_spectra.g.values
    refr_ind = scattering_spectra.n.values
    return musr, mus, g, refr_ind

def skin_epidermal_layer(absorption_spectra, scattering_spectra, muam694, musmr500=musmr500_bashkatov, musrr500=musrr500_bashkatov, include_blood=True, d_epidermis=100e-06):
    """
    Construct an epidermal layer.
    """

    #absorption
    mua_epidermis = absorption_spectra.melanin.values * muam694

    #epidermal blood to account for variation in papillary dermis
    if include_blood:
        mub_oxy = absorption_spectra.muabo.values
        mub_deoxy = absorption_spectra.muabd.values
        bvf_1 = 0.002
        oxy_1 = 0.6

        mua_epidermis += (mub_oxy*oxy_1 + mub_deoxy*(1-oxy_1))*bvf_1

    #scattering
    musr, mus, g, refr_ind = normal_skin_scattering(scattering_spectra, musrr500, musmr500)

    return {
        'mua': mua_epidermis,
        'musr': musr,
        'mus': mus,
        'd': d_epidermis,
        'n': refr_ind,
        'g': g
        }

def skin_dermal_layer(absorption_spectra, scattering_spectra, oxy, bvf, musmr500=musmr500_bashkatov, musrr500=musrr500_bashkatov, mua_other=mua_other_default, end_z_position=1):
    """
    Construct a standard dermal layer.
    """

    ##scattering
    musr, mus, g, refr_ind = normal_skin_scattering(scattering_spectra, musrr500, musmr500)

    ##absorption

    #blood
    muabo = absorption_spectra.muabo.values*oxy
    muabd = absorption_spectra.muabd.values*(1-oxy)
    mua = (muabo + muabd)*bvf

    #background absorption
    mua += mua_other

    return {
        'mua': mua,
        'mus': mus,
        'musr': musr,
        'd': end_z_position,
        'n': refr_ind,
        'g': g
    }


def normal_skin(absorption_spectra, scattering_spectra, muam694, oxy, bvf, musrr500=musrr500_bashkatov, musmr500=musmr500_bashkatov, mua_other_dermis_param=mua_other_default, d_epidermis=100e-06, include_blood_in_epidermis=True):
    """
    Generate a skin model structure containing constructed absorption and
    scattering spectra for simulation of reflectance.

    Return structure is a list containing two elements, one for epidermis and
    one for dermis. Each element is a structure containing fields mua, mus,
    musr, d, n and g.

    Parameters
    ----------
    absorption_spectra: pd.DataFrame
        Absorption spectra, as obtained from get_absorption_spectra(...)
    scattering_spectra: pd.DataFrame
        Scattering spectra, as obtained from get_scattering_spectra(...)
    muam694: float
        Melanin content in epidermis
    oxy: float
        Oxygenation in dermis
    bvf: float
        Blood volume fraction in dermis (NOTE: In inverse models, a formulation
        based on f_oxy and f_deoxy (where bvf = f_oxy + f_deoxy) can be more
        convenient)
    musrr: float
        Rayleigh scattering coefficient at 500 nm in cm-1
    musmr: float
        Mie scattering coefficient at 500 nm in cm-1
    mua_other_dermis_param: float
        Various background chromophores, expressed as a constant absorption coefficient
    d_epidermis: float
        Thickness of epidermis in meters
    include_blood_in_epidermis: boolean
        Whether blood should be included in epidermis
    Returns
    -------
    skin_model:
        Skin model, ready for input in get_reflectance(...).
    """

    skin_properties = [{}, {}]

    #epidermis
    skin_properties[0] = skin_epidermal_layer(absorption_spectra, scattering_spectra, muam694, musmr500, musrr500, include_blood_in_epidermis, d_epidermis)

    #dermis
    skin_properties[1] = skin_dermal_layer(absorption_spectra, scattering_spectra, oxy, bvf, musmr500=musmr500, musrr500=musrr500, mua_other=mua_other_dermis_param)

    return skin_properties

def normal_skin_threelayers(absorption, scattering, muam694, oxys, bvfs, thicknesses, mua_other=[mua_other_default, mua_other_default], musrr500=musrr500_bashkatov, musmr500=musmr500_bashkatov):
    """
    Create three-layered skin model.

    Parameters
    ----------
    absorption:
        Absorption spectra
    scattering:
        Scattering spectra
    muam694:
        Melanin absorption in epidermis
    oxys:
        Oxygenation in each layer
    bvfs:
        Blood volume fraction in each layer
    thicknesses:
        Thicknesses of the first and second layers
    mua_other:
        Constant absorption
    musrr500:
        Rayleigh coefficient, cm-1
    musmr500:
        Mie coefficient, cm-1
    Returns
    -------
    skin_model:
        Structure describing the absorption and scattering properties of the skin model, similar to normal_skin() but describing three layers instead of two.
    """

    #get absorption and scattering properties
    epidermis_and_upper_dermis = normal_skin(absorption, scattering, muam694, oxys[0], bvfs[0], d_epidermis=thicknesses[0], mua_other_dermis_param=mua_other[0], musrr500=musrr500, musmr500=musmr500)
    lower_dermis = normal_skin(absorption, scattering, muam694, oxys[1], bvfs[1], mua_other_dermis_param=mua_other[1], musrr500=musrr500, musmr500=musmr500)

    #construct skin model
    skin_properties = [{}, {}, {}]
    skin_properties[0] = epidermis_and_upper_dermis[0]
    skin_properties[1] = epidermis_and_upper_dermis[1]
    skin_properties[2] = lower_dermis[1]

    #set correct depth of third layer
    #NB: This is not the actual thickness of the second layer, but the
    #depth of the third layer.
    skin_properties[1]['d'] = thicknesses[0] + thicknesses[1]
    return skin_properties

def melanin_models(wlens):
    return pd.DataFrame({
                        'wlens': wlens,
                        'svaasand': (wlens/694)**(-3.46),
                        'eumelanin': np.exp(-2.429*(wlens - 694.0)/694.0), #Zonios et al., see master thesis for proper reference
                        'pheomelanin': np.exp(-4.780*(wlens - 694.0)/694.0), #Zonios et al., see above
                        'jacques': (wlens/694)**(-3.48)
                        }).set_index('wlens')

def jacques_background(wlens):
    """
    Jacques' baseline absorption in m-1.
    """

    return (0.244 + 85.3*np.exp(-(wlens - 154.0)/66.2))*100
