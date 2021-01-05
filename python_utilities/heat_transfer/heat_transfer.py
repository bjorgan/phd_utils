"""
Heat transfer functions for modeling heat transport and thermal damage through
a simple human tissue model.
"""

import numpy as np
from scipy.special import erfc

#human thermal conductivity, Log2017 (ref. to something else)
K_epidermis = 0.22 #W/(m K)
K_dermis = 0.40
K_subcutis = 0.20
K_muscle = 0.45

#Johnson2011
rho = 1200 #kg/m^3
c = 3600 #J/(kg K)

#ref here is https://www.engineeringtoolbox.com/water-liquid-gas-thermal-conductivity-temperature-pressure-d_2012.html at temp 100
K_water_liquid = 677.03*1.0e-03 #W/(m K)
K_water_gas = 24.57*1.0e-03

H_air = 10 #W/(K m^2), Johnson2011

def get_kappa(K, rho, c):
    return K/(rho*c)

def linear_heat_transfer(x, t, K_medium, rho_medium, c_medium, T_medium_initial, H_heat_transfer, T_external_applied):
    """
    Linear heat transfer H_heat_transfer at the boundary between medium at
    temperature T_external_applied and medium at initial temperature
    T_medium_initial and thermal conductivity K_medium.

    Parameters
    ----------
    x:
        Depth into the tissue.
    t:
        Time.
    K_medium:
        Thermal conductivity of the tissue.
    H_heat_transfer:
        Heat transfer coefficient from slab to the tissue.
    T_medium_initial:
        Initial tissue temperature.
    T_external_applied:
        Temperature of slab.
    """
    k = get_kappa(K_medium, rho_medium, c_medium)

    h = H_heat_transfer/K_medium
    erfc_factor_1 = erfc(x/(2*np.sqrt(k*t)))

    #combine factors in logdomain, since the exp-factors quickly approach
    #infinity while erfc-factor goes to zero
    log_exp_factor_1 = h*x
    log_exp_factor_2 = k*t*h**2
    log_erfc_factor_2 = np.log(erfc(x/(2*np.sqrt(k*t)) + h*np.sqrt(k*t)))
    exp_erfc_factor = np.exp(log_exp_factor_1 + log_exp_factor_2 + log_erfc_factor_2)

    return (erfc_factor_1 - exp_erfc_factor)*(T_external_applied - T_medium_initial) + T_medium_initial

def prescribed_surface_temperature(x, t, K_medium, rho_medium, c_medium, T_medium_initial, T_external_applied):
    """
    Medium at initial temp. T_medium_initial, held at constant temperature
    T_external_applied at the surface.  Same as H -> oo in
    linear_heat_transfer().
    """
    k = get_kappa(K_medium, rho_medium, c_medium)
    return (T_external_applied - T_medium_initial)*erfc(x/(2*np.sqrt(k*t))) + T_medium_initial

from scipy.integrate import quad
def damage_integral(T_func, t_end, delta_E = 6.28e05, P = 3.1e98):
    """
    Arrhenius damage integral, Henriques.

    Parameters
    ----------
    T_func: function
        Temperatures in degrees Celcius as a function of time only.
    t_end: float
        End time for damage integral calculation.
    delta_E: float
        Injury parameter, also given as E_a in the literature. Units J/mol,
        default value from Johnson2011, references Henriques (e08 J/_k_mol in
        the paper)
    P: float
        Injury parameter, also given as A in the literature. Units 1/s, default
        value from Log2017, references Henriques.
    Returns
    -------
    omega: float
        Damage measure. Expresses ln(C(0)/C(tau)), logarithm of the ratio of
        the original concentration of native tissue to the remaining native
        tissue at time tau.  Omega = 1 can be considered a threshold for
        coagulation/denaturation/irreversible damage.
    """

    kelvin_offset = 273.15
    R = 8.314 #J/mol K, molar gas constant
    damage_integrand = lambda t: np.exp(-delta_E/(R*(T_func(t) + kelvin_offset)) + np.log(P))
    return quad(damage_integrand, 0, t_end)[0]

def composite_solid_prescribed_surface_temperature(x, t, layer_1_props, layer_2_props, T_initial, T_external, num_terms=1):
    """
    Series solution to temperature development in composite solid, constant
    surface temperature.  Carslaw \cite{Carslaw} ch. 12.8, page 320: Laplace
    solution with series expansion and then inverse Laplace transform from
    tables. (Solution verified by hand, 2020-02-12.)
    """

    l = layer_1_props['d']
    kappa_1 = get_kappa(layer_1_props['K'], layer_1_props['rho'], layer_1_props['c'])
    kappa_2 = get_kappa(layer_2_props['K'], layer_2_props['rho'], layer_2_props['c'])

    k = np.sqrt(kappa_1/kappa_2)
    sigma = layer_2_props['K']*k/layer_1_props['K']
    alpha = (sigma-1)/(sigma+1)

    #translate coordinate system so that boundary is at x = 0
    x = x - l

    def v_1_term(n):
        return alpha**n * (erfc(((2*n+1)*l + x)/(2*np.sqrt(kappa_1*t))) - alpha*erfc(((2*n+1)*l - x)/(2*np.sqrt(kappa_1*t))))
    def v_2_term(n):
        return alpha**n * erfc(((2*n+1)*l + k*x)/(2*np.sqrt(kappa_1*t)))

    solution = 0
    if x < 0:
        solution = np.sum([v_1_term(n) for n in range(num_terms)])
    else:
        solution = 2/(1+sigma) * np.sum([v_2_term(n) for n in range(num_terms)])

    return solution*(T_external - T_initial) + T_initial

def layer(K, rho, c, d=None):
    return {'K': K, 'rho': rho, 'c': c, 'd': d}

def cool_down(x, t, layer_props, start_temperature_profile, H=H_air, T_ambient=20):
    """
    Green's function-solution: Initial temperature f(x), radiation at surface
    into medium at phi(t) (= constant here).  Carslaw, chap. 14.2, page 359.
    Integrate to t and x.

    Does not seem to work that well for t < 0.5 due to numerical problems.
    """

    kappa = get_kappa(layer_props['K'], layer_props['rho'], layer_props['c'])
    h = H/layer_props['K']

    integrand_x = lambda x_prime: 1/(2*np.sqrt(np.pi*kappa*t))*(np.exp(-(x-x_prime)**2/(4*kappa*t)) \
                                  + np.exp(-(x + x_prime)**2/(4*kappa*t))) \
                                  - h*np.exp(kappa*t*h**2 + h*(x+x_prime) + np.log(erfc((x+x_prime)/(2*np.sqrt(kappa*t)) + h*np.sqrt(kappa*t)))) #np.exp(... + log(...)) in order to avoid some numerical trouble

    integrand_t = lambda tau: kappa*h*(np.exp(-x**2/(4*kappa*(t-tau)))/np.sqrt(np.pi*kappa*(t-tau)) \
                                       - h*np.exp(kappa*h**2*(t - tau) + h*x)*erfc(x/(2*np.sqrt(kappa*(t-tau))) + h*np.sqrt(kappa*(t-tau))))


    x_integral = quad(lambda x: integrand_x(x)*start_temperature_profile(x), 0, np.inf)
    t_integral = quad(lambda t: integrand_t(t)*T_ambient, 0, t)
    return x_integral[0] + t_integral[0]
