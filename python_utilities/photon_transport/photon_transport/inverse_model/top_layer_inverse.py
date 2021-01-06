"""
Basic inverse model for estimating top layer properties given the reflectance
from a region without the top layer.

See https://doi.org/10.1364/BOE.399636.
"""

import scipy.optimize
from photon_transport.optimization.skin_model import get_musr_divided_by_mua, get_musr
from photon_transport.forward_model import ReflIsoL2, melanin_models, get_scattering_spectra, normal_skin_scattering, A_default
import numpy as np

class inverse_model:
    def __init__(self, wlens, fix_parameters=[], default_parameters={}, A = A_default):
        self.wlens = wlens
        self.fix_parameters = fix_parameters.copy()
        self.default_parameters = default_parameters.copy()
        self.A = A

    def add_epidermis(self, spectrum, muam694=300, musrr500=17.6, musmr500=18.78, b_mie=0.22, d_1=100e-06):
        #fit arbitrary onelayer model to wound spectrum: can be shown that
        #two-layer model yields same reflectance as long as mua/musr == constant
        #(see plot_[diffusion_model/monte_carlo]_scale_invariance_results.py)
        #Therefore need to find the same ratio, and use that in dermis of a new two-layer
        #model with epidermis added.
        C = get_musr_divided_by_mua(spectrum, A=self.A)[0]
        musr = get_musr(self.wlens, 1.1)
        mua = musr/C

        #construct twolayer model
        mua_epi = muam694*melanin_models(self.wlens).svaasand.values
        musr_epi = normal_skin_scattering(get_scattering_spectra(self.wlens, b_mie=b_mie), musrr500=musrr500, musmr500=musmr500)[0]

        refl = ReflIsoL2(mua_epi, mua, musr_epi, musr, d_1, A=self.A)
        return refl, mua_epi, mua, musr_epi, musr

    def fit_spectrum(self, wound, healed, start_ind=0, end_ind=60, pick_start_params_randomly=False):
        """
        Fit 'healed' spectrum by fitting an epidermis on top of 'wound'.

        Parameters
        ----------
        wound: ndarray
            Spectrum from wound region.
        healed: ndarray
            Spectrum from healed region.
        start_ind: int, optional
            Fit start index
        end_ind: int, optional
            Fit end index
        pick_start_params_randomly: boolean, optional
            Pick start parameters for fitting randomly.

        Returns
        -------
        ret_params: dict
            Fitted parameters.
        fitted_refl: ndarray
            Fitted reflectance
        mua_derm_fitted: ndarray
            Fitted dermal properties
        skin_model: array
            Skin model, as returned from forward_model.normal_skin().
        """

        def objective(params):
            refl = self.add_epidermis(wound, muam694=params[0], d_1=params[1], b_mie=params[2], musmr500=params[3], musrr500=params[4])[0]

            sse = (healed - refl)**2/refl
            rmse = np.mean(sse[start_ind:end_ind])
            return rmse, refl

        #fit parameters to epidermis added on top of dermis based on wound spectrum
        scaling = np.array([1000, 500e-06, 4, 50, 50])
        bounds = np.array([(0, 2000), (0, 500e-06), (0, 4.0), (0, 50), (1, 50)])

        index_mapping = {'muam694': 0, 'd_1': 1, 'b_mie': 2, 'musmr500': 3, 'musrr500': 4}
        start_params = np.array([300, 100e-06, 0.22, 18, 18])

        #modify start parameters to specified parameters
        for param_name in list(self.default_parameters):
            start_params[index_mapping[param_name]] = self.default_parameters[param_name]

        #modify bounds to around start parameter for fixed parameters
        for param_name in self.fix_parameters:
            bounds[index_mapping[param_name]] = [start_params[index_mapping[param_name]]]*2

        #pick start parameters within bounds
        if pick_start_params_randomly:
            start_params = np.random.random_sample(len(bounds))
            start_params = [start_params[i]*(bound[1]-bound[0]) + bound[0] for i, bound in enumerate(bounds)]

        bounds = bounds/scaling[:, None]

        res = scipy.optimize.minimize(fun=lambda params: objective(params*scaling)[0], x0=start_params/scaling, bounds=bounds)

        #fitted parameters
        params = res['x']*scaling
        ret_params = {}
        ret_params['muam694'] = params[0]
        ret_params['d_1'] = params[1]
        ret_params['b_mie'] = params[2]
        ret_params['musmr500'] = params[3]
        ret_params['musrr500'] = params[4]

        #get forward optical properties and reflectance
        fitted_refl, mua_epi, mua_derm, musr_epi, musr_derm = self.add_epidermis(wound, muam694=params[0], d_1=params[1], b_mie=params[2], musmr500=params[3], musrr500=params[4])

        #construct skin model out of this
        g = get_scattering_spectra(self.wlens)['g'].values
        skin_model = [{'mua': mua_epi, 'musr': musr_epi, 'mus': musr_epi/(1-g), 'g': g, 'n': np.repeat(1.4, len(self.wlens)), 'd': params[1]},
                      {'mua': mua_derm, 'musr': musr_derm, 'mus': musr_derm/(1-g), 'g': g, 'n': np.repeat(1.4, len(self.wlens)), 'd': 100000}]

        #get mua_epi when other parameters are kept fixed
        mua_derm_fitted = [scipy.optimize.minimize_scalar(lambda x: (ReflIsoL2(x, mua_derm[i], musr_epi[i], musr_derm[i], ret_params['d_1'], A=self.A) - healed[i])**2).x for i in range(len(mua_epi))]
        return ret_params, fitted_refl, mua_derm_fitted, skin_model
