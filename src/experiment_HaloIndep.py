# -*- coding: utf-8 -*-
"""
Created on Sun Mar  1 19:51:42 2015

@author: Andreea
"""

from __future__ import division
from __future__ import absolute_import
from __future__ import print_function
from experiment import *
import parallel_map as par


class Experiment_HaloIndep(Experiment):
    ''' This is the base class that implements the halo-independent analysis
    common to all experiments.
    Parameters:
        expername: string
            The name of the experiment.
        scattering_type: string
            The type of scattering. Can be
                - 'SI' (spin-independent)
                - 'SDAV' (spin-independent, axial-vector)
                - 'SDPS' (spin-independent, pseudo-scalar)
        mPhi: float, optional
            The mass of the mediator.
    '''
    def __init__(self, expername, scattering_type, mPhi=mPhiRef):
        Experiment.__init__(self, expername, scattering_type, mPhi)
        if self.energy_resolution_type == "Dirac":
            self.IntegratedResponse = self.IntegratedResponse_Dirac
        else:
            self.IntegratedResponse = self.IntegratedResponse_Other

    def DifferentialResponse(self, Eee, qER, const_factor):
        ''' Differential response function d**2 R / (d Eee d vmin) in [1/kg/keV/(km/s)]
            NOT including the velocity integral eta0
            Input:
                Eee: measured energy (electron equivalent)
                qER: q * ER for quenching factor q and recoil energy ER
                const_factor: factors entering the differential response that do not
                    depend on Eee
        '''
        self.count_diffresponse_calls += 1
        r_list = const_factor * self.Efficiency(Eee) * \
            np.array([self.ResolutionFunction(Eee, qer, self.EnergyResolution(qer))
                      for qer in qER])
        return r_list.sum()

    def ConstFactor(self, vmin, mx, fp, fn, delta, sign):
        ''' Collects the factors that don't depend on the measured energy Eee,
        so they only need to be computed once in Response function.
            Returns:
                a touple of ER, qER and const_factor
        '''
        ER = ERecoilBranch(vmin, self.mT, mx, delta, sign)
        q = self.QuenchingFactor(ER)
        qER = q * ER
        efficiencyER = self.Efficiency_ER(ER)
        const_factor = kilogram/SpeedOfLight**2 * \
            self.CrossSectionFactors(ER, mx, fp, fn, delta) * \
            np.abs(dERecoildVmin(vmin, self.mT, mx, delta, sign)) * efficiencyER
            # 1/kg/(km/s)
        return (ER, qER, const_factor)

    def DifferentialResponse_Full(self, vmin, Eee, mx, fp, fn, delta, sign):
        ''' Differential response function d**2 R / (d Eee d ER)
            NOT including the velocity integral eta0
            Same as DifferentialResponse, but computed given full input parameters,
        instead of the pre-computed const_factor.
        '''
        (ER, qER, const_factor) = self.ConstFactor(vmin, mx, fp, fn, delta, sign)
        return self.DifferentialResponse(Eee, qER, const_factor)

    def Response_Other(self, vmin, Eee1, Eee2, mx, fp, fn, delta):
        ''' Response function integral d**2 R / (d Eee d ER) between measured energies
        Eee1 and Eee2.
        NOT including eta0.
            For any finite resolution function (i.e. other than Dirac Delta).
        '''
        self.count_response_calls += 1
        if delta == 0:
            branches = [1]
        else:
            branches = [1, -1]
        result = 0
        for sign in branches:
            (ER, qER, const_factor) = self.ConstFactor(vmin, mx, fp, fn, delta, sign)
            result += integrate.quad(self.DifferentialResponse, Eee1, Eee2,
                                     args=(qER, const_factor),
                                     epsrel=PRECISSION, epsabs=0)[0]
        return result

    def Response_Dirac(self, vmin, Eee1, Eee2, mx, fp, fn, delta):
        ''' Response function integral d**2 R / (d Eee d ER) between measured energies
        Eee1 and Eee2,
        NOT including eta0.
            For Dirac Delta resolution function.
        '''
        self.count_response_calls += 1
        if delta == 0:
            branches = [1]
        else:
            branches = [1, -1]
        r_list_sum = 0
        for sign in branches:
            ER = ERecoilBranch(vmin, self.mT, mx, delta, sign)
            q = self.QuenchingFactor(ER)
            qER = q * ER
            integrated_delta = np.array([1. if Eee1 <= i < Eee2 else 0. for i in qER])
            efficiencyEee = self.Efficiency(Eee1, qER)
#            efficiencyER = self.Efficiency_ER(qER)
            efficiencyER = np.array(list(map(self.Efficiency_ER, qER)))
            r_list = kilogram/SpeedOfLight**2 * self.CrossSectionFactors(ER, mx, fp, fn, delta) * \
                np.abs(dERecoildVmin(vmin, self.mT, mx, delta, sign)) * \
                efficiencyEee * efficiencyER * integrated_delta
            r_list_sum += r_list.sum()
        return r_list_sum

    def IntegratedResponse_Other(self, vmin1, vmin2, Eee1, Eee2, mx, fp, fn, delta):
        ''' Integrated Response Function between measured energies Eee1 and Eee2,
        and all recoil energies ER.
        NOT including eta0.
            For any finite resolution function (i.e. other than Dirac Delta).
        '''
        midpoints = []
        integr = integrate.quad(self.Response_Other, vmin1, vmin2,
                                args=(Eee1, Eee2, mx, fp, fn, delta), points=midpoints,
                                epsrel=PRECISSION, epsabs=0)
        return integr[0]

    def IntegratedResponse_Dirac(self, vmin1, vmin2, Eee1, Eee2, mx, fp, fn, delta):
        ''' Integrated Response Function between measured energies Eee1 and Eee2,
        and all recoil energies ER.
        NOT including eta0.
            For Dirac Delta resolution function.
        '''
        # TODO! This is only valid for quenching factor 1!!! Extend to arbitrary q!
        E_delta = - delta * mx / (self.mT + mx)  # = muT * delta / self.mT for delta <= 0
        vmin_of_E1 = VMin(Eee1, self.mT, mx, delta)
        vmin_of_E2 = VMin(Eee2, self.mT, mx, delta)
        if Eee1 <= max(E_delta) and Eee2 >= min(E_delta):
            vmin_min = 0
        else:
            vmin_min = min(min(vmin_of_E1), min(vmin_of_E2))
        vmin_max = max(max(vmin_of_E1), max(vmin_of_E2))
        vmin1 = max(vmin_min, vmin1)
        vmin2 = min(vmin_max, vmin2)
        if vmin1 > vmin2:
            return 0

        integr = integrate.quad(self.Response_Dirac, vmin1, vmin2,
                                args=(Eee1, Eee2, mx, fp, fn, delta))  # , vec_func=False
#        print("Eee1, Eee2, integr = ", Eee1, " ", Eee2, " ", integr)
        return integr[0]


class MaxGapExperiment_HaloIndep(Experiment_HaloIndep):
    def __init__(self, expername, scattering_type, mPhi=mPhiRef):
        Experiment_HaloIndep.__init__(self, expername, scattering_type, mPhi)
        module = import_file(INPUT_DIR + expername + ".py")
        self.ERecoilList = module.ERecoilList
        self.ElistMaxGap = np.append(np.insert(
            np.array(list(filter(lambda x: self.Ethreshold < x < self.Emaximum,
                          self.ERecoilList))),
            0, self.Ethreshold), self.Emaximum)

    def TabulateMaximumGapLimit(self, vmin1, vmin2, mx, fp, fn, delta):
        print("vmin = ", vmin2)
        return np.array(list(map(lambda i, j:
                        self.IntegratedResponse(vmin1, vmin2, i, j, mx, fp, fn, delta),
                        self.ElistMaxGap[:-1], self.ElistMaxGap[1:])))

    def MaximumGapUpperBound(self, vmin_min, vmin_max, vmin_step, mx, fp, fn, delta,
                             output_file, processes=None):
        vmin_list = np.linspace(vmin_min, vmin_max, (vmin_max - vmin_min)/vmin_step + 1)
        vmin_list0 = np.insert(vmin_list, 0, 0.)
        xtable = np.zeros(self.ElistMaxGap.size - 1)
        upperlimit_table = np.array([])
        kwargs = ({'vmin1': vmin_list0[v_index], 'vmin2': vmin_list[v_index],
                   'mx': mx, 'fp': fp, 'fn': fn, 'delta': delta}
                  for v_index in range(vmin_list.size))
        xtable_list = par.parmap(self.TabulateMaximumGapLimit, kwargs, processes)
        for v_index in range(vmin_list.size):
            xtable += xtable_list[v_index]
            mu_scaled = xtable.sum()
            x_scaled = np.max(xtable)
            if x_scaled == 0:
                mu_over_x = np.inf
                result = [np.inf]
            else:
                mu_over_x = mu_scaled / x_scaled
                y_guess = np.real(-lambertw(-0.1 / mu_over_x, -1))
                y = fsolve(lambda x: MaximumGapC0scaled(x, mu_over_x) - ConfidenceLevel,
                           y_guess)
                result = y / x_scaled / self.Exposure
                result = result[0]
                print("vmin = ", vmin_list[v_index], "   mu_over_x = ", mu_over_x)
                print("xtable = ", xtable)
                print("result = ", result)
                to_print = np.log10(np.array([[mx, result]]))
                with open(output_file, 'ab') as f_handle:
                    np.savetxt(f_handle, to_print)
            upperlimit_table = np.append(upperlimit_table, [result])
        return upperlimit_table

    def UpperLimit(self, mx, fp, fn, delta, vmin_min, vmin_max, vmin_step,
                   output_file, processes=None, **unused_kwargs):
        upper_limit = self.MaximumGapUpperBound(vmin_min, vmin_max, vmin_step, mx,
                                                fp, fn, delta, output_file,
                                                processes=processes)
        vmin_list = np.linspace(vmin_min, vmin_max, (vmin_max - vmin_min)/vmin_step + 1)
        print("vmin_list = ", vmin_list)
        print("upper_limit = ", upper_limit)
        result = np.transpose([vmin_list, np.log10(upper_limit)])
        print("res = ", result)
        return result[result[:, 1] != np.inf]


class PoissonExperiment_HaloIndep(Experiment_HaloIndep):
    def __init__(self, expername, scattering_type, mPhi=mPhiRef):
        Experiment_HaloIndep.__init__(self, expername, scattering_type, mPhi)
        module = import_file(INPUT_DIR + expername + ".py")
        self.Expected_limit = module.Expected_limit

    def _PoissonUpperBound(self, vmin, mx, fp, fn, delta):
        muT = self.mT * mx / (self.mT + mx)
        Eee_max = max(2e6 * muT**2 * (vmin/SpeedOfLight)**2 / self.mT)
        print("self.Ethreshold =", self.Ethreshold)
        print("Eee_max =", Eee_max)
        int_response = self.IntegratedResponse(0, vmin, self.Ethreshold, Eee_max,
                                               mx, fp, fn, delta)
        print("int_response =", int_response)
        if int_response > 0:
            result = np.log10(self.Expected_limit / self.Exposure / int_response)
        else:
            result = np.inf
        return [vmin, result]

    def UpperLimit(self, mx, fp, fn, delta, vmin_min, vmin_max, vmin_step,
                   output_file, processes=None, **unused_kwargs):
        vmin_list = np.linspace(vmin_min, vmin_max, (vmin_max - vmin_min)/vmin_step + 1)
        kwargs = ({'vmin': vmin, 'mx': mx, 'fp': fp, 'fn': fn, 'delta': delta}
                  for vmin in vmin_list)
        upper_limit = np.array(par.parmap(self._PoissonUpperBound, kwargs, processes))
        upper_limit = upper_limit[upper_limit[:, 1] != np.inf]
        print("upper_limit = ", upper_limit)
        with open(output_file, 'ab') as f_handle:
            np.savetxt(f_handle, upper_limit)
        return upper_limit


class GaussianExperiment_HaloIndep(Experiment_HaloIndep):
    def __init__(self, expername, scattering_type, mPhi=mPhiRef, quenching_factor=None):
        Experiment_HaloIndep.__init__(self, expername, scattering_type, mPhi)
        module = import_file(INPUT_DIR + expername + ".py")
        self.BinEdges_left = module.BinEdges_left
        self.BinEdges_right = module.BinEdges_right
        self.BinData = module.BinData
        self.BinError = module.BinError
        self.BinSize = module.BinSize
        self.chiSquared = module.chiSquared[1]
#        self.Expected_limit = (np.sqrt(self.chiSquared) * self.BinError + self.BinData) * \
#            self.BinSize
        self.Expected_limit = module.Expected_limit * self.BinSize # 1/day/kg
        if quenching_factor is not None:
            self.QuenchingFactor = lambda e: quenching_factor

    def _GaussianUpperBound(self, vmin, mx, fp, fn, delta):
        int_response = \
            np.array(list(map(lambda i, j:
                              self.IntegratedResponse(0, vmin, i, j, mx, fp, fn, delta),
                              self.BinEdges_left, self.BinEdges_right)))
        result = np.min(self.Expected_limit / self.Exposure / int_response)
        if result > 0:
            result = np.log10(result)
        else:
            result = np.inf
        print("(vmin, result) =", (vmin, result))
        return [vmin, result]

    def UpperLimit(self, mx, fp, fn, delta, vmin_min, vmin_max, vmin_step,
                   output_file, processes=None, **unused_kwargs):
        vmin_list = np.linspace(vmin_min, vmin_max, (vmin_max - vmin_min)/vmin_step + 1)
        kwargs = ({'vmin': vmin, 'mx': mx, 'fp': fp, 'fn': fn, 'delta': delta}
                  for vmin in vmin_list)
        upper_limit = np.array(par.parmap(self._GaussianUpperBound, kwargs, processes))
        upper_limit = upper_limit[upper_limit[:, 1] != np.inf]
        print("upper_limit = ", upper_limit)
        with open(output_file, 'ab') as f_handle:
            np.savetxt(f_handle, upper_limit)
        return upper_limit


class Crosses_HaloIndep(Experiment_HaloIndep):
    def __init__(self, expername, scattering_type, mPhi=mPhiRef, quenching_factor=None):
        Experiment_HaloIndep.__init__(self, expername, scattering_type, mPhi)
        module = import_file(INPUT_DIR + expername + ".py")
        self.BinEdges = module.BinEdges
        self.BinData = module.BinData
        self.BinError = module.BinError
        self.QuenchingFactorOfEee = module.QuenchingFactorOfEee
        if quenching_factor is not None:
            self.QuenchingFactor = lambda e: quenching_factor

    def _VminRange(self, E1, E2, mT, mx, delta):
        E_delta = - delta * mx / (mT + mx)
        vmin_of_E1 = VMin(E1, mT, mx, delta)
        vmin_of_E2 = VMin(E2, mT, mx, delta)
        print(vmin_of_E1, vmin_of_E2)
        if E1 <= E_delta and E2 >= E_delta:
            vmin_min = 0
        else:
            vmin_min = min(vmin_of_E1, vmin_of_E2)
        vmin_max = max(vmin_of_E1, vmin_of_E2)
        return (vmin_min, vmin_max)

    def _AverageOverNuclides(self, quantity):
        return np.sum(quantity * self.mass_fraction) / np.sum(self.mass_fraction)

    def _Box(self, Eee1, Eee2, mT_avg, mx, fp, fn, delta, nsigma):
        E1 = Eee1 - nsigma * self.EnergyResolution(Eee1)
        E2 = Eee2 + nsigma * self.EnergyResolution(Eee2)
        ER1 = self._AverageOverNuclides(E1 / self.QuenchingFactorOfEee(E1))
        ER2 = self._AverageOverNuclides(E2 / self.QuenchingFactorOfEee(E2))
        (vmin1, vmin2) = self._VminRange(ER1, ER2, mT_avg, mx, delta)
        int_resp = self.IntegratedResponse(vmin1, vmin2, Eee1, Eee2, mx, fp, fn, delta)
        vmin_center = (vmin1 + vmin2)/2
        vmin_error = (vmin2 - vmin1)/2
        return (int_resp, vmin_center, vmin_error)

    def _Boxes(self, mx, fp, fn, delta, nsigma=1, processes=None):
        mT_avg = np.sum(self.mT * self.mass_fraction) / np.sum(self.mass_fraction)
        print("mT_avg =", mT_avg)
        kwargs = ({'Eee1': Eee1, 'Eee2': Eee2, 'mT_avg': mT_avg,
                   'mx': mx, 'fp': fp, 'fn': fn, 'delta': delta,
                   'nsigma': nsigma}
                  for Eee1, Eee2 in np.transpose([self.BinEdges[:-1], self.BinEdges[1:]]))
        return np.array(par.parmap(self._Box, kwargs, processes))

    def UpperLimit(self, mx, fp, fn, delta, vmin_min, vmin_max, vmin_step,
                   output_file, nsigma=1, processes=None):
        box_table = self._Boxes(mx, fp, fn, delta, nsigma=nsigma, processes=1)
        int_resp_list = box_table[:, 0]
        vmin_center_list = box_table[:, 1]
        vmin_error_list = box_table[:, 2]
        eta_list = self.BinData / int_resp_list
        eta_error_list = self.BinError / int_resp_list
        result = np.array([int_resp_list, vmin_center_list, vmin_error_list,
                           eta_list, eta_error_list])
        print(result)
        with open(output_file, 'ab') as f_handle:
            np.savetxt(f_handle, result)
        return result
