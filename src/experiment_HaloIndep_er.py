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


class Experiment_HaloIndep_ER(Experiment):
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
            self.ResolutionFunction(Eee, qER, self.EnergyResolution(qER))
        return r_list.sum()

    def ConstFactor(self, ER, ER1, ER2, mx, fp, fn, delta, sign):
        ''' Collects the factors that don't depend on the measured energy Eee,
        so they only need to be computed once in Response function.
            Returns:
                a touple of ER, qER and const_factor
        '''
        qER = ER * self.QuenchingFactor(ER)
        efficiencyER = self.Efficiency_ER(ER)
        inrange = np.array([1 if (er1 - ER) * (ER-er2) >= 0 else 0
                            for er1, er2 in zip(ER1, ER2)])
        const_factor = 1e-6 * kilogram * \
            self.CrossSectionFactors(ER, mx, fp, fn, delta) * efficiencyER * inrange
        return (qER, const_factor)

    def DifferentialResponse_Full(self, ER, Eee, ER1, ER2, mx, fp, fn, delta, sign):
        ''' Differential response function d**2 R / (d Eee d ER)
            NOT including the velocity integral eta0
            Same as DifferentialResponse, but computed given full input parameters,
        instead of the pre-computed const_factor.
        '''
        (qER, const_factor) = self.ConstFactor(ER, ER1, ER2, mx, fp, fn, delta, sign)
        return self.DifferentialResponse(Eee, qER, const_factor)

    def Response_Other(self, ER, Eee1, Eee2, ER1, ER2, sign, mx, fp, fn, delta):
        ''' Response function integral d**2 R / (d Eee d ER) between measured energies
        Eee1 and Eee2.
        NOT including eta0.
            For any finite resolution function (i.e. other than Dirac Delta).
        '''
        self.count_response_calls += 1
        (qER, const_factor) = self.ConstFactor(ER, ER1, ER2,
                                               mx, fp, fn, delta, sign)
        result = integrate.quad(self.DifferentialResponse, Eee1, Eee2,
                                args=(qER, const_factor),
                                epsrel=PRECISSION, epsabs=0)[0]
        return result

    def Response_Dirac(self, ER, Eee1, Eee2, ER1, ER2, sign, mx, fp, fn, delta):
        ''' Response function integral d**2 R / (d Eee d ER) between measured energies
        Eee1 and Eee2,
        NOT including eta0.
            For Dirac Delta resolution function.
        '''
        self.count_response_calls += 1
        qER = ER * self.QuenchingFactor(ER)
        integrated_delta = 1. if Eee1 <= qER < Eee2 else 0.
        efficiencyEee = self.Efficiency(Eee1, qER)
        efficiencyER = self.Efficiency_ER(qER)

        inrange = np.array([1 if (er1 - ER) * (ER-er2) >= 0 else 0
                            for er1, er2 in zip(ER1, ER2)])

        r_list = 1e-6 * kilogram * self.CrossSectionFactors(ER, mx, fp, fn, delta) * \
            efficiencyEee * efficiencyER * integrated_delta * inrange
        r_list_sum = r_list.sum()
        return r_list_sum

    def IntegratedResponse_Dirac(self, vmin1, vmin2, Eee1, Eee2, mx, fp, fn, delta):
        ''' Integrated Response Function between measured energies Eee1 and Eee2,
        and all recoil energies ER.
        NOT including eta0.
        '''
        vdelta = VminDelta(self.mT, mx, delta)
        if delta > 0:
            if vmin2 < min(vdelta):
                return 0
            vmin1 = max(min(vdelta), vmin1)
        if delta == 0:
            branches = [1]
        else:
            branches = [1, -1]

        result = 0
        for sign in branches:
            ER1 = [ERecoilBranch(vmin1, mT, mx, delta, sign) if vmin1 >= vd else 0
                   for mT, vd in zip(self.mT, vdelta)]
            ER2 = [ERecoilBranch(vmin2, mT, mx, delta, sign) if vmin2 >= vd else 0
                   for mT, vd in zip(self.mT, vdelta)]
            ERmin = min(min(ER1), min(ER2))
            ERmax = max(max(ER1), max(ER2))
            # TODO! This is only valid for quenching factor 1!!! Extend to arbitrary q!
            ERmin = max(ERmin, Eee1)
            ERmax = min(ERmax, Eee2)
            if ERmin < ERmax:
                integr = integrate.quad(self.Response_Dirac, ERmin, ERmax,
                                        args=(Eee1, Eee2, ER1, ER2, sign,
                                              mx, fp, fn, delta))
                result += integr[0]
        return result

    def IntegratedResponse_Other(self, vmin1, vmin2, Eee1, Eee2, mx, fp, fn, delta):
        ''' Integrated Response Function between measured energies Eee1 and Eee2,
        and all recoil energies ER.
        NOT including eta0.
        '''
        vdelta = VminDelta(self.mT, mx, delta)
        if delta > 0:
            if vmin2 < min(vdelta):
                return 0
            vmin1 = max(min(vdelta), vmin1)
        if delta == 0:
            branches = [1]
        else:
            branches = [1, -1]

        result = 0
        for sign in branches:
            ER1 = [ERecoilBranch(vmin1, mT, mx, delta, sign) if vmin1 >= vd else 0
                   for mT, vd in zip(self.mT, vdelta)]
            ER2 = [ERecoilBranch(vmin2, mT, mx, delta, sign) if vmin2 >= vd else 0
                   for mT, vd in zip(self.mT, vdelta)]
            ERmin = min(min(ER1), min(ER2))
            ERmax = max(max(ER1), max(ER2))
            midpoints = []
            if ERmin < Eee1 < ERmax:
                midpoints += [Eee1]
            if ERmin < Eee2 < ERmax:
                midpoints += [Eee2]
            if delta > 0:
                E_delta = delta * mx / (self.mT + mx)
                midpoints += [Ed for Ed in E_delta if ERmin < Ed < ERmax]
            if ERmin < ERmax:
                midpoints = sorted(midpoints)
                print('midpoints =', midpoints)
                integr = integrate.quad(self.Response_Other, ERmin, ERmax,
                                        args=(Eee1, Eee2, ER1, ER2, sign,
                                              mx, fp, fn, delta),
                                        points=midpoints, epsrel=PRECISSION, epsabs=0)
                result += integr[0]
        print('result =', result)
        return result
