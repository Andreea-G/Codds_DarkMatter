# -*- coding: utf-8 -*-
"""
Created on Sun Mar  1 19:51:42 2015

@author: Andreea
"""

from experiment import *

class Experiment_HaloIndep(Experiment):
    def __init__(self, expername, scattering_type, mPhi = mPhiRef):
        Experiment.__init__(self, expername, scattering_type, mPhi)
#        module = import_file(INPUT_DIR + expername + ".py")
        if self.energy_resolution_type == "Dirac":
            self.IntegratedResponse = self.IntegratedResponse_Dirac
        else:
            self.IntegratedResponse = self.IntegratedResponse_Other

    def DifferentialResponse(self, Eee, qER, const_factor): 
        ''' Differential response function d**2 R / (d Eee d ER), as a function of measured energy Eee and recoil energy ER,
         NOT including eta0
        '''
        self.count_diffresponse_calls += 1
        r_list = const_factor * self.Efficiency(Eee) * self.ResolutionFunction(Eee, qER, self.EnergyResolution(qER))
        return r_list.sum()

    def ConstFactor(self, vmin, mx, fp, fn, delta, sign):
        ''' Collects the factors that don't depend on Eee, so they only need to be computed once in Response fnc
            Returns a touple of qER and const_factor
        '''
        ER = ERecoilBranch(vmin, self.mT, mx, delta, sign)
        q = self.QuenchingFactor(ER)
        qER = q * ER
        efficiencyER = self.Efficiency_ER(ER)
        const_factor = kilogram/SpeedOfLight**2 * self.CrossSectionFactors(ER, mx, fp, fn, delta) * \
            dERecoildVmin(vmin, self.mT, mx, delta, sign) * efficiencyER
        return (qER, const_factor)

    def Response_Other(self, vmin, Eee1, Eee2, mx, fp, fn, delta):
        ''' Response function integral d**2 R / (d Eee d ER) between measured energies Eee1 and Eee2, 
        as a function of recoil energy ER, NOT including eta0.
            For any finite resolution function (i.e. other than Dirac Delta).
        '''
        self.count_response_calls += 1
        if delta == 0:
            branches = [1]
        else:
            branches = [1, -1]
        result = 0
        for sign in branches:
            (qER, const_factor) = _ConstFactor(self, vmin, mx, fp, fn, delta, sign)
            result += integrate.quad(self.DifferentialResponse, Eee1, Eee2, \
                args=(qER, const_factor), epsrel = PRECISSION, epsabs = 0)[0]
        return result

    def Response_Dirac(self, vmin, Eee1, Eee2, mx, fp, fn, delta): 
        ''' Response function integral d**2 R / (d Eee d ER) between measured energies Eee1 and Eee2, 
        as a function of recoil energy ER, NOT including eta0.
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
#            print("ER, q, qER, intdelta, effEee, effER = ", ER, "; ", q, "; ", qER, "; ", integrated_delta, "; ", efficiencyEee, "; ", efficiencyER)
#            print("CSfactors = ", self.CrossSectionFactors(ER, mx, fp, fn, delta))
#            print("deriv = ", dERecoildVmin(vmin, self.mT, mx, delta, sign))
            r_list = kilogram/SpeedOfLight**2 * self.CrossSectionFactors(ER, mx, fp, fn, delta) * \
                dERecoildVmin(vmin, self.mT, mx, delta, sign) * \
                efficiencyEee * efficiencyER * integrated_delta
            r_list_sum += r_list.sum()
        return r_list_sum

    def IntegratedResponse_Other(self, vmin1, vmin2, Eee1, Eee2, mx, fp, fn, delta):
        ''' Integrated Response Function between measured energies Eee1 and Eee2, and all of recoil energies ER,
        NOT including eta0.
            For any finite resolution function (i.e. other than Dirac Delta).
        '''
        midpoints = []
        integr = integrate.quad(self.Response_Other, vmin1, vmin2, \
            args=(Eee1, Eee2, mx, fp, fn, delta), points = midpoints, epsrel = PRECISSION, epsabs = 0)
        print("Eee1, Eee2, integr = ", Eee1, " ", Eee2, " ", integr)
        return integr[0]

    def IntegratedResponse_Dirac(self, vmin1, vmin2, Eee1, Eee2, mx, fp, fn, delta):
        ''' Integrated Response Function between measured energies Eee1 and Eee2, and all of recoil energies ER,
        NOT including eta0.
            For Dirac Delta resolution function.
        '''
        #TODO! This is only valid for quenching factor 1!!! Extend to arbitrary q!
        E_delta = - delta * mx / (self.mT + mx)  #  = muT * delta / self.mT for delta <= 0
        vmin_of_E1 = VMin(Eee1, self.mT, mx, delta)
        vmin_of_E2 = VMin(Eee2, self.mT, mx, delta)
        if Eee1 <= max(E_delta) or Eee2 >= min(E_delta):
            vmin_min = 0
        else:
            vmin_min = min(min(vmin_of_E1), min(vmin_of_E2))
        vmin_max = max(max(vmin_of_E1), max(vmin_of_E2))
        vmin1 = max(vmin_min, vmin1)
        vmin2 = min(vmin_max, vmin2)
        if vmin1 > vmin2:
            return 0
        
        integr = integrate.quad(self.Response_Dirac, vmin1, vmin2, \
            args=(Eee1, Eee2, mx, fp, fn, delta)) #, vec_func=False
        print("Eee1, Eee2, integr = ", Eee1, " ", Eee2, " ", integr)
        return integr[0]


class MaxGapExperiment_HaloIndep(Experiment_HaloIndep):
    def __init__(self, expername, scattering_type, mPhi = mPhiRef):
        Experiment_HaloIndep.__init__(self, expername, scattering_type, mPhi)
        module = import_file(INPUT_DIR + expername + ".py")
        self.ERecoilList = module.ERecoilList
        self.ElistMaxGap = np.append( np.insert( \
            np.array(list(filter(lambda x: self.Ethreshold < x < self.Emaximum, \
            self.ERecoilList))), \
            0, self.Ethreshold), self.Emaximum)
        print("Emin, Emax = ", self.Ethreshold, " ", self.Emaximum)
        print("elist = ", self.ElistMaxGap)
        

    def TabulateMaximumGapLimit(self, vmin_min, vmin_max, vmin_step, mx, fp, fn, delta, output_file):
        vmin_list = np.linspace(vmin_min, vmin_max, (vmin_max - vmin_min)/vmin_step + 1)
        vmin_list0 = np.insert(vmin_list,0,0.)
        xtable = np.zeros(self.ElistMaxGap.size - 1)
        upperlimit_table = np.array([])
        for v_index in range(vmin_list.size):
            print("vmin = ", vmin_list[v_index])
            xtable += np.array(list(map(lambda i, j: \
                self.IntegratedResponse(vmin_list0[v_index], vmin_list[v_index], i, j, mx, fp, fn, delta), \
                self.ElistMaxGap[:-1], self.ElistMaxGap[1:])))
            mu_scaled = xtable.sum()
            x_scaled = np.max(xtable)
            if x_scaled == 0:
                mu_over_x = np.inf
                result = [np.inf]
            else:
                mu_over_x = mu_scaled / x_scaled
                y_guess = np.real(-lambertw(-0.1 / mu_over_x, -1))
                y = fsolve(lambda x: MaximumGapC0scaled(x, mu_over_x) - ConfidenceLevel, y_guess)
                result =  y / x_scaled / self.Exposure
                result = result[0]
                print("vmin = ", vmin_list[v_index], "   mu_over_x = ", mu_over_x)
                print("xtable = ", xtable)
                print("result = ", result)
                to_print = np.log10(np.array([[mx, result]]))
                with open(output_file,'ab') as f_handle:
                    np.savetxt(f_handle, to_print)
            upperlimit_table = np.append(upperlimit_table, [result])
        return upperlimit_table
        
    def UpperLimit(self, mx, fp, fn, delta, vmin_min, vmin_max, vmin_step, output_file):
        upper_limit = self.TabulateMaximumGapLimit(vmin_min, vmin_max, vmin_step, mx, fp, fn, delta, output_file)
        vmin_list = np.linspace(vmin_min, vmin_max, (vmin_max - vmin_min)/vmin_step + 1)
        print("vmin_list = ", vmin_list)
        print("upper_limit = ", upper_limit)
        result = np.transpose([vmin_list, np.log10(upper_limit)])
        print("res = ", result)
        return result[result[:, 1] != np.inf]
        
        
        upper_limit = np.array([self.TabulateMaximumGapLimit(vmin_list0[i], vmin_list[i], mx, fp, fn, delta, output_file) \
            for i in range(vmin_list.size)]).flatten()
