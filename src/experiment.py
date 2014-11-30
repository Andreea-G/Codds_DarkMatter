# -*- coding: utf-8 -*-
"""
Created on Thu Nov 20 22:52:11 2014

@author: Andreea
"""
from __future__ import absolute_import
from __future__ import print_function
from __future__ import division


INPUT_DIR = "Data/"
OUTPUT_MAIN_DIR = "Output/"
PRECISSION = 1e-3

def import_file(full_path_to_module):
    import os, sys
    directory, module_name = os.path.split(full_path_to_module)
    module_name = os.path.splitext(module_name)[0]
        
    path = list(sys.path)
    sys.path.insert(0, directory)
    try:
        module = __import__(module_name)
        return module
    finally:
        sys.path[:] = path # restore 
        
        
from globalfnc import *
from haxtonFF import *
import numpy as np
from numpy import pi
from scipy import integrate
from scipy.optimize import fsolve
from scipy.special import lambertw

class Experiment:
    def __init__(self, expername, scattering_type, mPhi = mPhiRef):
        module = import_file(INPUT_DIR + expername + ".py")
        self.name = expername
        self.scattering_type = scattering_type
        self.energy_resolution_type = module.energy_resolution_type
        self.EnergyResolution = module.EnergyResolution
        self.ResolutionFunction = Gaussian if \
            self.energy_resolution_type == "Gaussian" else GPoisson  
        
        self.mPhi = mPhi
        self.numT = module.num_target_nuclides
        self.mT = module.target_nuclide_mass_list
        self.A = module.target_nuclide_AZC_list[:,0]
        self.Z = module.target_nuclide_AZC_list[:,1]
        self.mass_fraction = module.target_nuclide_AZC_list[:,2]
        self.J = module.target_nuclide_JSpSn_list[:,0]
        self.SpScaled = module.target_nuclide_JSpSn_list[:,1]
        self.SnScaled = module.target_nuclide_JSpSn_list[:,2]

        self.bsq = 41.467/(45. * self.A**(-1./3) - 25. * self.A**(-2./3)) * fermiGeV**2
        self.FF66_function_list = np.array(map(lambda a, z: \
            FFSigmaPPJ.get((np.trunc(a), np.trunc(z)), np.array([[lambda y: 0]*2]*2)), \
            self.A, self.Z))

        self.FF = FF_options[self.scattering_type][module.FF[scattering_type]]

        CrossSectionFactors_options = {'SI' : self.CrossSectionFactors_SI,
           'SD66' : self.CrossSectionFactors_SD66,
           'SD44' : self.CrossSectionFactors_SD44,
        }
        self.CrossSectionFactors = CrossSectionFactors_options[self.scattering_type]
        if self.energy_resolution_type == "Dirac":
            self.IntegratedResponseSHM = self.IntegratedResponseSHM_Dirac
        else:
            self.IntegratedResponseSHM = self.IntegratedResponseSHM_Other
            
        self.QuenchingFactor = module.QuenchingFactor  
        self.Efficiency = module.Efficiency
        self.Efficiency_ER = module.Efficiency_ER
        self.Ethreshold = module.Ethreshold
        self.Emaximum = module.Emaximum
        self.ERmaximum = module.ERmaximum
        self.ERecoilList = module.ERecoilList
        self.ElistMaxGap = np.append( np.insert( \
            np.array(filter(lambda x: self.Ethreshold < x < self.Emaximum, self.ERecoilList)), \
            0, self.Ethreshold), self.Emaximum)
        self.Exposure = module.Exposure
        
        self.count_diffresponse_calls = 0
        self.count_response_calls = 0
        
    def FormFactor(self, ER):
        result = self.FF(ER, self.A, self.mT)
        return result
    
    def CrossSectionFactors_SI(self, ER, mx, fp, fn, delta):
        mu_p = ProtonMass * mx / (ProtonMass + mx)
        return self.mass_fraction * 1./(2.*mu_p**2) * \
            mPhiRef**4 / (4. * self.mT**2 * (ER + self.mPhi**2/(2. * self.mT))**2) * \
            ((self.Z + (self.A - self.Z) * fn/fp)**2) * self.FormFactor(ER) 

    def FF66normlalized(self, ER):
        y = 2.e-6 * self.mT * ER * self.bsq / 4.
        return np.array([np.array([[self.FF66_function_list[i,N1,N2](y[i]) \
            for N2 in [0,1]] for N1 in [0,1]]) \
            * np.exp(-2. * y[i]) for i in range(self.numT)])
#        l = np.array([self.FF66_function_list[:, N1, N2][i](y[i]) \
#            for i in range(self.numT)])
#        return l * np.exp(-2. * y)

    '''
    def FF66normlalized(self, ER, N1, N2):
        y = 2.e-6 * self.mT * ER * self.bsq / 4.
        l = np.array([self.FF66_function_list[:, N1, N2][i](y[i]) \
            for i in range(self.numT)])
        return l * np.exp(-2. * y)
    '''
    
    def CrossSectionFactors_SD66(self, ER, mx, fp, fn, delta):
        #ffelemQ = FFElementQ(self.Z)
        mu_p = ProtonMass * mx / (ProtonMass + mx)
        ff = self.FF66normlalized(ER)
        return self.mass_fraction * 3./(8.*mu_p**6) * self.mT**2 * 1e-12 * ER**2 * \
            (SpeedOfLight/v0bar)**4 * \
            mPhiRef**4 / (4. * self.mT**2 * (ER + self.mPhi**2/(2 * self.mT))**2) * \
            (ff[:,0,0] + ff[:,0,1] * 2 * fn/fp + ff[:,1,1] * (fn/fp)**2) 
        '''        
        return self.mass_fraction * 3./(8.*mu_p**6) * self.mT**2 * 1e-12 * ER**2 * \
            (SpeedOfLight/v0bar)**4 * \
            mPhiRef**4 / (4. * self.mT**2 * (ER + self.mPhi**2/(2 * self.mT))**2) * \
            (self.FF66normlalized(ER, 0, 0) + \
            self.FF66normlalized(ER, 0, 1) * 2 * fn/fp + \
            self.FF66normlalized(ER, 1, 1) * (fn/fp)**2) 
        return self.mass_fraction * 3./(8.*mu_p**6) * self.mT**2 * 1e-12 * ER**2 * \
            (SpeedOfLight/v0bar)**4 * \
            mPhiRef**4 / (4. * self.mT**2 * (ER + self.mPhi**2/(2 * self.mT))**2) * \
            (ffelemQ * (FF66normlalized(ER, self.A, self.Z, self.mT, 0, 0) + \
            FF66normlalized(ER, self.A, self.Z, self.mT, 0, 1) * 2 * fn/fp + \
            FF66normlalized(ER, self.A, self.Z, self.mT, 1, 1) * (fn/fp)**2) + \
            (1 - ffelemQ) * (4./3. * (4. * pi)/(2 * self.J + 1.)) * \
            (self.SpScaled + self.SnScaled * fn/fp)**2 * self.FormFactor(ER))
        '''
        
    #TODO
    def CrossSectionFactors_SD44(self, ER, mx, fp, fn, delta):
        #print(exper.name, "SD44")
        return 0        
        
    def Resolution(self, Eee, qER):
        return self.ResolutionFunction(Eee, qER, self.EnergyResolution(qER))
        
    def DifferentialResponseSHM(self, Eee, ER, mx, fp, fn, delta): 
        self.count_diffresponse_calls += 1
        q = self.QuenchingFactor(ER)
        qER = q * ER
        vmin = VMin(ER, self.mT, mx, delta)
        r_list = 1.e-6 * kilogram * self.CrossSectionFactors(ER, mx, fp, fn, delta) * \
            self.Efficiency(Eee) * self.Efficiency_ER(ER) * \
            self.ResolutionFunction(Eee, qER, self.EnergyResolution(qER)) * \
            eta0Maxwellian(vmin, vobs, v0bar, vesc)
#        print("ER, Eee, resp = ", ER, " ", Eee, " ", r_list.sum())
        return r_list.sum()
        
    def ResponseSHM_Dirac(self, ER, Eee1, Eee2, mx, fp, fn, delta): 
        q = self.QuenchingFactor(ER)
        qER = q * ER
        vmin = VMin(ER, self.mT, mx, delta)
        efficiency_ER = self.Efficiency_ER(qER)
        integrated_delta = 1. if Eee1 <= qER < Eee2 else 0.
        r_list = 1.e-6 * kilogram * self.CrossSectionFactors(ER, mx, fp, fn, delta) * \
            efficiency_ER * \
            integrated_delta * eta0Maxwellian(vmin, vobs, v0bar, vesc)
        self.count_response_calls += 1
        return r_list.sum()
#       print(ER, " ", r)
#        return r
        
    def ResponseSHM_Other(self, ER, Eee1, Eee2, mx, fp, fn, delta):
        self.count_response_calls += 1
        result = integrate.quad(self.DifferentialResponseSHM, Eee1, Eee2, \
            args=(ER, mx, fp, fn, delta), epsrel = PRECISSION, epsabs = 0)[0]
#        print(ER, " ", result)
        return result

    def IntegratedResponseSHM_Dirac(self, Eee1, Eee2, mx, fp, fn, delta):
        vmax = vesc + vobs
        muT = self.mT * mx / (self.mT + mx)
        vdelta = SpeedOfLight / 500. * np.sqrt(delta / 2. / muT) if delta > 0 \
            else np.array([0] * self.numT)
        ER_plus_list = map(lambda i, j: ERecoilBranch(vmax, i, mx, delta, 1) \
            if j < vmax else 0., self.mT, vdelta)
        ER_minus_list = map(lambda i, j: ERecoilBranch(vmax, i, mx, delta, -1) \
            if j < vmax else 1.e6, self.mT, vdelta)
        ER_plus = min(np.max(ER_plus_list), self.ERmaximum)
        ER_minus = np.min(ER_minus_list)
        if ER_minus < ER_plus:
            integr = integrate.quad(self.ResponseSHM_Dirac, ER_minus, ER_plus, \
                args=(Eee1, Eee2, mx, fp, fn, delta)) #, vec_func=False
            print("Eee1, Eee2, integr = ", Eee1, " ", Eee2, " ", integr)
            return integr[0]
        else:
            return 0.
            
    def IntegratedResponseSHM_Other(self, Eee1, Eee2, mx, fp, fn, delta):
        vmax = vesc + vobs
        muT = self.mT * mx / (self.mT + mx)
        vdelta = SpeedOfLight / 500. * np.sqrt(delta / 2. / muT) if delta > 0 \
            else np.array([0] * self.numT)
        ER_plus_list = map(lambda i, j: ERecoilBranch(vmax, i, mx, delta, 1) \
            if j < vmax else 0., self.mT, vdelta)
        ER_minus_list = map(lambda i, j: ERecoilBranch(vmax, i, mx, delta, -1) \
            if j < vmax else 1.e6, self.mT, vdelta)
        ER_plus = min(np.max(ER_plus_list), self.ERmaximum)
        ER_minus = np.min(ER_minus_list)
        if ER_minus < ER_plus:
            integr = integrate.quad(self.ResponseSHM_Other, ER_minus, ER_plus, \
                args=(Eee1, Eee2, mx, fp, fn, delta), epsrel = PRECISSION, epsabs = 0)
            '''
            integr = integrate.dblquad(self.DifferentialResponseSHM, ER_minus, ER_plus, \
                lambda Eee: Eee1, lambda Eee: Eee2, \
                args=(mx, fp, fn, delta), epsrel = PRECISSION, epsabs = 0)
            '''
            print("Eee1, Eee2, integr = ", Eee1, " ", Eee2, " ", integr,)
            return integr[0]
        else:
            return 0.
    
            
    def MaximumGapUpperBoundSHM(self, mx, fp, fn, delta, output_file):
        print("mx = ", mx)
        xtable = np.array(map(lambda i, j: \
            self.IntegratedResponseSHM(i, j, mx, fp, fn, delta), \
            self.ElistMaxGap[:-1], self.ElistMaxGap[1:]))
        mu_scaled = xtable.sum()
        x_scaled = np.max(xtable)
        if x_scaled == 0:
            mu_over_x = np.inf
            result = [np.inf]
        else:
            mu_over_x = mu_scaled / x_scaled
            y_guess = np.real(-lambertw(-0.1 / mu_over_x, -1))
#            print("y_guess = ", y_guess)
            y = fsolve(lambda x: MaximumGapC0scaled(x, mu_over_x) - 0.9, y_guess)
            result =  y / x_scaled / self.Exposure
        print("mx = ", mx, "   mu_over_x = ", mu_over_x)
        print("xtable = ", xtable)
        print("result = ", result[0])
        to_print = np.log10(np.array([[mx, result[0]]]))
        with open(output_file,'a') as f_handle:
            np.savetxt(f_handle, to_print)
        return result
        
    def MaximumGapLimit(self, fp, fn, delta, mx_min, mx_max, num_steps, output_file):
        mx_list = np.logspace(np.log10(mx_min), np.log10(mx_max), num_steps)
        upper_limit = np.array(map(lambda mx: \
            self.MaximumGapUpperBoundSHM(mx, fp, fn, delta, output_file), mx_list))
        print(mx_list)
        print(upper_limit)
        return np.log10(np.transpose([mx_list, upper_limit.flatten()]))
        
