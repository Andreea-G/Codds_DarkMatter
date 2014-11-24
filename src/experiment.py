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
        
        self.mPhi = mPhi
        self.numT = module.num_target_nuclides
        self.mT = module.target_nuclide_mass_list
        self.A = module.target_nuclide_AZC_list[:,0]
        self.Z = module.target_nuclide_AZC_list[:,1]
        self.mass_fraction = module.target_nuclide_AZC_list[:,2]
        self.J = module.target_nuclide_JSpSn_list[:,0]
        self.SpScaled = module.target_nuclide_JSpSn_list[:,1]
        self.SnScaled = module.target_nuclide_JSpSn_list[:,2]

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
            
        self.QuenchingFactorList = module.QuenchingFactorList  
        self.Efficiency = module.Efficiency
        self.Ethreshold = module.Ethreshold
        self.Emaximum = module.Emaximum
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


    def CrossSectionFactors_SD66(self, ER, mx, fp, fn, delta):
        #ffelemQ = FFElementQ(self.Z)
        mu_p = ProtonMass * mx / (ProtonMass + mx)
        return self.mass_fraction * 3./(8.*mu_p**6) * self.mT**2 * 1e-12 * ER**2 * \
            (SpeedOfLight/v0bar)**4 * \
            mPhiRef**4 / (4. * self.mT**2 * (ER + self.mPhi**2/(2 * self.mT))**2) * \
            (FF66normlalized(ER, self.A, self.Z, self.mT, 0, 0) + \
            FF66normlalized(ER, self.A, self.Z, self.mT, 0, 1) * 2 * fn/fp + \
            FF66normlalized(ER, self.A, self.Z, self.mT, 1, 1) * (fn/fp)**2) 
        '''
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

    def QuenchingFactor(self,ER):
        return [self.QuenchingFactorList[i](ER) for i in range(self.numT)]        
        
    #TODO extend to other than Gaussian
    def Resolution(self, Eee, qER):
        return map(lambda mu: Gaussian(Eee, mu, self.EnergyResolution(mu)), qER)
        
    def DifferentialResponseSHM(self, Eee, ER, mx, fp, fn, delta): 
        self.count_diffresponse_calls += 1
        q = np.array(self.QuenchingFactor(ER))
        qER = q * ER
        vmin = VMin(ER, self.mT, mx, delta)
        efficiency = np.array(map(self.Efficiency, qER))
        r_list = 1e-6 * kilogram * self.CrossSectionFactors(ER, mx, fp, fn, delta) * \
            efficiency * \
            self.Resolution(Eee, qER) * eta0Maxwellian(vmin, vobs, v0bar, vesc)
        return r_list.sum()
        
    def ResponseSHM_Dirac(self, ER, Eee1, Eee2, mx, fp, fn, delta): 
        q = np.array(self.QuenchingFactor(ER))
        qER = q * ER
        vmin = VMin(ER, self.mT, mx, delta)
        efficiency = map(self.Efficiency,qER)
        integrated_delta = map(lambda e: 1. if Eee1 <= e < Eee2 else 0., qER)
        r_list = 1.e-6 * kilogram * self.CrossSectionFactors(ER, mx, fp, fn, delta) * \
            efficiency * \
            integrated_delta * eta0Maxwellian(vmin, vobs, v0bar, vesc)
        self.count_response_calls += 1
        return r_list.sum()
#       print(ER, " ", r)
#        return r
        
    def ResponseSHM_Other(self, ER, Eee1, Eee2, mx, fp, fn, delta):
        self.count_response_calls += 1
        return integrate.quad(self.DifferentialResponseSHM, Eee1, Eee2, \
            args=(ER, mx, fp, fn, delta), epsrel = 1e-3, epsabs = 0)[0]
#        print(ER, " ", r)
#        return r
  
    def IntegratedResponseSHM_Dirac(self, Eee1, Eee2, mx, fp, fn, delta):
        vmax = vesc + vobs
        muT = self.mT * mx / (self.mT + mx)
        vdelta = SpeedOfLight / 500. * np.sqrt(delta / 2. / muT) if delta > 0 \
            else np.array([0] * self.numT)
        ER_plus_list = map(lambda i, j: ERecoilBranch(vmax, i, mx, delta, 1) \
            if j < vmax else 0., self.mT, vdelta)
        ER_minus_list = map(lambda i, j: ERecoilBranch(vmax, i, mx, delta, -1) \
            if j < vmax else 1.e6, self.mT, vdelta)
        ER_plus = np.max(ER_plus_list)
        ER_minus = np.max(ER_minus_list)
        if ER_minus < ER_plus:
            return integrate.quad(self.ResponseSHM_Dirac, ER_minus, ER_plus, \
                args=(Eee1, Eee2, mx, fp, fn, delta))[0]
        else:
            return 0.
            
    '''
    def IntegratedResponseSHM_Other(self, Eee1, Eee2, mx, fp, fn, delta):
        vmax = vesc + vobs
        muT = self.mT * mx / (self.mT + mx)
        vdelta = SpeedOfLight / 500. * np.sqrt(delta / 2. / muT) if delta > 0 \
            else np.array([0] * self.numT)
        ER_plus_list = map(lambda i, j: ERecoilBranch(vmax, i, mx, delta, 1) \
            if j < vmax else 0., self.mT, vdelta)
        ER_minus_list = map(lambda i, j: ERecoilBranch(vmax, i, mx, delta, -1) \
            if j < vmax else 1.e6, self.mT, vdelta)
        ER_plus = np.max(ER_plus_list)
        ER_minus = np.max(ER_minus_list)
        if ER_minus < ER_plus:
            integr = integrate.dblquad(self.DifferentialResponseSHM, ER_minus, ER_plus, \
                lambda Eee: Eee1, lambda Eee: Eee2, \
                args=(mx, fp, fn, delta), epsrel = 1e-3, epsabs = 0)
            print(integr)
            return integr[0]
        else:
            return 0.
    '''      
    def IntegratedResponseSHM_Other(self, Eee1, Eee2, mx, fp, fn, delta):
        vmax = vesc + vobs
        muT = self.mT * mx / (self.mT + mx)
        vdelta = SpeedOfLight / 500. * np.sqrt(delta / 2. / muT) if delta > 0 \
            else np.array([0] * self.numT)
        ER_plus_list = map(lambda i, j: ERecoilBranch(vmax, i, mx, delta, 1) \
            if j < vmax else 0., self.mT, vdelta)
        ER_minus_list = map(lambda i, j: ERecoilBranch(vmax, i, mx, delta, -1) \
            if j < vmax else 1.e6, self.mT, vdelta)
        ER_plus = np.max(ER_plus_list)
        ER_minus = np.max(ER_minus_list)
        if ER_minus < ER_plus:
            integr = integrate.quad(self.ResponseSHM_Other, ER_minus, ER_plus, \
                args=(Eee1, Eee2, mx, fp, fn, delta), epsrel = 1e-2, epsabs = 0)
            print(integr)
            return integr[0]
        else:
            return 0.
    
            
    def MaximumGapUpperBoundSHM(self, mx, fp, fn, delta, output_file):
        print("mx = ", mx)
        xtable = np.array(map(lambda i, j: self.IntegratedResponseSHM(i, j, mx, fp, fn, delta), \
            self.ElistMaxGap[:-1], self.ElistMaxGap[1:]))
        mu_scaled = xtable.sum()
        x_scaled = np.max(xtable)
        mu_over_x = mu_scaled / x_scaled
        print("mx = ", mx, "   mu_over_x = ", mu_over_x)
        print("xtable = ", xtable)
        y_guess = np.real(-lambertw(-0.1 / mu_over_x, -1))
#        print("y_guess = ", y_guess)
        y = fsolve(lambda x: MaximumGapC0scaled(x, mu_over_x) - 0.9, y_guess)
        result =  y / x_scaled / self.Exposure
        to_print = np.log10(np.array([[mx, result[0]]]))
        with open(output_file,'a') as f_handle:
            np.savetxt(f_handle, to_print)
        return result
        
    def MaximumGapLimit(self, fp, fn, delta, mx_min, mx_max, num_steps, output_file):
        mx_list = np.logspace(np.log10(mx_min), np.log10(mx_max), num_steps)
        upper_limit = np.array(map(lambda mx: self.MaximumGapUpperBoundSHM(mx, fp, fn, delta, output_file), mx_list))
        print(mx_list)
        print(upper_limit)
        return np.log10(np.transpose([mx_list, upper_limit.flatten()]))
        