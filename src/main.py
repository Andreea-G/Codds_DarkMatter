# -*- coding: utf-8 -*-
"""
Created on Thu Nov 20 22:52:11 2014

@author: Andreea
"""
from __future__ import absolute_import
from __future__ import print_function
from __future__ import division


INPUT_DIR = "Data/"
OUTPUT_DIR = "Output/"

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
import numpy as np
from numpy import pi
from scipy import integrate
from scipy.optimize import fsolve
from scipy.special import lambertw

def CrossSectionFactors_SI(exper, ER, mx, fp, fn, delta):
    #print(exper.name, "SI")
    mu_p = ProtonMass * mx / (ProtonMass + mx)
    return exper.mass_fraction * 1./(2.*mu_p**2) * \
        mPhiRef**4 / (4. * exper.mT**2 * (ER + exper.mPhi**2/(2. * exper.mT))**2) * \
        ((exper.Z + (exper.A - exper.Z) * fn/fp)**2) * exper.FormFactor(ER) 


def CrossSectionFactors_SD66(exper, ER, mx, fp, fn, delta):
    #print(exper.name, "SID66")
    mu_p = ProtonMass * mx / (ProtonMass + mx)
    return exper.mass_fraction * 3./(8.*mu_p**6) * exper.mT**2 * 1e-12 * ER**2 * \
        (SpeedOfLight/v0bar)**4 * \
        mPhiRef**4 / (4. * exper.mT**2 * (ER + exper.mPhi**2/(2 * exper.mT))**2) * \
        (4./3. * (4. * pi)/(2 * exper.J + 1.)) * (exper.SpScaled + exper.SnScaled * fn/fp)**2 * \
        exper.FormFactor(ER) 

#TODO
def CrossSectionFactors_SD44(exper, ER, mx, fp, fn, delta):
    #print(exper.name, "SD44")
    return 0

CrossSectionFactors_options = {'SI' : CrossSectionFactors_SI,
           'SD66' : CrossSectionFactors_SD66,
           'SD44' : CrossSectionFactors_SD44,
}

class Experiment:
    def __init__(self, expername, scattering_type, mPhi = mPhiRef):
        module = import_file(INPUT_DIR + expername + ".py")
        self.name = expername
        self.scattering_type = scattering_type
        self.energy_resolution_type = module.energy_resolution_type       
        self.FF = module.FF[scattering_type]
        
        self.mPhi = mPhi
        self.numT = module.num_target_nuclides
        self.mT = module.target_nuclide_mass_list
        self.A = module.target_nuclide_AZC_list[:,0]
        self.Z = module.target_nuclide_AZC_list[:,1]
        self.mass_fraction = module.target_nuclide_AZC_list[:,2]
        self.J = module.target_nuclide_JSpSn_list[:,0]
        self.SpScaled = module.target_nuclide_JSpSn_list[:,1]
        self.SnScaled = module.target_nuclide_JSpSn_list[:,2]

        self.QuenchingFactorList = module.QuenchingFactorList  
        self.Efficiency = module.Efficiency
        self.Ethreshold = module.Ethreshold
        self.Emaximum = module.Emaximum
        self.ERecoilList = module.ERecoilList
        self.ElistMaxGap = np.append( np.insert( \
            np.array(filter(lambda x: self.Ethreshold < x < self.Emaximum, self.ERecoilList)), \
            0, self.Ethreshold), self.Emaximum)
        self.Exposure = module.Exposure
        
    def FormFactor(self, ER):
        result = FF_options[self.scattering_type][self.FF](ER, self.A, self.mT)
        #print("FF ", result)
        return result

    def CrossSectionFactors(self, ER, mx, fp, fn, delta):
        result = CrossSectionFactors_options[self.scattering_type](self, ER, mx, fp, fn, delta)
        #print("CS ", result)
        return result

    def QuenchingFactor(self,ER):
        return [self.QuenchingFactorList[i](ER) for i in range(self.numT)]        
        
    def DifferentialResponseSHM(self, ER, Eee, mx, fp, fn, delta): 
        q = np.array(self.QuenchingFactor(ER))
        qER = q * ER
        vmin = VMin(ER, self.mT, mx, delta)
        #print("vmin = ",vmin)
        efficiency = map(self.Efficiency,qER)
        #print("eta0 = ", eta0Maxwellian(vmin, vobs, v0bar, vesc))
        r_list = 1e-6 * kilogram * self.CrossSectionFactors(ER, mx, fp, fn, delta) * efficiency * \
            eta0Maxwellian(vmin, vobs, v0bar, vesc)
        return r_list.sum()
        
    def ResponseSHM(self, ER, Eee1, Eee2, mx, fp, fn, delta): 
        q = np.array(self.QuenchingFactor(ER))
        qER = q * ER
        vmin = VMin(ER, self.mT, mx, delta)
        #print("vmin = ",vmin)
        efficiency = map(self.Efficiency,qER)
        #print("eta0 = ", eta0Maxwellian(vmin, vobs, v0bar, vesc))
        integrated_delta = map(lambda e: 1. if Eee1 <= e < Eee2 else 0., qER)
        r_list = 1.e-6 * kilogram * self.CrossSectionFactors(ER, mx, fp, fn, delta) * efficiency * \
            integrated_delta * eta0Maxwellian(vmin, vobs, v0bar, vesc)
        return r_list.sum()
        
    def IntegratedResponseSHM(self, Eee1, Eee2, mx, fp, fn, delta):
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
            return integrate.quad(self.ResponseSHM, ER_minus, ER_plus, args=(Eee1, Eee2, mx, fp, fn, delta))[0]
        else:
            return 0.
            
    def MaximumGapUpperBoundSHM(self, mx, fp, fn, delta):
        xtable = np.array(map(lambda i, j: self.IntegratedResponseSHM(i, j, mx, fp, fn, delta), \
            self.ElistMaxGap[:-1], self.ElistMaxGap[1:]))
        mu_scaled = xtable.sum()
        x_scaled = np.max(xtable)
        mu_over_x = mu_scaled / x_scaled
        print("mx = ", mx, "   mu_over_x = ", mu_over_x)
        print("xtable = ", xtable)
        y_guess = np.real(-lambertw(-0.1 / mu_over_x, -1))
        print("y_guess = ", y_guess)
        y = fsolve(lambda x: MaximumGapC0scaled(x, mu_over_x) - 0.9, y_guess)
        return y / x_scaled / self.Exposure
        
    def MaximumGapLimit(self, fp, fn, delta, mx_min, mx_max, num_steps):
        mx_list = np.logspace(np.log10(mx_min), np.log10(mx_max), num_steps)
        upper_limit = np.array(map(lambda mx: self.MaximumGapUpperBoundSHM(mx, fp, fn, delta), mx_list))
        return np.log10(np.transpose([mx_list, upper_limit.flatten()]))
        
def main():
    exper = "superCDMS"
    scattering_type = 'SD66'
    mPhi = 1000.
    
    cl1 = Experiment(exper, scattering_type, mPhi)
    print('name = ', cl1.name)
    print(cl1.mass_fraction)
    print("A = ", cl1.A)
    print(cl1.J)
    print(cl1.SpScaled)
    print(cl1.SnScaled)
    
    #ER = 6.
    #mx = 10.
    fp = 1.
    fn = 0.
    delta = 0.
    #Eee1 = 2
    #Eee2 = 3
    #print("diff response = ", cl1.DifferentialResponseSHM(ER, Eee1, mx, fp, fn, delta))
    #print("response = ", cl1.ResponseSHM(ER, Eee1, Eee2, mx, fp, fn, delta))
    #print("int response = ", cl1.IntegratedResponseSHM(Eee1, Eee2, mx, fp, fn, delta))
    #print("max gap = ", cl1.MaximumGapUpperBoundSHM(mx, fp, fn, delta))

    mx_min = 6.
    mx_max = 100.
    num_steps = 30
    print("max gap = ", cl1.MaximumGapLimit(fp, fn, delta, mx_min, mx_max, num_steps))
    
if __name__ == '__main__':
    main()

