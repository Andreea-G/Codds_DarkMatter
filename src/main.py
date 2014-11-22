# -*- coding: utf-8 -*-
"""
Created on Thu Nov 20 22:52:11 2014

@author: Andreea
"""
from __future__ import absolute_import

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


def CrossSectionFactors_SI(exper, ER, mx, fp, fn, delta):
    print exper.name, "SI"
    mu_p = ProtonMass*mx/(ProtonMass + mx)
    return exper.mass_fraction * 1./(2.*mu_p**2) * \
        mPhiRef**4 / (4 * exper.mT**2 * (ER + exper.mPhi**2/(2 * exper.mT))**2) * \
        ((exper.Z + (exper.A - exper.Z) * fn/fp)**2) * exper.FormFactor(ER) 


def CrossSectionFactors_SD66(exper, ER, mx, fp, fn, delta):
    print exper.name, "SID66"
    mu_p = ProtonMass*mx/(ProtonMass + mx)
    return exper.mass_fraction * 3./(8.*mu_p**6) * exper.mT**2 * 1e-12 * ER**2 * \
        (SpeedOfLight/v0bar)**4 * \
        mPhiRef**4 / (4 * exper.mT**2 * (ER + exper.mPhi**2/(2 * exper.mT))**2) * \
        (4./3. * (4. * pi)/(2 * exper.J + 1.)) * (exper.SpScaled + exper.SnScaled * fn/fp)**2 * \
        exper.FormFactor(ER) 

#TODO
def CrossSectionFactors_SD44(exper, ER, mx, fp, fn, delta):
    print exper.name, "SD44"
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
        self.mT = module.target_nuclide_mass_list
        self.A = module.target_nuclide_AZC_list[:,0]
        self.Z = module.target_nuclide_AZC_list[:,1]
        self.mass_fraction = module.target_nuclide_AZC_list[:,2]
        self.J = module.target_nuclide_JSpSn_list[:,0]
        self.SpScaled = module.target_nuclide_JSpSn_list[:,1];
        self.SnScaled = module.target_nuclide_JSpSn_list[:,2];
        
        
    def FormFactor(self, ER):
        result = FF_options[self.scattering_type][self.FF](ER, self.A, self.mT)
        print "FF ", result
        return result

    def CrossSectionFactors(self, ER, mx, fp, fn, delta):
        result = CrossSectionFactors_options[self.scattering_type](self, ER, mx, fp, fn, delta)
        print "CS ", result
        return result
                  
        
def main():
    exper = "superCDMS"
    scattering_type = 'SD66'

    mPhi = 0.01
    cl1 = Experiment(exper, scattering_type, mPhi)
    print 'name = ', cl1.name
    print cl1.mass_fraction
    print "A = ", cl1.A
    print cl1.J
    print cl1.SpScaled 
    print cl1.SnScaled 
    
    ER = 10
    mx = 1
    fp = 1
    fn = 0
    delta = 0
    print cl1.CrossSectionFactors(ER, mx, fp, fn, delta)
    
if __name__ == '__main__':
    main()

