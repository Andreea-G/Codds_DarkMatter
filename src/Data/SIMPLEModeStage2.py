# -*- coding: utf-8 -*-
"""
Created on Wed Nov 19 00:18:55 2014

@author: Andreea
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import numpy as np
pi = np.pi
#from scipy.interpolate import interp1d
from interp import interp1d
from globalfnc import ConfidenceLevel

name = "SIMPLEModeStage2"
modulated = False

energy_resolution_type = "Dirac"
EnergyResolution = lambda e: 1.
FFSD = 'GaussianFFSD'
FFSI = 'HelmFF'
FF = {'SI' : FFSI,
      'SD66' : FFSD, 
      'SD44' : FFSD,
}
target_nuclide_AZC_list = np.array([[12, 6, 0.153648], [13, 6, 0.00186884], [35, 17, 0.171531], \
    [37, 17, 0.0579854], [19, 9, 0.614966]])
target_nuclide_JSpSn_list = np.array([[0, 0, 0], [1./2, -0.026 * np.sqrt(3./2 / pi),  0.115 * np.sqrt(3./2 / pi)], \
    [3./2, -0.051 * np.sqrt(5./3 / pi),  0.0088 * np.sqrt(5./3 / pi)], \
    [3./2, -0.051 * np.sqrt(5./3 / pi),  0.0088 * np.sqrt(5./3 / pi)], \
    [1./2, 0.4751 * np.sqrt(3./2 / pi),  -0.0087 * np.sqrt(3./2 / pi)]])
target_nuclide_mass_list = np.array([11.1779, 12.1125, 32.5733, 34.4335, 17.6969])
num_target_nuclides = target_nuclide_mass_list.size

QuenchingFactor = lambda e: 1.

Ethreshold = 8
Emaximum = np.inf
ERmaximum = np.inf

Efficiency = lambda e, er: 1.

Efficiency_ER = lambda e: np.array(0.97 * (1 - np.exp(-4.2 * (e/8 - 1)))) if Ethreshold <= e else np.array(0.)

Exposure = 6.71
ERecoilList = np.array([])

if ConfidenceLevel == 0.9:
    Expected_limit = 2.39 # upper limit of expected events for 2 observed events at 90 % CL
elif ConfidenceLevel == 0.99:
    Expected_limit = 4.715
else:
    Expected_limit = 2.39 # upper limit of expected events for 2 observed events at 90 % CL
    print("Warning! You asked for a ConfidenceLevel that is not given! Using 0.9 instead.")
