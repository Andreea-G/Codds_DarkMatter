# -*- coding: utf-8 -*-
"""
Created on Wed Nov 19 00:18:55 2014

@author: Andreea
"""
from __future__ import absolute_import
from __future__ import division
import numpy as np
pi = np.pi
#from scipy.interpolate import interp1d
from interp import interp1d

name = "KIMS2012"

energy_resolution_type = "Gaussian"
EnergyResolution = lambda e: 0.582 * np.sqrt(e) + 0.0021 * e
FFSD = 'GaussianFFSD'
FFSI = 'HelmFF'
FF = {'SI' : FFSI,
      'SD66' : FFSD, 
      'SD44' : FFSD,
}
target_nuclide_AZC_list = np.array([[127, 53, 2.29522], [133, 55, 2.40375]])
target_nuclide_JSpSn_list = np.array([[5./2, 0.309 * np.sqrt(21./10 / pi), .075 * np.sqrt(21./10 / pi)], \
    [7./2, -.370 * np.sqrt(18./7 / pi), .003 * np.sqrt(18./7 / pi)]])
target_nuclide_mass_list = np.array([118.211, 123.801])
num_target_nuclides = target_nuclide_mass_list.size

QuenchingFactor = lambda e: 0.05

Ethreshold = 3
Emaximum = 11
ERmaximum = 300.

Efficiency_ER = lambda er: 1.

Efficiency_interp = interp1d(np.array([3., 10., 11.]), np.array([0.3, 1.1, 1.1]))
Efficiency = lambda e: Efficiency_interp(e) if Ethreshold <= e < Emaximum \
    else np.array(0.)

Exposure = 24524.3
ERecoilList = np.array([])

BinSize = 1;
BinEdges = np.array(range(Ethreshold, Emaximum + BinSize))
BinData = np.array([0., 0., 0., 0., 0., 0.006678, 0.007997, 0.005212])
BinError = np.array([0.008209, 0.007622, 0.007036, 0.003078, 0.007622, 0.002345, \
    0.011286, 0.010114])

chiSquared = 13.3616