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

name = "PICASSO"
modulated = False

energy_resolution_type = "Dirac" 
# actually Bubble Nucleation, but similar enough to implement like Dirac
EnergyResolution = lambda e: 0.5
FFSD = 'GaussianFFSD'
FFSI = 'HelmFF'
FF = {'SI' : FFSI,
      'SDPS' : FFSD, 
      'SDAV' : FFSD,
}
target_nuclide_AZC_list = np.array([[19, 9, 0.7981563864573104]])
target_nuclide_JSpSn_list = np.array([[1./2, 0.4751 * np.sqrt(3./2 / pi), -0.0087 * np.sqrt(3./2 / pi)]])
target_nuclide_mass_list = np.array([17.6969])
num_target_nuclides = target_nuclide_mass_list.size

QuenchingFactor = lambda e: 1.

Ethreshold = 1.7
Emaximum = 100
ERmaximum = np.inf

Efficiency_ER = lambda er: 1.
alpha = 5.
Efficiency = lambda e, er: 1. - np.exp(alpha * (1. - er/e))

Exposure = 1.  # not needed
ERecoilList = np.array([])

BinSize = 1;
# BinEdges are actually threshold energies
BinEdges_left = np.array([1.723498, 2.900465, 4.098237, 5.813693, 6.896901, 16.334835, 38.841959, 54.882078])
BinEdges_right = 100 * np.ones(BinEdges_left.size)
# BinData are rate values
BinData = np.array([-6.027919, -0.317259, 1.586294, -0.190355, 0., 1.395939, -0.253807, 1.332487])
# BinError are rate errors
BinError = np.array([7.170051, 1.77665, 9.073604, 9.200507, 1.269036, 1.649746, 1.77665, 4.63198])

if ConfidenceLevel == 0.9:
    chiSquared = 13.3616
elif ConfidenceLevel == 0.99:
    chiSquared = 20.0902
else:
    chiSquared = 13.3616
    print("Warning! You asked for a ConfidenceLevel that is not given! Using 0.9 instead.")
