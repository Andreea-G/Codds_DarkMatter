# -*- coding: utf-8 -*-
"""
Created on Wed Nov 19 00:18:55 2014

@author: Andreea
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np
from interp import interp1d
pi = np.pi

name = "KIMS2012"
modulated = False

energy_resolution_type = "Gaussian"

def EnergyResolution(e):
    return 0.582 * np.sqrt(e) + 0.0021 * e

FFSD = 'GaussianFFSD'
FFSI = 'HelmFF'
FF = {'SI': FFSI,
      'SDPS': FFSD,
      'SDAV': FFSD,
      }
target_nuclide_AZC_list = np.array([[127, 53, 2.29522], [133, 55, 2.40375]])
target_nuclide_JSpSn_list = \
    np.array([[5./2, 0.309 * np.sqrt(21./10 / pi), .075 * np.sqrt(21./10 / pi)],
              [7./2, -.370 * np.sqrt(18./7 / pi), .003 * np.sqrt(18./7 / pi)]])
target_nuclide_mass_list = np.array([118.211, 123.801])
num_target_nuclides = target_nuclide_mass_list.size

def QuenchingFactor(e):
    return 0.05 * np.ones_like(e)

Ethreshold = 3
Emaximum = 11
ERmaximum = 300.

def Efficiency_ER(er):
    return np.ones_like(er)

Efficiency_interp = interp1d(np.array([3., 10., 11.]), np.array([0.3, 1.1, 1.1]))
def Efficiency(e):
    return Efficiency_interp(e) if Ethreshold <= e < Emaximum \
            else np.array(0.)

Exposure = 24524.3
ERecoilList = np.array([])

BinSize = 1
BinEdges = np.array(range(Ethreshold, Emaximum + BinSize))
BinEdges_left = BinEdges[:-1]
BinEdges_right = BinEdges[1:]
BinData = np.array([0., 0., 0., 0., 0., 0.006678, 0.007997, 0.005212])
BinError = np.array([0.008209, 0.007622, 0.007036, 0.003078, 0.007622, 0.002345,
                     0.011286, 0.010114])
Expected_limit = np.array([0.016759777, 0.016871508, 0.014636871, 0.008156424,
                           0.015810056, 0.016145252, 0.028268157, 0.026368715])
# from Fig 4 of 1204.2646v2
