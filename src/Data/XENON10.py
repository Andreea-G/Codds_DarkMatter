# -*- coding: utf-8 -*-
"""
Created on Wed Dec  3 21:06:05 2014

@author: Andreea
"""

from __future__ import absolute_import
from __future__ import division
import numpy as np
pi = np.pi
#from scipy.interpolate import interp1d
#from interp import interp1d

name = "XENON10"

energy_resolution_type = "Poisson"
EnergyResolution = lambda e: 0.5
FFSD = 'GaussianFFSD'
FFSI = 'HelmFF'
FF = {'SI' : FFSI,
      'SD66' : FFSD, 
}
target_nuclide_AZC_list = np.array([[124, 54, 0.0008966], [126, 54, 0.0008535], [128, 54, 0.018607], \
    [129, 54, 0.25920], [130, 54, 0.040280], [131, 54, 0.21170], [132, 54, 0.27035], [134, 54, 0.10644], [136, 54, 0.09168]])
target_nuclide_JSpSn_list = np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0], \
    [1./2, 0.010 * np.sqrt(3./2 / pi), .329 * np.sqrt(3./2 / pi)], [0, 0, 0], \
    [3./2, -0.009 * np.sqrt(5./2 / pi), -.272 * np.sqrt(5./2 / pi)], [0, 0, 0], [0, 0, 0], [0, 0, 0]])
target_nuclide_mass_list = np.array([115.418, 117.279, 119.141, 120.074, 121.004, 121.937, 122.868, \
    124.732, 126.597])
num_target_nuclides = target_nuclide_mass_list.size

def QuenchingFactor(e):
    #Absolutely no clue what this is...
    k = 0.110
    NexNi = 1.09
    Xi4Ni = 0.032
    eps = 13.8e-3
    Z = 54.
    eps_e = 11.5 * Z**(-7./3) * e
    g = 3 * eps_e**0.15 + 0.7 * eps_e**0.6 + eps_e
    fn = k * g / (1 + k * g)
    Xi = Xi4Ni / (1 + NexNi) * fn * e / 4. / eps
    nenen_gamma = np.log(1 + Xi) / Xi / (1 + NexNi)
    return nenen_gamma * fn / eps
    
Ethreshold = 1.4 * QuenchingFactor(1.4)
Emaximum = 10. * QuenchingFactor(10.)
ERmaximum = 10.

Efficiency_ER = lambda er: np.array(1.) if er >= 1.4 else np.array(0.)
Efficiency = lambda e: 0.94 if e >= Ethreshold else 0.

Exposure = 15.
ERecoilList = np.array(map(lambda e: e * QuenchingFactor(e), [1.86, 1.90, 2.29, 2.38, 2.44, 3.20, 3.77, 4.15, 5.31, 5.32, 6.27, \
6.57, 6.58, 7.08, 7.28, 7.46, 7.75, 7.84, 7.94, 8.37, 9.19, 9.46, 9.98]))

            