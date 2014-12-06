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

name = "LUX2013zero"

energy_resolution_type = "Gaussian"
EnergyResolution = lambda e: np.sqrt(0.085849 + 0.003136 * e)
FFSD = 'GaussianFFSD'
FFSI = 'HelmFF'
FF = {'SI' : FFSI,
      'SD66' : FFSD, 
      'SD44' : FFSD,       
}
target_nuclide_AZC_list = np.array([[28, 14, 0.918663943428171], [29, 14, 0.04833558589888038], \
    [30, 14, 0.03300047067294847]])
target_nuclide_JSpSn_list = np.array([[0, 0, 0], [1./2, -0.0019 * np.sqrt(3./(2 * pi)), \
    .1334 * np.sqrt(3./(2 * pi))], [0, 0, 0]])
target_nuclide_mass_list = np.array([26.0603, 26.9914, 27.9204])
num_target_nuclides = target_nuclide_mass_list.size

QuenchingFactor = lambda e: 1

Ethreshold = 7.
Emaximum = 12.3
ERmaximum = 12.3

Efficiency_ER = lambda er: 1.

Efficiency_interp = interp1d(np.array([6.84211, 7.89474, 8.94737, 10., 11.0526, 12.6316, 14.0351, 14.9123, \
    15.7895, 17.193, 19.1228, 20.1754, 22.1053, 23.3333, 25.614, 27.3684, 29.1228, 30., 31.7544, 35.2632, \
    37.8947, 40.7018, 43.5088, 47.0175, 50., 52.9825, 57.0175, 60.3509, 63.3333, 66.1404, 70., 74.2105, \
    77.5439, 81.9298, 85.7895, 89.6491, 93.3333, 98.0702, 100.351]), \
    np.array([0.084286, 0.13, 0.167143, 0.19, 0.204286, 0.215714, 0.221429, 0.247143, 0.258571, 0.272857, \
    0.275714, 0.324286, 0.327143, 0.344286, 0.355714, 0.361429, 0.372857, 0.392857, 0.398571, 0.41, 0.424286, \
    0.432857, 0.435714, 0.444286, 0.452857, 0.455714, 0.461429, 0.47, 0.472857, 0.478571, 0.478571, 0.487143, \
    0.487143, 0.492857, 0.498571, 0.501429, 0.512857, 0.512857, 0.518571]))
Efficiency = lambda e: Efficiency_interp(e) if e >= Ethreshold else np.array(0.)

Exposure = 140.2
ERecoilList = np.array([8.2, 9.5, 12.3])

            