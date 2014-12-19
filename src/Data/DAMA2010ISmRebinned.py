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

name = "DAMA2010ISmRebinned"
modulated = True

energy_resolution_type = "Gaussian"
EnergyResolution = lambda e: 0.448 * np.sqrt(e) + 0.0091 * e
FFSD = 'GaussianFFSD'
FFSI = 'HelmFF'
FF = {'SI' : FFSI,
      'SD66' : FFSD, 
      'SD44' : FFSD,
}
target_nuclide_AZC_list = np.array([[127, 53, 0.848684]])
target_nuclide_JSpSn_list = np.array([[5./2, 0.309  * np.sqrt(21./10 / pi), .075 * np.sqrt(21./10 / pi)]])
target_nuclide_mass_list = np.array([118.211])
num_target_nuclides = target_nuclide_mass_list.size

QuenchingFactor = lambda e: 0.09

Ethreshold = 2.
Emaximum = 1000.
ERmaximum = 10000.

Efficiency_ER = lambda er: 1.

Efficiency = lambda e: np.array(1.) if Ethreshold <= e < Emaximum else np.array(0.)

Exposure = 1.33 * 1000 * 365.25
ERecoilList = np.array([])

BinEdges = np.array([2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6., 6.5, 20.])
BinData = np.array([0.00883333, 0.0126667, 0.0108333, 0.00422222, 0.00438889, 0.00288889, \
    0.00411111, 0.00238889, 0.000222222, 0.000518722])
BinError = np.array([0.00183333, 0.002, 0.00194444, 0.00177778, 0.00155556, 0.00144444, \
    0.00138889, 0.00122222, 0.00138889, 0.00707238])

