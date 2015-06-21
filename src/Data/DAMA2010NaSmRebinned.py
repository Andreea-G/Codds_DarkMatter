# -*- coding: utf-8 -*-
"""
Created on Wed Nov 19 00:18:55 2014

@author: Andreea
"""
from __future__ import absolute_import
from __future__ import division
import numpy as np
pi = np.pi

name = "DAMA2010NaSmRebinned"
modulated = True

energy_resolution_type = "Gaussian"

def EnergyResolution(e):
    return 0.448 * np.sqrt(e) + 0.0091 * e

FFSD = 'GaussianFFSD'
FFSI = 'HelmFF'
FF = {'SI': FFSI,
      'SDPS': FFSD,
      'SDAV': FFSD,
      }

target_nuclide_AZC_list = np.array([[23, 11, 0.153373]])
target_nuclide_JSpSn_list = \
    np.array([[3./2, 0.2477 * np.sqrt(5./3 / pi), .0198 * np.sqrt(5./3 / pi)]])
target_nuclide_mass_list = np.array([21.4148])
num_target_nuclides = target_nuclide_mass_list.size

def QuenchingFactor(e):
    return 0.4 * np.ones_like(e)

def QuenchingFactorOfEee(e):
    return QuenchingFactor(e)  # since it's a constant function

Ethreshold = 2.
Emaximum = 1000.
ERmaximum = 2500.

def Efficiency_ER(er):
    return np.ones_like(er)

def Efficiency(e):
    return np.array(1.) if Ethreshold <= e < Emaximum else np.array(0.)

Exposure = 1.33 * 1000 * 365.25
ERecoilList = np.array([])

#BinEdges = np.array([2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6., 6.5, 20.])
#BinData = np.array([0.00883333, 0.0126667, 0.0108333, 0.00422222, 0.00438889, 0.00288889,
#                    0.00411111, 0.00238889, 0.000222222, 0.000518722])
#BinError = np.array([0.00183333, 0.002, 0.00194444, 0.00177778, 0.00155556, 0.00144444,
#                     0.00138889, 0.00122222, 0.00138889, 0.00707238])

BinEdges = np.array([2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6., 6.5, 7., 7.5, 8., 8.5, 9., \
                     9.5, 10., 10.5, 11., 11.5, 12., 12.5, 13., 13.5, 14., 14.5, 15., \
                     15.5, 16., 16.5, 17., 17.5, 18., 18.5, 19., 19.5, 20.])
BinData = np.array([0.00883333, 0.0126667, 0.0108333, 0.00422222, 0.00438889, 0.00288889, \
                    0.00411111, 0.00238889, 0.000222222, 0.0000555556, 0.000722222, \
                    -0.001, 0.0005, 0.00244444, -0.001, -0.00266667, 0.00372222, \
                    -0.00111111, -0.00144444, -0.00127778, 0.00144444, -0.00166667, \
                    0.000555556, 0.00138889, -0.000388889, -0.000555556, 0.00133333, \
                    0.000777778, 0.000555556, 0.00205556, 0.00105556, 0.000277778, \
                    0.00288889, 0.000222222, 0.00166667, 0.00188889])
BinError = np.array([0.00183333, 0.002, 0.00194444, 0.00177778, 0.00155556, 0.00144444, \
                     0.00138889, 0.00122222, 0.00138889, 0.00138889, 0.00133333, \
                     0.00127778, 0.00127778, 0.00127778, 0.00144444, 0.00138889, \
                     0.00144444, 0.00144444, 0.00138889, 0.00155556, 0.0015, 0.00138889, \
                     0.0015, 0.00144444, 0.00133333, 0.00138889, 0.00133333, 0.00133333, \
                     0.00133333, 0.00138889, 0.00127778, 0.00133333, 0.00127778, \
                     0.00122222, 0.00116667, 0.00122222])
