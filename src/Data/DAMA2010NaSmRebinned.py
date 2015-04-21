# -*- coding: utf-8 -*-
"""
Created on Wed Nov 19 00:18:55 2014

@author: Andreea
"""
from __future__ import absolute_import
from __future__ import division
import numpy as np
from interp import interp1d
pi = np.pi

name = "DAMA2010NaSmRebinned"
modulated = True

energy_resolution_type = "Gaussian"

def EnergyResolution(e): return 0.448 * np.sqrt(e) + 0.0091 * e

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
    try:
        return 0.4 * np.ones(len(e))
    except TypeError:
        return np.array(0.4)

def QuenchingFactorOfEee(e):
    return QuenchingFactor(e)  # since it's a constant function

Ethreshold = 2.
Emaximum = 1000.
ERmaximum = 2500.

def Efficiency_ER(er): return 1.

def Efficiency(e): return np.array(1.) if Ethreshold <= e < Emaximum else np.array(0.)

Exposure = 1.33 * 1000 * 365.25
ERecoilList = np.array([])

BinEdges = np.array([2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6., 6.5, 20.])
BinData = np.array([0.00883333, 0.0126667, 0.0108333, 0.00422222, 0.00438889, 0.00288889,
                    0.00411111, 0.00238889, 0.000222222, 0.000518722])
BinError = np.array([0.00183333, 0.002, 0.00194444, 0.00177778, 0.00155556, 0.00144444,
                     0.00138889, 0.00122222, 0.00138889, 0.00707238])
