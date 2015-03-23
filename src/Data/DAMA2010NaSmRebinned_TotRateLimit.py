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

name = "DAMA2010NaSmRebinned_TotRateLimit"
modulated = False

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

def QuenchingFactor(e): return 0.4

Ethreshold = 2.
Emaximum = 1000.
ERmaximum = 2500.

def Efficiency_ER(er): return 1.

def Efficiency(e): return np.array(1.) if Ethreshold <= e < Emaximum else np.array(0.)

Exposure = 1.33 * 1000 * 365.25
ERecoilList = np.array([])

BinEdges = np.array([2.00139, 2.2527, 2.50331, 2.74973, 3.00174, 3.25235, 3.50366,
                     3.75008, 4.00069, 4.252, 4.50401, 4.74973, 5.00104, 5.25235, 5.50366,
                     5.75008, 6.])
BinData = np.array([0.248111, 0.268326, 0.29697, 0.31601, 0.327703, 0.315959, 0.296845,
                    0.271759, 0.249216, 0.227695, 0.214179, 0.213116, 0.213536, 0.218709,
                    0.213038, 0.219467])
