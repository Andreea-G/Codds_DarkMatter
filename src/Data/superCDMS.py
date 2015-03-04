# -*- coding: utf-8 -*-
"""
Created on Wed Nov 19 00:18:55 2014

@author: Andreea
"""
from __future__ import absolute_import
from __future__ import division
import numpy as np
pi = np.pi
from scipy.interpolate import interp1d

name = "superCDMS"
modulated = False

energy_resolution_type = "Dirac"
EnergyResolution = lambda e: np.ones(e.size)
#energy_resolution_type = "Gaussian"
#EnergyResolution = lambda e: np.sqrt(0.293**2 + 0.056**2 * e)
#EnergyResolution = 0.001
FFSD = 'GaussianFFSD'
FFSI = 'HelmFF'
FF = {'SI' : FFSI,
      'SDPS' : FFSD, 
      'SDAV' : FFSD,
}
target_nuclide_AZC_list = np.array([[70., 32., 0.19608], [72., 32., 0.27040], [73., 32., 0.07790], \
    [74., 32., 0.37378], [76., 32., 0.08184]])
target_nuclide_JSpSn_list = np.array([[0., 0., 0.], [0., 0., 0.], \
    [9./2, 0.0392517 * np.sqrt(((2*9./2 + 1)*(9./2 + 1))/(4*pi*9./2)), .375312 * np.sqrt(((2*9./2 + 1)*(9./2 + 1))/(4*pi*9./2))], \
    [0., 0., 0.], [0., 0., 0.]])
target_nuclide_mass_list = np.array([65.134, 66.995, 67.9278, 68.8571, 70.7203])

num_target_nuclides = target_nuclide_mass_list.size

QuenchingFactor = lambda e: np.ones(e.size)
Ethreshold = 1.63799
Emaximum = 10.0011
ERmaximum = 10.0011

Efficiency_interp = interp1d(np.array([1.63799, 1.93525, 2.35928, 2.37871, 3.12938, 3.15831, 3.8895, 3.90877, \
    4.2841, 4.30358, 4.63016, 4.64942, 5.38539, 5.4095, 5.78968, 6.15036, 6.16481, 6.8911, 6.92511, \
    9.16257, 9.18213, 10.0011]), \
    np.array([0.044225, 0.071339, 0.086737, 0.105692, 0.112107, 0.196045, 0.19975, 0.260222, 0.26388, \
    0.268395, 0.275658, 0.339739, 0.366008, 0.43731, 0.44819, 0.459066, 0.506, 0.514216, 0.543101, \
    0.544292, 0.529854, 0.532668]))
    
Efficiency = lambda e, er: np.ones(er.size)
Efficiency_ER = lambda e: Efficiency_interp(e) if Ethreshold <= e < Emaximum else np.array(0.)
#Efficiency_ER = lambda ER: np.array([Efficiency_interp(e) if Ethreshold <= e < Emaximum else 0. for e in np.array([1]) * ER])

Exposure = 577.
ERecoilList = np.array([1.7, 1.8, 1.9, 1.9, 2.3, 2.7, 3.0, 5.8, 7.0, 7.8, 9.4])

            