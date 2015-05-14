# -*- coding: utf-8 -*-
"""
Created on Fri May  8 15:04:52 2015

@author: Andreea
"""

from __future__ import absolute_import
from __future__ import division
import numpy as np
#from scipy.interpolate import interp1d
from interp import interp1d
#from scipy.optimize import curve_fit
#import matplotlib.pyplot as plt

pi = np.pi

name = "CDMSSiGe_artif"
modulated = False

energy_resolution_type = "Gaussian"

def EnergyResolution(e):
    return np.sqrt(0.085849 + 0.003136 * e)

FFSD = 'GaussianFFSD'
FFSI = 'HelmFF'
FF = {'SI': FFSI,
      'SDPS': FFSD,
      'SDAV': FFSD,
      }
target_nuclide_AZC_list = np.array([[28, 14, 0.918663943428171], [29, 14, 0.04833558589888038],
                                    [30, 14, 0.03300047067294847],
                                    [70., 32., 0.19608], [72., 32., 0.27040], [73., 32., 0.07790],
                                    [74., 32., 0.37378], [76., 32., 0.08184]])
target_nuclide_JSpSn_list = \
    np.array([[0, 0, 0], [1./2, -0.0019 * np.sqrt(3./(2 * pi)), .1334 * np.sqrt(3./(2 * pi))],
              [0, 0, 0],
              [0., 0., 0.], [0., 0., 0.],
              [9./2, 0.0392517 * np.sqrt(((2*9./2 + 1)*(9./2 + 1))/(4*pi*9./2)),
               .375312 * np.sqrt(((2*9./2 + 1)*(9./2 + 1))/(4*pi*9./2))],
              [0., 0., 0.], [0., 0., 0.]])
target_nuclide_mass_list = np.array([26.0603, 26.9914, 27.9204,
                                     65.134, 66.995, 67.9278, 68.8571, 70.7203])
num_target_nuclides = target_nuclide_mass_list.size

def QuenchingFactor(e):
    return np.ones_like(e)

Ethreshold = 7.
Emaximum = 100.
ERmaximum = 100.

x = np.array([6.88073000e+00, 6.93807000e+00, 7.33945000e+00, 7.91284000e+00, 8.25688000e+00,
              8.48624000e+00, 9.11697000e+00, 9.74771000e+00, 1.05505000e+01, 1.14106000e+01,
              1.23853000e+01, 1.33028000e+01, 1.41628000e+01, 1.47362000e+01, 1.49083000e+01,
              1.49083000e+01, 1.49656000e+01, 1.56537000e+01, 1.63991000e+01, 1.72592000e+01,
              1.82913000e+01, 1.91514000e+01, 1.97821000e+01, 1.98968000e+01, 2.12156000e+01,
              2.26491000e+01, 2.97592000e+01, 3.00459000e+01, 3.19954000e+01, 3.51491000e+01,
              3.83028000e+01, 4.18578000e+01, 4.67890000e+01, 5.36697000e+01, 1.00057000e+02])
y = np.array([4.17071000e-02, 8.43841000e-02, 1.10572000e-01, 1.29971000e-01, 1.41610000e-01,
              1.55189000e-01, 1.70708000e-01, 1.84287000e-01, 1.97866000e-01, 2.07565000e-01,
              2.16295000e-01, 2.22114000e-01, 2.28904000e-01, 2.32784000e-01, 2.38603000e-01,
              2.49273000e-01, 2.54122000e-01, 2.60912000e-01, 2.66731000e-01, 2.73521000e-01,
              2.79340000e-01, 2.84190000e-01, 2.86130000e-01, 3.23957000e-01, 3.32687000e-01,
              3.40446000e-01, 3.80213000e-01, 3.90883000e-01, 4.02522000e-01, 4.15131000e-01,
              4.23860000e-01, 4.31620000e-01, 4.41319000e-01, 4.52958000e-01, 5.19884000e-01])

Efficiency_interp = interp1d(x, y)

def Efficiency(e):
    return Efficiency_interp(e) if e >= Ethreshold else 0.

def Efficiency_ER(er):
    return np.ones_like(er)

Exposure = 140.2
ERecoilList = np.array([8.2, 9.5, 12.3])

mu_BKG_i = np.array([0.0185176, 0.0203303, 0.0217355])
NBKG = 0.7
