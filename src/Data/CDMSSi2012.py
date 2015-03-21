# -*- coding: utf-8 -*-
"""
Created on Wed Nov 19 00:18:55 2014

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

name = "CDMSSi2012"
modulated = False

energy_resolution_type = "Gaussian"
EnergyResolution = lambda e: np.sqrt(0.085849 + 0.003136 * e)
FFSD = 'GaussianFFSD'
FFSI = 'HelmFF'
FF = {'SI': FFSI,
      'SDPS': FFSD,
      'SDAV': FFSD,
      }
target_nuclide_AZC_list = np.array([[28, 14, 0.918663943428171], [29, 14, 0.04833558589888038],
                                    [30, 14, 0.03300047067294847]])
target_nuclide_JSpSn_list = np.array([[0, 0, 0], [1./2, -0.0019 * np.sqrt(3./(2 * pi)),
                                      .1334 * np.sqrt(3./(2 * pi))], [0, 0, 0]])
target_nuclide_mass_list = np.array([26.0603, 26.9914, 27.9204])
num_target_nuclides = target_nuclide_mass_list.size

QuenchingFactor = lambda e: np.ones(e.size)

Ethreshold = 7.
Emaximum = 100.
ERmaximum = 100.

#x0 = np.array([6.84211, 7.89474, 8.94737, 10., 11.0526, 12.6316, 14.0351, 14.9123,
#    15.7895, 17.193, 19.1228, 20.1754, 22.1053, 23.3333, 25.614, 27.3684, 29.1228, 30., 31.7544, 35.2632,
#    37.8947, 40.7018, 43.5088, 47.0175, 50., 52.9825, 57.0175, 60.3509, 63.3333, 66.1404, 70., 74.2105,
#    77.5439, 81.9298, 85.7895, 89.6491, 93.3333, 98.0702, 100.351])
#
#y0 = np.array([0.084286, 0.13, 0.167143, 0.19, 0.204286, 0.215714, 0.221429, 0.247143, 0.258571, 0.272857,
#    0.275714, 0.324286, 0.327143, 0.344286, 0.355714, 0.361429, 0.372857, 0.392857, 0.398571, 0.41, 0.424286,
#    0.432857, 0.435714, 0.444286, 0.452857, 0.455714, 0.461429, 0.47, 0.472857, 0.478571, 0.478571, 0.487143,
#    0.487143, 0.492857, 0.498571, 0.501429, 0.512857, 0.512857, 0.518571])

x = np.array([6.88073000e+00, 6.93807000e+00, 7.33945000e+00, 7.91284000e+00, 8.25688000e+00,
              8.48624000e+00, 9.11697000e+00, 9.74771000e+00, 1.05505000e+01, 1.14106000e+01, 1.23853000e+01, 1.33028000e+01,
              1.41628000e+01, 1.47362000e+01, 1.49083000e+01, 1.49083000e+01, 1.49656000e+01, 1.56537000e+01, 1.63991000e+01,
              1.72592000e+01, 1.82913000e+01, 1.91514000e+01, 1.97821000e+01, 1.98968000e+01, 2.12156000e+01, 2.26491000e+01,
              2.97592000e+01, 3.00459000e+01, 3.19954000e+01, 3.51491000e+01, 3.83028000e+01, 4.18578000e+01, 4.67890000e+01,
              5.36697000e+01, 1.00057000e+02])
y = np.array([4.17071000e-02, 8.43841000e-02, 1.10572000e-01, 1.29971000e-01, 1.41610000e-01, 1.55189000e-01,
              1.70708000e-01, 1.84287000e-01, 1.97866000e-01, 2.07565000e-01, 2.16295000e-01, 2.22114000e-01, 2.28904000e-01,
              2.32784000e-01, 2.38603000e-01, 2.49273000e-01, 2.54122000e-01, 2.60912000e-01, 2.66731000e-01, 2.73521000e-01,
              2.79340000e-01, 2.84190000e-01, 2.86130000e-01, 3.23957000e-01, 3.32687000e-01, 3.40446000e-01, 3.80213000e-01,
              3.90883000e-01, 4.02522000e-01, 4.15131000e-01, 4.23860000e-01, 4.31620000e-01, 4.41319000e-01, 4.52958000e-01, 5.19884000e-01])

Efficiency_interp = interp1d(x, y)
Efficiency = lambda e: Efficiency_interp(e) if e >= Ethreshold else 0.

Efficiency_ER = lambda er: np.ones(er.size)

Exposure = 140.2
ERecoilList = np.array([8.2, 9.5, 12.3])

mu_BKG_i = np.array([0.0185176, 0.0203303, 0.0217355])
NBKG = 0.7
