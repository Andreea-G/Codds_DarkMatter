# -*- coding: utf-8 -*-
"""
Created on Wed Nov 19 00:18:55 2014

@author: Andreea
"""
from __future__ import absolute_import
import numpy as np
pi = np.pi

myvar = 'me'

FFSD = 'GaussianFFSD'
FFSI = 'HelmFF'
FF = {'SI' : FFSI,
      'SD66' : FFSD, 
}

energy_resolution_type = "Dirac"
target_nuclide_AZC_list = np.array([[70., 32., 0.19608], [72., 32., 0.27040], [73., 32., 0.07790], \
    [74., 32., 0.37378], [76., 32., 0.08184]])
target_nuclide_JSpSn_list = np.array([[0., 0., 0.], [0., 0., 0.], \
    [9./2, 0.0392517 * np.sqrt(((2*9./2 + 1)*(9./2 + 1))/(4*pi*9./2)), .375312 * np.sqrt(((2*9./2 + 1)*(9./2 + 1))/(4*pi*9./2))], \
    [0., 0., 0.], [0., 0., 0.]])
target_nuclide_mass_list = np.array([65.134, 66.995, 67.9278, 68.8571, 70.7203])
