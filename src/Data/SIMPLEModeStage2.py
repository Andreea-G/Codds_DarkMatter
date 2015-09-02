"""
Copyright (c) 2015 Andreea Georgescu

Created on Wed Nov 19 00:18:55 2014

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import numpy as np
from globalfnc import ConfidenceLevel
pi = np.pi

name = "SIMPLEModeStage2"
modulated = False

energy_resolution_type = "Dirac"

def EnergyResolution(e):
    return np.ones_like(e)

FFSD = 'GaussianFFSD'
FFSI = 'HelmFF'
FF = {'SI': FFSI,
      'SDPS': FFSD,
      'SDAV': FFSD,
      }

target_nuclide_AZC_list = \
    np.array([[12, 6, 0.153648], [13, 6, 0.00186884], [35, 17, 0.171531],
              [37, 17, 0.0579854], [19, 9, 0.614966]])
target_nuclide_JSpSn_list = \
    np.array([[0, 0, 0],
              [1./2, -0.026 * np.sqrt(3./2 / pi),  0.115 * np.sqrt(3./2 / pi)],
              [3./2, -0.051 * np.sqrt(5./3 / pi),  0.0088 * np.sqrt(5./3 / pi)],
              [3./2, -0.051 * np.sqrt(5./3 / pi),  0.0088 * np.sqrt(5./3 / pi)],
              [1./2, 0.4751 * np.sqrt(3./2 / pi),  -0.0087 * np.sqrt(3./2 / pi)]])
target_nuclide_mass_list = np.array([11.1779, 12.1125, 32.5733, 34.4335, 17.6969])
num_target_nuclides = target_nuclide_mass_list.size

def QuenchingFactor(e):
    return np.ones_like(e)

Ethreshold = 8
Emaximum = np.inf
ERmaximum = np.inf

def Efficiency(e, er):
    return 1.

def Efficiency_ER(e):
    return np.array(0.97 * (1 - np.exp(-4.2 * (e/8 - 1)))) \
            if Ethreshold <= e else np.array(0.)

Exposure = 6.71
ERecoilList = np.array([])

if ConfidenceLevel == 0.9:
    Expected_limit = 2.39  # upper limit of expected events for 2 observed events at 90 % CL
elif ConfidenceLevel == 0.99:
    Expected_limit = 4.715
else:
    Expected_limit = 2.39  # upper limit of expected events for 2 observed events at 90 % CL
    print("Warning! You asked for a ConfidenceLevel that is not given! Using 0.9 instead.")
