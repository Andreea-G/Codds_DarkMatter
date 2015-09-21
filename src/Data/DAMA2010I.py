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
import numpy as np
pi = np.pi

name = "DAMA2010I"
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

target_nuclide_AZC_list = np.array([[127, 53, 0.848684]])
target_nuclide_JSpSn_list = \
    np.array([[5./2, 0.309 * np.sqrt(21./10 / pi), .075 * np.sqrt(21./10 / pi)]])
target_nuclide_mass_list = np.array([118.211])
num_target_nuclides = target_nuclide_mass_list.size

def QuenchingFactor(e):
    return 0.09 * np.ones_like(e)

def QuenchingFactorOfEee(e):
    return QuenchingFactor(e)  # since it's a constant function

Ethreshold = 2.
Emaximum = 1000.
ERmaximum = 10000.

def Efficiency_ER(er):
    return np.ones_like(er)

def Efficiency(e):
    return np.array(1.) if Ethreshold <= e < Emaximum else np.array(0.)

Exposure = 1.33 * 1000 * 365.25
ERecoilList = np.array([])

BinEdges = np.array([2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6.,
                     6.5, 7., 7.5, 8., 8.5, 9., 9.5, 10., 10.5,
                     11., 11.5, 12., 12.5, 13., 13.5, 14., 14.5, 15.,
                     15.5, 16., 16.5, 17., 17.5, 18., 18.5, 19., 19.5, 20.])
BinData = np.array([0.00870655, 0.0126448, 0.0107336, 0.00413128, 0.00432433,
                    0.00281854, 0.00405405, 0.00235521, 0.00021236, 0.00011583,
                    0.00075289, -0.00106178, 0.00052123, 0.00247104, -0.00102316,
                    -0.0026641, 0.00376448, -0.0011776, -0.00142857, -0.00121622,
                    0.00146718, -0.00175676, 0.00063707, 0.00135135, -0.00050193,
                    -0.00055985, 0.00127413, 0.00081081, 0.00061776, 0.00204633,
                    0.00098456, 0.00030888, 0.00279923, 0.00021236, 0.00166023,
                    0.00185328])
BinError = np.array([0.00182433, 0.00201738, 0.00203667, 0.0018726, 0.00163127,
                     0.00143823, 0.00144788, 0.00137066, 0.00132239, 0.00128378,
                     0.00135135, 0.00132239, 0.001361, 0.00133205, 0.00132239,
                     0.00137066, 0.00140927, 0.00139962, 0.00140927, 0.00146718,
                     0.00145753, 0.00142857, 0.00147683, 0.00143822, 0.00145753,
                     0.00142857, 0.00139961, 0.00141892, 0.00138996, 0.0013417,
                     0.00132239, 0.00129343, 0.0013224, 0.00123552, 0.00126448,
                     0.00123552])

