# -*- coding: utf-8 -*-
"""
Created on Sun Mar  1 21:28:44 2015

@author: Andreea
"""
from globalfnc import VMin


def Vmin_range(exper_name, mx, delta, mPhi=1000., quenching=None, FOX_METHOD=False):
    vmin_range_options = {}
    if exper_name == "CDMSSi2012" and FOX_METHOD:
        vmin_step = vmin_min = 1
        vmin_max = 1000
    elif "LUX" in exper_name:
        vmin_step = vmin_min = 1
        vmin_max = 1000
        vmin_range_options = {(9., 0, 1000.): (450, vmin_max, vmin_step),
                              (3.5, -50, 1000.): (550, vmin_max, vmin_step),
                              (1.3, -200, 1000.): (750, vmin_max, vmin_step)
                              }
    elif "DAMA" in exper_name:
        vmin_step = vmin_min = 1
        if "Na" in exper_name:
            mT = 21.4148
            qT = 0.3
        else:
            mT = 118.211
            qT = 0.06
        Emax = 20
        vmin_max = round(VMin(Emax/qT, mT, mx, delta)) + 200
        if delta == 0:
            Ethreshold = 2.
            qT = 0.09
            vmin_step = VMin(Ethreshold/qT, mT, mx, delta)
        vmin_range_options = {(30.14, 0, 1000.): (vmin_step, vmin_max, 1000),
                              (47.35, 0, 1000.): (vmin_step, vmin_max, 600)
                              }
    elif "KIMS" in exper_name:
        vmin_min = vmin_step = 1
        vmin_max = 1000
    else:
        vmin_step = vmin_min = 1
        vmin_max = 1000
    default = (vmin_min, vmin_max, vmin_step)
    return vmin_range_options.get((mx, delta, mPhi), default)


def Steepness(exper_name, mx, delta, mPhi=1000.):
    if exper_name != "CDMSSi2012":
        return None
    steepness_options = {(9., 0, 1000.): (1., 2.5, 1),
                         (3.5, -50, 1000.): (0.2, 1, 1),
                         (1.3, -200, 1000.): (0.1, 0.6, 1),
                         }
    default = (1.5, 2.5, 1)
    return steepness_options.get((mx, delta, mPhi), default)


def Logeta_guess(exper_name, mx, delta, mPhi=1000.):
    if exper_name != "CDMSSi2012":
        return None
    logeta_options = {(9., 0, 1000.): -21,
                      (3.5, -50, 1000.): -24,
                      (1.3, -200, 1000.): -22,
                      }
    return logeta_options[(mx, delta, mPhi)]


def Vmin_EHIBand_range(exper_name, mx, delta, mPhi=1000.):
    if exper_name != "CDMSSi2012":
        return None
    options = {(9., 0, 1000.): (0, 1000, 100),
               (3.5, -50, 1000.): (0, 1000, 80),
               (1.3, -200, 1000.): (0, 1000, 80),
               }
    return options[(mx, delta, mPhi)]


# input of the form (mx, fn, delta, mPhi)
input_list = [(9., 1, 0., 1000.), (9., -0.8, 0., 1000.), (9., -0.7, 0., 1000.),  # 0 - 2
              (3.5, 1, -50, 1000.), (3.5, -0.8, -50, 1000.), (3.5, -0.7, -50, 1000.),  # 3 - 5
              (1.3, 1, -200, 1000), (1.3, -0.8, -200, 1000), (1.3, -0.7, -200, 1000),  # 6 - 8
              (9., 0, 0., 1000.), (33., 0, 0, 1000), (47., 0, 0, 1000),  # 9 - 11
              (7, 0, 0, 1000), (30.14, 0, 0, 1000), (47.35, 0, 0, 1000),  # 12 - 14, SDPS
              (38, 0, 50, 1000.), (45, 0, 100, 1000),  # 15 - 16, SDPS
              (40, 0, 50, 1000), (52, 0, 100, 1000), (80, 0, 100, 0)]  # 17 - 19, SDAV
logeta_EHIBand_percent_range = (0.2, 0.2, 50)
