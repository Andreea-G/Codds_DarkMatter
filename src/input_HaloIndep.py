# -*- coding: utf-8 -*-
"""
Created on Sun Mar  1 21:28:44 2015

@author: Andreea
"""


def Vmin_range(exper_name, mx, delta, mPhi=1000., quenching=None, FOX_METHOD=False):
    if exper_name == "CDMSSi2012" and FOX_METHOD:
        vmin_step = 0.5
        vmin_range_options = {(9., 0, 1000.): (vmin_step, 2000, vmin_step),
                              (3.5, -50, 1000.): (vmin_step, 2000, vmin_step),
                              (1.3, -200, 1000.): (vmin_step, 2000, vmin_step)
                              }
    elif "LUX" in exper_name:
        vmin_step = 5
        vmin_range_options = {(9., 0, 1000.): (450, 1000, vmin_step),
                              (3.5, -50, 1000.): (550, 1000, vmin_step),
                              (1.3, -200, 1000.): (750, 1000, vmin_step)
                              }
    else:
        vmin_step = 5
        vmin_range_options = {}
    default = (vmin_step, 1000, vmin_step)
    return vmin_range_options.get((mx, delta, mPhi), default)


def Steepness(exper_name, mx, delta, mPhi=1000.):
    # defaults are
    # (steepness_vmin, steepness_vmin_center, steepness_logeta) = (1.5, 2.5, 1)
    if exper_name != "CDMSSi2012":
        return None
    steepness_options = {(9., 0, 1000.): (1., 2.5, 1),
                         (3.5, -50, 1000.): (0.2, 1, 1),
                         (1.3, -200, 1000.): (0.1, 0.6, 1),
                         }
    return steepness_options[(mx, delta, mPhi)]


def Logeta_guess(exper_name, mx, delta, mPhi=1000.):
    # defaults are
    # (steepness_vmin, steepness_vmin_center, steepness_logeta) = (1.5, 2.5, 1)
    if exper_name != "CDMSSi2012":
        return None
    logeta_options = {(9., 0, 1000.): -21,
                      (3.5, -50, 1000.): -24,
                      (1.3, -200, 1000.): -22,
                      }
    return logeta_options[(mx, delta, mPhi)]


def Vmin_FoxBand_range(exper_name, mx, delta, mPhi=1000.):
    # defaults are
    # (steepness_vmin, steepness_vmin_center, steepness_logeta) = (1.5, 2.5, 1)
    if exper_name != "CDMSSi2012":
        return None
    options = {(9., 0, 1000.): (0, 1000, 100),
               (3.5, -50, 1000.): (0, 1000, 80),
               (1.3, -200, 1000.): (0, 1000, 80),
               }
    return options[(mx, delta, mPhi)]


# input of the form (mx, fn, delta, mPhi)
input_list = [(9., -0.8, 0., 1000.), (3.5, -0.8, -50, 1000.), (1.3, -0.8, -200, 1000),  # 0 - 2
              (9., 1, 0, 1000.), (9., -0.7, 0, 1000.), (9, 0, 0, 1000),  # 3 - 5
              (3.5, 1, -50, 1000.), (33., 0, 0, 1000), (47., 0, 0, 1000),  # 6 - 8
              (7, 0, 0, 1000), (34.541, 0, 0, 1000), (57.35, 0, 0, 1000),  # 9 - 11, SDPS
              (38, 0, 50, 1000.), (45, 0, 100, 1000),  # 12 - 13, SDPS
              (58, 0, 50, 1000), (52, 0, 100, 1000), (80, 0, 100, 0)]  # 14 - 16, SDAV
logeta_FoxBand_percent_range = (0.2, 0.2, 50)
