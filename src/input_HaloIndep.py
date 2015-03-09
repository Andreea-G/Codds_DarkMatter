# -*- coding: utf-8 -*-
"""
Created on Sun Mar  1 21:28:44 2015

@author: Andreea
"""

def Vmin_range(exper_name, mx, delta, mPhi = 1000., quenching = None):
    if exper_name == "superCDMS":
        vmin_step = 10
        vmin_range_options = {(9., 0, 1000.): (200, 1000, vmin_step),
                              (3.5, -50, 1000.): (vmin_step, 1000, vmin_step),
        }
    if "LUX" in exper_name:
        vmin_step = 50
        vmin_range_options = {(9., 0, 1000.): (200, 1000, vmin_step),
                              (3.5, -50, 1000.): (600, 1000, vmin_step),
        }
    else:
        vmin_step = 0.5
        vmin_range_options = {(9., 0, 1000.): (vmin_step, 2000, vmin_step),
                              (3.5, -50, 1000.): (vmin_step, 2000, vmin_step),
                              (1.3, -200, 1000.): (vmin_step, 2000, vmin_step),
        }
    return vmin_range_options[(mx, delta, mPhi)]

def Steepness(exper_name, mx, delta, mPhi = 1000.):
    # defaults are (steepness_vmin, steepness_vmin_center, steepness_logeta) = (1.5, 2.5, 1)
    if exper_name != "CDMSSi2012":
        return None
    steepness_options = { (9., 0, 1000.): (1., 2.5, 1.5),
                          (3.5, -50, 1000.): (0.2, 1, 1.5),
                          (1.3, -200, 1000.): (0.1, 0.6, 1.5),
    }
    return steepness_options[(mx, delta, mPhi)]

def Logeta_guess(exper_name, mx, delta, mPhi = 1000.):
    # defaults are (steepness_vmin, steepness_vmin_center, steepness_logeta) = (1.5, 2.5, 1)
    if exper_name != "CDMSSi2012":
        return None
    logeta_options = { (9., 0, 1000.): -21,
                       (3.5, -50, 1000.): -24,
                       (1.3, -200, 1000.): -22,
    }
    return logeta_options[(mx, delta, mPhi)]


# input of the form (mx, fn, delta, mPhi)
input_list = [(9., -0.8, 0., 1000.), (3.5, -0.8, -50, 1000.), (1.3, -0.8, -200, 1000), \
    (9., 1, 0, 1000.), (3.5, 1, -50, 1000.)]
vmin_FoxBand_range = (0, 1000, 80)
logeta_FoxBand_percent_range = (0.2, 0.2, 30)
#vmin_FoxBand_range = (300, 800, 20)
#logeta_FoxBand_percent_range = (0.2, 0.2, 10)

