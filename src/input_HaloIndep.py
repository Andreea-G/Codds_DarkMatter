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
        vmin_range_options = {(9., 0, 1000.): (vmin_step, 1000, vmin_step),
                              (3.5, -50, 1000.): (100, 1000, vmin_step),
        }
    return vmin_range_options[(mx, delta, mPhi)]

input_list = [(9., 1, 0, 1000.), (3.5, 1, -50, 1000.)]
vmin_FoxBand_range = (300, 800, 30)
logeta_FoxBand_percent_range = (0.2, 0.2, 15)
steepness = None