# -*- coding: utf-8 -*-
"""
Created on Sat Feb 28 21:41:09 2015

@author: andreea
"""

def DM_mass_range(exper_name, delta, mPhi = 1000., quenching = None):
    if exper_name == "superCDMS":
        num_steps = 60
        mx_range_options = {(0, 1000.): (5, 100, num_steps),
#        mx_range_options = {(0, 1000.): (7, 100, num_steps),
                            (0, 0.): (4, 130, num_steps),
                            (-30, 1000.): (2.3, 100, num_steps),
                            (-30, 0.): (2, 100, num_steps),
                            (-50, 1000.): (1.8, 50, num_steps),
                            (50, 1000.): (20, 100, num_steps),
                            (100, 1000.): (30, 100, num_steps),
        }
    elif "LUX" in exper_name:
        num_steps = 30
        mx_range_options = {(0, 1000.): (5.85, 100, num_steps),
#        mx_range_options = {(0, 1000.): (7.6, 100, num_steps),
                            (0, 0.): (5.80, 130, num_steps),
                            (-30, 1000.): (3.95, 100, num_steps),
                            (-30, 0.): (3.95, 100, num_steps),
                            (-50, 1000.): (3.197, 50, num_steps),
                            (50, 1000.): (17.66, 100, num_steps),
                            (100, 1000.): (40, 100, num_steps),
                            (100, 0.): (40, 300, num_steps),
        }
    elif exper_name == "KIMS2012":
        num_steps = 40
        mx_range_options = {(0, 1000.): (10, 100, num_steps),
                            (0, 0.): (10, 130, num_steps),
                            (-30, 1000.): (8, 100, num_steps),
                            (-30, 0.): (10, 100, num_steps),
                            (-50, 1000.): (6, 50, num_steps),
                            (50, 1000.): (17, 100, num_steps),
                            (50, 0.): (17, 100, num_steps),
                            (100, 1000.): (41, 100, num_steps),
                            (100, 0.): (41.08, 300, num_steps),
        }
    elif exper_name == "SIMPLEModeStage2":
        num_steps = 200
        mx_range_options = {(0, 1000.): (4, 100, num_steps),
                            (0, 0.): (4, 130, num_steps),
                            (-30, 1000.): (2, 100, num_steps),
                            (-30, 0.): (2, 100, num_steps),
                            (-50, 1000.): (1.5, 50, num_steps),
                            (50, 1000.): (18, 100, num_steps),
                            (100, 1000.): (30, 100, num_steps),
        }
    elif exper_name == "DAMA2010NaSmRebinned":
        num_steps = 60
        if quenching == 0.4:
            mx_range_options = {(0, 1000.): (5, 15, num_steps),
#            mx_range_options = {(0, 1000.): (8, 15, num_steps),
                                (0, 0.): (6, 20, num_steps),
                                (-30, 1000.): (2, 4, num_steps),
                                (-30, 0.): (2, 3, num_steps),
                                (-50, 1000.): (1.5, 2.5, num_steps),
            }
        else:
            mx_range_options = {(0, 1000.): (6, 20, num_steps),
#            mx_range_options = {(0, 1000.): (10, 20, num_steps),
                                (0, 0.): (7, 30, num_steps),
                                (-30, 1000.): (3, 5, num_steps),
                                (-30, 0.): (2.5, 4, num_steps),
                                (-50, 1000.): (2, 3.5, num_steps),
            }
    elif exper_name == "DAMA2010ISmRebinned":
        num_steps = 60
        if quenching == 0.09:
            mx_range_options = {(0, 1000.): (25, 65, num_steps),
#            mx_range_options = {(0, 1000.): (40, 65, num_steps),
                                 (0, 0.): (35, 90, num_steps),
                                (-30, 1000.): (20, 50, num_steps),
                                (-30, 0.): (30, 50, num_steps),
                                (-50, 1000.): (20, 40, num_steps),
                                (50, 1000.): (30, 80, num_steps),
                                (50, 0.): (40, 100, num_steps),
                                (100, 1000.): (42, 65, num_steps),
                                (100, 0.): (50, 300, num_steps),
#                                (100, 0.): (42, 65, num_steps),
            }
        else:
            mx_range_options = {(0, 1000.): (35, 90, num_steps),
#            mx_range_options = {(0, 1000.): (55, 85, num_steps),
                                (0, 0.): (40, 130, num_steps),
                                (-30, 1000.): (30, 60, num_steps),
                                (-30, 0.): (40, 70, num_steps),
                                (-50, 1000.): (30, 50, num_steps),
                                (50, 1000.): (40, 100, num_steps),
                                (50, 0.): (55, 100, num_steps),
                                (100, 1000.): (50, 90, num_steps),
                                (100, 0.): (100, 300, num_steps),
#                                (100, 0.): (50, 90, num_steps),
            }
    elif exper_name == "DAMA2010NaSmRebinned_TotRateLimit":
        num_steps = 60
        mx_range_options = {(0, 1000.): (3, 20, num_steps),
                            (0, 0.): (3, 20, num_steps),
                            (-30, 1000.): (1, 10, num_steps),
                            (-30, 0.): (1, 10, num_steps),
                            (-50, 1000.): (1, 10, num_steps),
        }
    else:
        num_steps = 30
        mx_range_options = {(0, 1000.): (3, 100, num_steps),
                            (0, 0.): (3, 130, num_steps),
                            (-30, 1000.): (1, 100, num_steps),
                            (-30, 0.): (1, 100, num_steps),
                            (-50, 1000.): (1, 50, num_steps),
                            (50, 1000.): (10, 100, num_steps),
                            (100, 1000.): (30, 100, num_steps),
        }
    return mx_range_options[(delta, mPhi)]

input_list = [(-1/16.4, 0, 1000.), \
        (0, 0, 1000.), (0, -30, 1000.),  (0, -50, 1000.), (0, 50, 1000.), (0, 0, 0.), \
        (0, -30, 0.), (0, 100, 1000.), (0, 50, 0.), (0, 100, 0.)]