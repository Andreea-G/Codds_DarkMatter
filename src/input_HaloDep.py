# -*- coding: utf-8 -*-
"""
Created on Sat Feb 28 21:41:09 2015

@author: Andreea
"""


def DM_mass_range(exper_name, delta, mPhi=1000., quenching=None):
    """ Range and number of steps for the DM mass.
    Input:
        exper_name: string
            Name of experiment.
        delta: float
            DM mass split.
        mPhi: float, optional
            Mass of mediator.
        quenching: float, optional
            quenching factor, needed for experiments that can have multiple options.
    Returns:
        (vmin_min, vmin_max, vmin_step): tuple (float, float, int)
             DM mass range and number or steps
    """
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
            mx_range_options = {(0, 1000.): (5, 25, num_steps),
#            mx_range_options = {(0, 1000.): (8, 15, num_steps),
                                (0, 0.): (5.5, 25, num_steps),
                                (-30, 1000.): (2, 4, num_steps),
                                (-30, 0.): (2, 3, num_steps),
                                (-50, 1000.): (1.4, 2.5, num_steps),
                                }
        else:
            mx_range_options = {(0, 1000.): (6, 30, num_steps),
#            mx_range_options = {(0, 1000.): (10, 25, num_steps),
#                                (0, 0.): (7, 22, num_steps),
                                (0, 0.): (7, 29, num_steps),
                                (-30, 1000.): (3, 5, num_steps),
                                (-30, 0.): (2.5, 4, num_steps),
                                (-50, 1000.): (1.9, 3.5, num_steps),
                                }
    elif exper_name == "DAMA2010ISmRebinned":
        num_steps = 80
        if quenching == 0.06:
            mx_range_options = {(0, 1000.): (30, 100, num_steps),
#            mx_range_options = {(0, 1000.): (55, 85, num_steps),
#                                (0, 0.): (22, 130, num_steps),
                                (0, 0.): (22, 100, num_steps),
                                (-30, 1000.): (25, 100, num_steps),
                                (-30, 0.): (40, 100, num_steps),
                                (-50, 1000.): (25, 50, num_steps),
                                (50, 1000.): (35, 100, num_steps),
                                (50, 0.): (55, 100, num_steps),
                                (100, 1000.): (45, 100, num_steps),
                                (100, 0.): (50, 500, num_steps),
#                                (100, 0.): (50, 90, num_steps),
                                }
        else:
            mx_range_options = {(0, 1000.): (22, 80, num_steps),
#            mx_range_options = {(0, 1000.): (40, 65, num_steps),
#                                (0, 0.): (29, 90, num_steps),
                                (0, 0.): (22, 100, num_steps),
                                (-30, 1000.): (20, 80, num_steps),
                                (-30, 0.): (30, 100, num_steps),
                                (-50, 1000.): (18, 40, num_steps),
                                (50, 1000.): (30, 100, num_steps),
                                (50, 0.): (40, 100, num_steps),
                                (100, 1000.): (41, 100, num_steps),
                                (100, 0.): (30, 300, num_steps),
#                                (100, 0.): (42, 65, num_steps),
                                }
    elif exper_name == "DAMA2010NaSmRebinned_TotRateLimit":
        num_steps = 60
#        mx_range_options = {(0, 1000.): (3, 15, num_steps),
        mx_range_options = {(0, 1000.): (3, 20, num_steps),
                            (0, 0.): (3, 15, num_steps),
                            (-30, 1000.): (1, 10, num_steps),
                            (-30, 0.): (1, 10, num_steps),
                            (-50, 1000.): (1, 10, num_steps),
                            }
    else:
        num_steps = 60
        if 'SI' in scattering[0]:
            mx_range_options = {(0, 1000.): (3, 30, num_steps),
                                (-50, 1000.): (2, 10, num_steps),
                                (-200, 1000.): (0.4, 2, num_steps),
                                (-500, 1000.): (0.2, 1, num_steps)
                                }
        else:
            mx_range_options = {(0, 1000.): (3, 100, num_steps),
                                (0, 0.): (3, 130, num_steps),
                                (-30, 1000.): (1, 100, num_steps),
                                (-30, 0.): (1, 100, num_steps),
                                (-50, 1000.): (1, 50, num_steps),
                                (50, 1000.): (10, 100, num_steps),
                                (100, 1000.): (30, 100, num_steps),
                                }
    return mx_range_options[(delta, mPhi)]


""" List of input values of the form (fn, delta, mPhi).
"""
scattering = 'SI'
if scattering == 'SI':
    input_list = [(1, 0, 1000.), (1, -50, 1000.), (1, -200, 1000.), (1, -500, 1000.),  # 0 - 3
                  (-0.8, 0, 1000.), (-0.8, -50, 1000.), (-0.8, -200, 1000.), (-0.8, -500, 1000.),  # 4 - 7
                  (-0.7, 0, 1000.), (-0.7, -50, 1000.), (-0.7, -200, 1000.)]  # 8 - 10
else:
    input_list = [(-1/16.4, 0, 1000.),  # 0
                  (0, 0, 1000.), (0, -30, 1000.),  (0, -50, 1000.), (0, 0, 0.), (0, -30, 0.),  # 1 - 5
                  (0, 50, 1000.), (0, 100, 1000.), (0, 100, 0.)]  # 6 - 8
