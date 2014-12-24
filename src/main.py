# -*- coding: utf-8 -*-
"""
Created on Thu Nov 20 22:52:11 2014

@author: Andreea
"""
from __future__ import print_function
from __future__ import division
#import profile
from runprogram import *

def main():
    implemented_exper = ["superCDMS", \
        "LUX2013zero", "LUX2013one", "LUX2013three", "LUX2013five", "LUX2013many", \
        "KIMS2012", "PICASSO", "DAMA2010NaSmRebinned", "DAMA2010ISmRebinned", "SIMPLEModeStage2", \
        "XENON10", "CDMSlite2013CoGeNTQ", "CDMSSi2012"]
    scattering_type = 'SD66'
    mPhi = 1000.
    fp = 1.
    fn = 0.
    delta = 0.

    mx_min = 3.18
    mx_max = 100.
    num_steps = 60

    plot_dots = F

    inputs = [(0, -30, 1., 100), (0, -30, 2, 5), (0, -30, 20, 60), (0, -30, 7.8, 100), \
        (0, 100, 40, 65), (0, 100, 55, 80), (0, 100, 20, 100)]

    inputs = [(0, 0, 3, 130), (0, 0, 10, 30), (0, 0, 70, 130), (0, 0, 10, 130), (0, 0, 5.6, 130), \
        (0, -30, 2, 3.5), (0, -30, 30, 70), (0, -30, 1.5, 100), (0, -30, 7, 100), (0, -30, 4, 100), \
        (0, 100, 55, 300), (0, 100, 100, 500), (0, 100, 40, 300)]
    inputs = [(0, 0, 3, 100), (0, -30, 2.6, 100)]

    RUN_PROGRAM = T
    MAKE_PLOT = F

    qKIMS_list = [0.05, 0.1]
    qDAMANa_list = [0.4, 0.3]
    qDAMAI_list = [0.09, 0.06]
    exper_list = [implemented_exper[i] for i in [0]]
    filename_tail_list = [""]
    plt.close()
    for exper_name in exper_list:
        for filename_tail in filename_tail_list:
            for (fn, delta, mx_min, mx_max) in inputs[0:1]:
                if exper_name == "KIMS2012":
                    for quenching in qKIMS_list:
                        run_program(exper_name, scattering_type, mPhi, fp, fn, delta, \
                            mx_min, mx_max, num_steps, RUN_PROGRAM, MAKE_PLOT, \
                            filename_tail, plot_dots = plot_dots, quenching = quenching)
                elif exper_name == "DAMA2010NaSmRebinned":
                    for quenching in qDAMANa_list:
                        run_program(exper_name, scattering_type, mPhi, fp, fn, delta, \
                            mx_min, mx_max, num_steps, RUN_PROGRAM, MAKE_PLOT, \
                            filename_tail, plot_dots = plot_dots, quenching = quenching)
                elif exper_name == "DAMA2010ISmRebinned":
                    for quenching in qDAMAI_list:
                        run_program(exper_name, scattering_type, mPhi, fp, fn, delta, \
                            mx_min, mx_max, num_steps, RUN_PROGRAM, MAKE_PLOT, \
                            filename_tail, plot_dots = plot_dots, quenching = quenching)
                else:
                    run_program(exper_name, scattering_type, mPhi, fp, fn, delta, \
                        mx_min, mx_max, num_steps, RUN_PROGRAM, MAKE_PLOT, \
                        filename_tail, plot_dots = plot_dots)

    plt.show()
    
if __name__ == '__main__':
    main()
#    profile.run("main()")

