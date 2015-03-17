# -*- coding: utf-8 -*-
"""
Created on Thu Nov 20 22:52:11 2014

@author: Andreea
"""
from __future__ import print_function
from __future__ import division
import profile
from runprogram import *
import os   # for speaking
#from input_DAMApaper import *
from input_HaloIndep import *


def main():
    implemented_exper = ["superCDMS", \
        "LUX2013zero", "LUX2013one", "LUX2013three", "LUX2013five", "LUX2013many", \
        "SIMPLEModeStage2", "PICASSO", "KIMS2012", "DAMA2010NaSmRebinned", "DAMA2010ISmRebinned", \
        "DAMA2010NaSmRebinned_TotRateLimit", 
        "XENON10", "CDMSlite2013CoGeNTQ", "CDMSSi2012"]
    scattering_type = 'SI'
    fp = 1.

    plot_dots = T
    RUN_PROGRAM = T
    MAKE_PLOT = F
    HALO_DEP = F
    HALO_INDEP = not HALO_DEP
#    FOX_METHOD = [T, T, T, F, F, F, F, F]     # Multiple
#    FOX_METHOD = [T, F, F, F, F, F, F, F]     # ResponseTables
    FOX_METHOD = [F, T, F, F, F, F, F, F]     # OptimalLikelihood
#    FOX_METHOD = [F, F, T, F, F, F, F, F]     # ImportOptimalLikelihood
#    FOX_METHOD = [F, F, F, T, F, F, F, F]     # ConstrainedOptimalLikelihood
#    FOX_METHOD = [F, F, F, F, T, F, F, F]     # VminLogetaSamplingTable
#    FOX_METHOD = [F, F, F, F, F, T, F, F]     # LogLikelihoodList
#    FOX_METHOD = [F, F, F, F, F, F, T, F]     # FoxBand
#    FOX_METHOD = [F, F, F, F, F, F, F, T]     # FoxBandPlot
    
    qKIMS_list = [0.05, 0.1]
    qDAMANa_list = [0.4, 0.3]
    qDAMAI_list = [0.09, 0.06]
    qDAMANa_Rate_list = [0.4]
    quenching_list = {"KIMS2012": qKIMS_list,
                      "DAMA2010NaSmRebinned": qDAMANa_list,
                      "DAMA2010ISmRebinned": qDAMAI_list,
                      "DAMA2010NaSmRebinned_TotRateLimit": qDAMANa_Rate_list,
    }

    exper_list = [implemented_exper[i] for i in [14]]
#    exper_list = implemented_exper
    filename_tail_list = ["_jump_smooth"]
    OUTPUT_MAIN_DIR = "../Output1/"

    try:
        plt.show()
        
        if HALO_INDEP:
            for exper_name in exper_list:
                for filename_tail in filename_tail_list:
                    for (mx, fn, delta, mPhi) in input_list[0:1]:
                        for quenching in quenching_list.get(exper_name, [None]):
                            (vmin_start, vmin_end, vmin_step) = Vmin_range(exper_name, mx, delta, mPhi, quenching)
                            print(vmin_start, " ", vmin_end, " ", vmin_step)
                            run_program(exper_name, scattering_type, mPhi, fp, fn, delta, \
                                RUN_PROGRAM, MAKE_PLOT, HALO_DEP, FOX_METHOD, \
                                mx = mx, vmin_range = (vmin_start, vmin_end, vmin_step), \
                                vmin_FoxBand_range = vmin_FoxBand_range, \
                                logeta_FoxBand_percent_range = logeta_FoxBand_percent_range, \
                                steepness = Steepness(exper_name, mx, delta, mPhi), \
                                logeta_guess = Logeta_guess(exper_name, mx, delta, mPhi), \
                                filename_tail = filename_tail, OUTPUT_MAIN_DIR = OUTPUT_MAIN_DIR, \
                                plot_dots = plot_dots, quenching = quenching)

        if HALO_DEP:
            for exper_name in exper_list:
                for filename_tail in filename_tail_list:
                    for (fn, delta, mPhi) in input_list[0:1]:
                        for quenching in quenching_list.get(exper_name, [None]):
                            (mx_min, mx_max, num_steps) = DM_mass_range(exper_name, delta, mPhi, quenching)
                            print(mx_min, " ", mx_max, " ", num_steps)
                            run_program(exper_name, scattering_type, mPhi, fp, fn, delta, \
                                RUN_PROGRAM, MAKE_PLOT, HALO_DEP, FOX_METHOD, \
                                mx_range = (mx_min, mx_max, num_steps), \
                                filename_tail = filename_tail, OUTPUT_MAIN_DIR = OUTPUT_MAIN_DIR, \
                                plot_dots = plot_dots, quenching = quenching)
                                
        plt.show()
        
    finally:
        if RUN_PROGRAM:
#            os.system("say 'Finished running program'")
            None
    
if __name__ == '__main__':
    main()
#    profile.run("main()")

