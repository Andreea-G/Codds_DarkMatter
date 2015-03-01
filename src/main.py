# -*- coding: utf-8 -*-
"""
Created on Thu Nov 20 22:52:11 2014

@author: Andreea
"""
from __future__ import print_function
from __future__ import division
#import profile
from runprogram import *
import os   # for speaking
from input_DAMApaper import *



def main():
    implemented_exper = ["superCDMS", \
        "LUX2013zero", "LUX2013one", "LUX2013three", "LUX2013five", "LUX2013many", \
        "SIMPLEModeStage2", "PICASSO", "KIMS2012", "DAMA2010NaSmRebinned", "DAMA2010ISmRebinned", \
        "DAMA2010NaSmRebinned_TotRateLimit", 
        "XENON10", "CDMSlite2013CoGeNTQ", "CDMSSi2012"]
    scattering_type = 'SD66'
    fp = 1.

    plot_dots = F
    RUN_PROGRAM = F
    MAKE_PLOT = T

    qKIMS_list = [0.05, 0.1]
    qDAMANa_list = [0.4, 0.3]
    qDAMAI_list = [0.09, 0.06]
    qDAMANa_Rate_list = [0.4]
    quenching_list = {"KIMS2012": qKIMS_list,
                      "DAMA2010NaSmRebinned": qDAMANa_list,
                      "DAMA2010ISmRebinned": qDAMAI_list,
                      "DAMA2010NaSmRebinned_TotRateLimit": qDAMANa_Rate_list,
    }

    exper_list = [implemented_exper[i] for i in [0,1,2,3,4,5,6,7,8,9]]
#    exper_list = implemented_exper
    filename_tail_list = [""]
    plt.close()
    for exper_name in exper_list:
        for filename_tail in filename_tail_list:
            for (fn, delta, mPhi) in input_list[0:1]:
                for quenching in quenching_list.get(exper_name, [None]):
                    (mx_min, mx_max, num_steps) = DM_mass_range(exper_name, delta, mPhi, quenching)
                    print(mx_min, " ", mx_max, " ", num_steps)
                    run_program(exper_name, scattering_type, mPhi, fp, fn, delta, \
                        mx_min, mx_max, num_steps, RUN_PROGRAM, MAKE_PLOT, \
                        filename_tail, plot_dots = plot_dots, quenching = quenching)
    plt.show()
    
    if RUN_PROGRAM:
        os.system("say 'Finished running program'")
    
if __name__ == '__main__':
    main()
#    profile.run("main()")

