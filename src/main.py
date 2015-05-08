# -*- coding: utf-8 -*-
"""
Created on Thu Nov 20 22:52:11 2014

@author: Andreea
"""
import profile
import os   # for speaking
from input_main import *


def main():
    implemented_exper_list = \
        ["superCDMS",
         "LUX2013zero", "LUX2013one", "LUX2013three", "LUX2013five", "LUX2013many",
         "SIMPLEModeStage2", "PICASSO", "KIMS2012", "DAMA2010NaSmRebinned", "DAMA2010ISmRebinned",
         "DAMA2010NaSmRebinned_TotRateLimit",
         "XENON10", "XENON100", "CDMSlite2013CoGeNTQ", "CDMSSi2012"]

    FOX_METHOD = {}
#    FOX_METHOD['ResponseTables'] = T
#    FOX_METHOD['OptimalLikelihood'] = T
#    FOX_METHOD['ImportOptimalLikelihood'] = T
#    FOX_METHOD['ConstrainedOptimalLikelihood'] = T
#    FOX_METHOD['VminLogetaSamplingTable'] = T
#    FOX_METHOD['LogLikelihoodList'] = T
#    FOX_METHOD['FoxBand'] = T
#    FOX_METHOD['FoxBandPlot'] = T

    HALO_DEP = F
    plot_dots = T
    RUN_PROGRAM = F
    MAKE_PLOT = T

    inp = Input(HALO_DEP, implemented_exper_list=implemented_exper_list,
                RUN_PROGRAM=RUN_PROGRAM, FOX_METHOD=FOX_METHOD,
                MAKE_PLOT=MAKE_PLOT, plot_dots=plot_dots)

    inp.SetScattering_type(['SI', 'SDAV', 'SDPS'])
    inp.SetExperList([0])
    inp.SetInputList(slice(0, 1))
    inp.OUTPUT_MAIN_DIR = "../Output1/"
    inp.filename_tail_list = ["_delete"]
    inp.extra_tail = "_mix"

    try:
        plt.close()
        inp.RunProgram()
        if MAKE_PLOT:
            plt.show()
    finally:
        if inp.RUN_PROGRAM:
            os.system("say 'Finished running program'")
#            None


if __name__ == '__main__':
    main()
#    profile.run("main()")
