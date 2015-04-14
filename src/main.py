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
#    inp.FOX_METHOD = [F] * 8
#    inp.FOX_METHOD = [T, T, T, F, F, F, F, F]     # Multiple
#    inp.FOX_METHOD = [T, F, F, F, F, F, F, F]     # ResponseTables
#    inp.FOX_METHOD = [F, T, F, F, F, F, F, F]     # OptimalLikelihood
#    inp.FOX_METHOD = [F, F, T, F, F, F, F, F]     # ImportOptimalLikelihood
#    inp.FOX_METHOD = [F, F, F, T, F, F, F, F]     # ConstrainedOptimalLikelihood
#    inp.FOX_METHOD = [F, F, F, F, T, F, F, F]     # VminLogetaSamplingTable
#    inp.FOX_METHOD = [F, F, F, F, F, T, F, F]     # LogLikelihoodList
#    inp.FOX_METHOD = [F, F, F, F, F, F, T, F]     # FoxBand
#    inp.FOX_METHOD = [F, F, F, F, T, F, F, T]     # FoxBandPlot
    HALO_DEP = T
    plot_dots = T
    RUN_PROGRAM = F
    MAKE_PLOT = T

    inp = Input(HALO_DEP, implemented_exper_list=implemented_exper_list,
                RUN_PROGRAM=RUN_PROGRAM, MAKE_PLOT=MAKE_PLOT, plot_dots=plot_dots)

    inp.SetScattering_type(['SI', 'SDAV', 'SDPS'])
    inp.SetExperList([0])
    inp.SetInputList(slice(0, 1))
    inp.OUTPUT_MAIN_DIR = "../Output1/"
    inp.filename_tail_list = ["_delete"]
    inp.extra_tail = "_mix"

    try:
        plt.close()
        inp.RunProgram()
        plt.show()
    finally:
        if inp.RUN_PROGRAM:
            os.system("say 'Finished running program'")
#            None


if __name__ == '__main__':
    main()
#    profile.run("main()")
