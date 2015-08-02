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
        ["SuperCDMS",  # 0
         "LUX2013zero", "LUX2013one", "LUX2013three", "LUX2013five", "LUX2013many",  # 1 - 5
         "SIMPLEModeStage2", "PICASSO", "KIMS2012", "XENON10", "XENON100",  # 6 - 10
         "DAMA2010NaSmRebinned", "DAMA2010ISmRebinned", "DAMA2010NaSmRebinned_TotRateLimit",  # 11 - 13
         "DAMA2010NaSmRebinned DAMA2010ISmRebinned", "DAMA2010ISmRebinned DAMA2010NaSmRebinned",  # 14 - 15
         "CDMSlite2013CoGeNTQ", "CDMSSi2012", "CDMSSiGeArtif", "CDMSSiArtif"]  # 16 - 19

    EHI_METHOD = {}
    # EHI_METHOD['ResponseTables'] = T
    # EHI_METHOD['OptimalLikelihood'] = T
    # EHI_METHOD['ImportOptimalLikelihood'] = T
    # EHI_METHOD['ConstrainedOptimalLikelihood'] = T
    # EHI_METHOD['VminLogetaSamplingTable'] = T
    # EHI_METHOD['LogLikelihoodList'] = T
    # EHI_METHOD['ConfidenceBand'] = T
    # EHI_METHOD['ConfidenceBandPlot'] = T

    HALO_DEP = F
    plot_dots = F
    RUN_PROGRAM = F
    MAKE_REGIONS = F
    MAKE_PLOT = F

    inp = Input(HALO_DEP, implemented_exper_list,
                RUN_PROGRAM=RUN_PROGRAM, MAKE_REGIONS=MAKE_REGIONS, MAKE_PLOT=MAKE_PLOT,
                EHI_METHOD=EHI_METHOD, plot_dots=plot_dots)

    inp.SetScattering_type(['SI'])
    inp.SetInputList([0])
    inp.SetExperList([17])
    inp.OUTPUT_MAIN_DIR = "../Output_Band/"
    inp.filename_tail_list = [""]
    inp.extra_tail = "_mix"
    # inp.initial_energy_bin = [3, 6]
    # inp.confidence_levels.extend([confidence_level(s) for s in [3, 5, 7]])

    try:
        plt.close()
        inp.RunProgram()
        if MAKE_PLOT or EHI_METHOD.get('ConfidenceBandPlot', F):
            # plt.ylim([-32, -20])
            # plt.ylim([-28, -13])
            plt.show()
    finally:
        if inp.RUN_PROGRAM or inp.MAKE_REGIONS:
            os.system("say 'Finished running program'")
            # pass


if __name__ == '__main__':
    main()
    # profile.run("main()")
