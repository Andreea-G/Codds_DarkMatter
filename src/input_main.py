# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 03:40:17 2015

@author: Andreea
"""

import numpy as np
import inspect
from itertools import product
from runprogram import *
from collections import namedtuple


input_filename_list = {True: "input_DAMApaper",
                       False: "input_HaloIndep",
                       }

FoxBools = namedtuple('FoxBools', ['ResponseTables', 'OptimalLikelihood',
                                   'ImportOptimalLikelihood',
                                   'ConstrainedOptimalLikelihood',
                                   'VminLogetaSamplingTable', 'LogLikelihoodList',
                                   'FoxBand', 'FoxBandPlot'])
len_FoxBools = len(FoxBools._fields)
FoxBools.__new__.__defaults__ = tuple([F] * len_FoxBools)


class Input:
    def __init__(self, HALO_DEP,
                 implemented_exper_list=None, index_list=[0], input_slice=slice(None),
                 scattering_type='SI', scattering_type_list=None,
                 filename_tail_list=[""], extra_tail="", OUTPUT_MAIN_DIR="../Output/",
                 MAKE_PLOT=False, RUN_PROGRAM=False, FOX_METHOD={},
                 plot_dots=False,
                 delta_logL=[1, 4]):
        print(HALO_DEP)
        print(input_filename_list)
        module = import_file(input_filename_list[HALO_DEP] + ".py")

        if implemented_exper_list is not None:
            self.implemented_exper_list = np.array(implemented_exper_list)
        else:
            self.implemented_exper_list = \
                np.array(["superCDMS",
                          "LUX2013zero", "LUX2013one", "LUX2013three", "LUX2013five", "LUX2013many",
                          "SIMPLEModeStage2", "PICASSO", "KIMS2012", "DAMA2010NaSmRebinned", "DAMA2010ISmRebinned",
                          "DAMA2010NaSmRebinned_TotRateLimit",
                          "XENON10", "XENON100", "CDMSlite2013CoGeNTQ", "CDMSSi2012"])
        self.SetExperList(index_list)
        self.scattering_type_list = scattering_type_list \
            if scattering_type_list is not None \
            else [scattering_type]
        self.filename_tail_list = filename_tail_list
        self.extra_tail = extra_tail  # for Fox method
        self.input_list = np.array(module.input_list)[input_slice]

        self.OUTPUT_MAIN_DIR = OUTPUT_MAIN_DIR
        self.MAKE_PLOT = MAKE_PLOT
        self.plot_dots = plot_dots
        self.RUN_PROGRAM = RUN_PROGRAM
        self.HALO_DEP = HALO_DEP
        self.FOX_METHOD = FoxBools(**FOX_METHOD)

        self.qKIMS_list = [0.05, 0.1]
        self.qDAMANa_list = [0.4, 0.3]
        self.qDAMAI_list = [0.09, 0.06]
        self.qDAMANa_Rate_list = [0.4]

        self.fp = 1

        self.delta_logL = delta_logL  # for Fox method

    def SetExperList(self, index_list):
        self.exper_list = self.implemented_exper_list[index_list]

    def SetScattering_type(self, scattering_type_list):
        try:
            len(scattering_type_list)
        except TypeError:
            self.scattering_type_list = [scattering_type_list]
        else:
            self.scattering_type_list = scattering_type_list

    def SetInputList(self, input_slice):
        module = import_file(input_filename_list[self.HALO_DEP] + ".py")
        self.input_list = np.array(module.input_list)[input_slice]

    def QuenchingList(self):
        quenching_list = {"KIMS2012": self.qKIMS_list,
                          "DAMA2010NaSmRebinned": self.qDAMANa_list,
                          "DAMA2010ISmRebinned": self.qDAMAI_list,
                          "DAMA2010NaSmRebinned_TotRateLimit": self.qDAMANa_Rate_list,
                          }
        return quenching_list.get(self.exper_name, [None])

    def _GetKwargs(self):
        attributes = inspect.getmembers(self, lambda a: not(inspect.isroutine(a)))
        kwargs = dict([a for a in attributes
                       if '__' not in a[0] and '_list' not in a[0]])
#        print(kwargs)
        return kwargs

    def Run_HaloIndep(self):
        module = import_file(input_filename_list[self.HALO_DEP] + ".py")
        for self.exper_name, self.scattering_type, self.filename_tail, \
                (self.mx, self.fn, self.delta, self.mPhi) \
                in product(self.exper_list, self.scattering_type_list,
                           self.filename_tail_list, self.input_list):
            for self.quenching in self.QuenchingList():
                self.vmin_range = \
                    module.Vmin_range(self.exper_name, self.mx, self.delta, mPhi=self.mPhi,
                                      quenching=self.quenching,
                                      FOX_METHOD=np.any(self.FOX_METHOD))
                if np.any(self.FOX_METHOD):
                    self.vmin_FoxBand_range = \
                        module.Vmin_FoxBand_range(self.exper_name, self.mx, self.delta,
                                                  self.mPhi)
                print(self.vmin_range)
                kwargs = self._GetKwargs()
                run_program(**kwargs)
        return

    def Run_HaloDep(self):
        module = import_file(input_filename_list[self.HALO_DEP] + ".py")
        for self.exper_name, self.scattering_type, self.filename_tail, \
                (self.fn, self.delta, self.mPhi) \
                in product(self.exper_list, self.scattering_type_list,
                           self.filename_tail_list, self.input_list):
            for self.quenching in self.QuenchingList():
                self.mx_range = \
                    module.DM_mass_range(self.exper_name, self.delta, self.mPhi,
                                         self.quenching)
                print(self.mx_range)
                kwargs = self._GetKwargs()
                run_program(**kwargs)
        return

    def RunProgram(self):
        if self.HALO_DEP:
            self.Run_HaloDep()
        else:
            self.Run_HaloIndep()
        return


if __name__ == "__main__":
    inp = Input()
    inp.RunProgram()
