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


input_filename_list = {True: "input_HaloDep",
                       False: "input_HaloIndep",
                       }

EHIBools = namedtuple('EHIBools', ['ResponseTables', 'OptimalLikelihood',
                                   'ImportOptimalLikelihood',
                                   'ConstrainedOptimalLikelihood',
                                   'VminLogetaSamplingTable', 'LogLikelihoodList',
                                   'ConfidenceBand', 'ConfidenceBandPlot'])
len_EHIBools = len(EHIBools._fields)
EHIBools.__new__.__defaults__ = tuple([F] * len_EHIBools)


class Input:
    ''' Class implementing the input parameters to be passed to run_program().
    Input:
        HALO_DEP: bool, optional
            Whether the analysis is halo-dependent or halo-independent.
        implemented_exper_list: list of strings
            List of names of implemented expemeriments.
        exper_indices: slice or list of ints, optional
            Indices in the implemented_exper_list that the program will be run for.
        input_indices: slice or list of ints, optional
            Indices in the list of input values that the program will be run for.
        scattering_types: string or list of strings, optional
            Type of scattering, or list of types of scattering. Can be
                - 'SI' (spin-independent)
                - 'SDAV' (spin-independent, axial-vector)
                - 'SDPS' (spin-independent, pseudo-scalar).
        RUN_PROGRAM: bool, optional
            Whether the data should be (re-)computed.
        MAKE_REGIONS: bool, optional
            Whether the regions should be (re-)computed in the case of halo-dependent
            analysis and experiments with potential DM signals.
        MAKE_PLOT: bool, optional
            Whether the data should be plotted.
        EHI_Method: ndarray of bools, optional
            Whether each step of the EHI Method is to be performed.
        OUTPUT_MAIN_DIR: string, optional
            Name of main output directory.
        filename_tail_list: list of strings, optional
            List of tags to be added to the file name.
        extra_tail: string, optional
            Additional tail to be added to filenames for the EHI confidence band.
        plot_dots: bool, optional
            Whether the plot should show the data points or just the interpolation.
        CL_list: list, optional
            List of confidence levels.
        sigma_dev_list: list, optional
            List of how many sigma deviations away. This then gets converted into
            confidence levels.
    '''
    def __init__(self, HALO_DEP,
                 implemented_exper_list, exper_indices=slice(None),
                 input_indices=slice(None),
                 scattering_types='SI',
                 RUN_PROGRAM=False, MAKE_REGIONS=False, MAKE_PLOT=False, EHI_METHOD={},
                 OUTPUT_MAIN_DIR="../Output/", filename_tail_list=[""], extra_tail="",
                 plot_dots=False,
                 CL_list=[0.9], sigma_dev_list=[1]):
        print('HALO_DEP =', HALO_DEP)
        module = import_file(input_filename_list[HALO_DEP] + ".py")

        self.implemented_exper_list = np.array(implemented_exper_list)
        self.SetExperList(exper_indices)
        self.scattering_type_list = scattering_types \
            if isinstance(scattering_types, list) \
            else [scattering_types]
        self.filename_tail_list = filename_tail_list
        self.extra_tail = extra_tail  # for EHI method
        self.input_list = np.array(module.input_list)[input_indices]

        self.OUTPUT_MAIN_DIR = OUTPUT_MAIN_DIR
        self.MAKE_PLOT = MAKE_PLOT
        self.plot_dots = plot_dots
        self.RUN_PROGRAM = RUN_PROGRAM
        self.MAKE_REGIONS = MAKE_REGIONS
        self.HALO_DEP = HALO_DEP
        self.EHI_METHOD = EHIBools(**EHI_METHOD)

        self.qKIMS_list = [0.05, 0.1]
        self.qDAMANa_list = [0.4, 0.3]
        self.qDAMAI_list = [0.09, 0.06]
        self.qDAMANa_Rate_list = [0.4]

        self.fp = 1

        if CL_list is not None:
            self.confidence_levels = CL_list
        else:
            self.confidence_levels = []
        if sigma_dev_list is not None:
            self.confidence_levels.extend([confidence_level(s) for s in sigma_dev_list])
        self.confidence_levels.sort()

    def SetExperList(self, exper_indices):
        self.exper_list = self.implemented_exper_list[exper_indices]

    def SetScattering_type(self, scattering_type_list):
        try:
            len(scattering_type_list)
        except TypeError:
            self.scattering_type_list = [scattering_type_list]
        else:
            self.scattering_type_list = scattering_type_list

    def SetInputList(self, input_indices):
        module = import_file(input_filename_list[self.HALO_DEP] + ".py")
        self.input_list = np.array(module.input_list)[input_indices]

    def QuenchingList(self):
        quenching_list = {"KIMS2012": self.qKIMS_list,
                          "DAMA2010NaSmRebinned": self.qDAMANa_list,
                          "DAMA2010ISmRebinned": self.qDAMAI_list,
                          "DAMA2010NaSmRebinned_TotRateLimit": self.qDAMANa_Rate_list,
                          }
        q = [quenching_list.get(exp, [None]) for exp in self.exper_name.split()]
        print('q =', repr(q))
        return zip(*q)

    def _GetKwargs(self):
        ''' Collects the input parameters that are to be passed to run_program() are the
        member variables of this class, that are not ending in '_list' and are not hidden.
        Returns:
            kwargs: dictionary
                Keyward arguments that will be passed as input to run_program().
        '''
        attributes = inspect.getmembers(self, lambda a: not(inspect.isroutine(a)))
        kwargs = dict([a for a in attributes
                       if '__' not in a[0] and '_list' not in a[0]])
        return kwargs

    def _Run_HaloIndep(self):
        ''' Run the program for halo-independent analysis.
        '''
        module = import_file(input_filename_list[self.HALO_DEP] + ".py")
        for self.exper_name, self.scattering_type, self.filename_tail, \
                (self.mx, self.fn, self.delta, self.mPhi) \
                in product(self.exper_list, self.scattering_type_list,
                           self.filename_tail_list, self.input_list):
            if self.exper_name == "CDMSSi2012" and np.any(self.EHI_METHOD):
                self.vmin_EHIBand_range = \
                    module.Vmin_EHIBand_range(self.exper_name, self.mx,
                                              self.delta, self.mPhi)
                self.logeta_EHIBand_percent_range = module.logeta_EHIBand_percent_range
                self.steepness = module.Steepness(self.exper_name, self.mx,
                                                  self.delta, self.mPhi)
                self.logeta_guess = module.Logeta_guess(self.exper_name, self.mx,
                                                        self.delta, self.mPhi)
            for self.quenching in self.QuenchingList():
                self.vmin_range = \
                    module.Vmin_range(self.exper_name.split()[0], self.mx,
                                      self.delta, mPhi=self.mPhi,
                                      quenching=self.quenching[0],
                                      EHI_METHOD=np.any(self.EHI_METHOD))
                if np.any(self.EHI_METHOD):
                    self.vmin_EHIBand_range = \
                        module.Vmin_EHIBand_range(self.exper_name.split()[0], self.mx,
                                                  self.delta, self.mPhi)
                print(self.vmin_range)
                if len(self.quenching) == 1:
                    self.quenching = self.quenching[0]
                kwargs = self._GetKwargs()
                run_program = RunProgram()
                run_program(**kwargs)
        return

    def _Run_HaloDep(self):
        ''' Run the program for halo-dependent analysis.
        '''
        module = import_file(input_filename_list[self.HALO_DEP] + ".py")
        for self.exper_name, self.scattering_type, self.filename_tail, \
                (self.fn, self.delta, self.mPhi) \
                in product(self.exper_list, self.scattering_type_list,
                           self.filename_tail_list, self.input_list):
            for self.quenching in self.QuenchingList():
                try:
                    self.mx_range = \
                        module.DM_mass_range(self.exper_name.split()[0], self.delta,
                                             self.mPhi, self.quenching[0])
                except KeyError as key_error:
                    print('KeyError:', key_error)
                    continue
                print(self.mx_range)
                if len(self.quenching) == 1:
                    self.quenching = self.quenching[0]
                kwargs = self._GetKwargs()
                run_program = RunProgram()
                try:
                    run_program(**kwargs)
                except FileNotFoundError as file_error:
                    print('FileNotFoundError:', file_error)
        return

    def RunProgram(self):
        ''' Main run of the program.
        '''
        if self.HALO_DEP:
            self._Run_HaloDep()
        else:
            self._Run_HaloIndep()
        return


if __name__ == "__main__":
    inp = Input()
    inp.RunProgram()
