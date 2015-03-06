# -*- coding: utf-8 -*-
"""
Created on Sat Dec  6 15:32:34 2014

@author: Andreea
"""

from __future__ import print_function
from __future__ import division
from experiment import *
from experiment_HaloIndep import *
from experiment_FoxMethod import *
import matplotlib.pyplot as plt

def Plot_Upper_Limit(exper_name, upper_limit, HALO_DEP, plot_dots = True, plot_close = True, plot_show = True):
    ''' Make plots for the upper limits.
        Input:
            exper_name = name of experiment
            upper_limit = list of x and y coordinates for points represnting the upper limit on the plot
            HALO_DEP = True or False, whether the analysis is halo-dependent or halo-independent
            plot_dots = True or False, whether the plot should show the data points or just the interpolation
            plot_close = True or False, whether the plot should be cleared and started from scratch,
                or new limits should be added to a previous plot
            plot_show = True or False, whether the plot should be shown or not.

    '''
    from scipy.interpolate import interp1d

    if plot_close:
        plt.close()

    # make a list of the x and y coordinates of the plots, and plot them
    if upper_limit.size == 0:   # nothing to plot
        print("upper_limit is empty!")
    elif upper_limit.ndim == 1: # only one point, so no interpolation
        x = [upper_limit[0]]
        y = [upper_limit[1]]
        plt.plot(x, y, "o")
    else:   # more than one point, so decide on the interpolation order and plot
        x = upper_limit[:,0]
        y = upper_limit[:,1]
        num_points = x.size
        if num_points == 2:
            interp_kind = "linear"
        elif num_points == 3:
            interp_kind = "quadratic"
        else:
            interp_kind = "cubic"
        interp = interp1d(x, y, kind = interp_kind)
        x1 = np.linspace(x[0], x[-1], 50)
        if plot_dots:
            plt.plot(x, y, "o")
        plt.plot(x1, interp(x1))

    # set axis labels, depending on whether it is for halo-dependent or not
    if HALO_DEP:
        plt.xlabel('Log10(m [GeV])')
        plt.ylabel('Log10(sigma)')
    else:
        plt.xlabel('vmin [km/s]')
        plt.ylabel('eta rho sigma / m [days^-1]')

    # show plot
    if plot_show:
        plt.show()


def run_program(exper_name, scattering_type, mPhi, fp, fn, delta, \
    RUN_PROGRAM, MAKE_PLOT, HALO_DEP, FOX_METHOD, \
    mx = None, mx_range = None, vmin_range = None, \
    vmin_FoxBand_range = None, logeta_FoxBand_percent_range = None, steepness = None, \
    filename_tail = "", OUTPUT_MAIN_DIR = "Output/", plot_dots = True, quenching = None):
    ''' Main run of the program.
        Input:
            exper_name = name of experiment
            scattering_type: 'SI' for spin-dependent, 'SDPS' for pseudo-scalar, 'SDAV' for avial-vector
            mPhi = mass of mediator
            fp and fn = couplings to proton and neutron
            delta = DM mass split
            RUN_PROGRAM = True or False, whether the data should be (re-)computed
            MAKE_PLOT = True or False, whether the data should be plotted
            HALO_DEP = True or False, whether the analysis is halo-dependent or halo-independent
            FOX_METHOD = array of True or False, whether each step of the Fox Method is to be performed
            mx = DM mass, only give for halo-independent analysis
            mx_range = (mx_min, mx_max, num_steps) = DM mass range and number or steps, only for halo-dependent
            vmin_range = (vmin_min, vmin_max, vmin_step) = vmin range and step size, only for halo-independent
            filename_tail = optional tag to be added to the file name 
            OUTPUT_MAIN_DIR = name of main output directory, if different from "Output/"
            plot_dots = True or False, whether the plot should show the data points or just the interpolation
            quenching = quenching factor, needed for experiments that can have multiple options (such as KIMS or DAMA).
    '''

    print('name = ', exper_name)
    # Select which experiment class we must use, depending on what statistical analysis we need.
    if HALO_DEP:
        print('Halo Dependent')
        if exper_name in MaximumGapLimit_exper:
            exper = MaxGapExperiment(exper_name, scattering_type, mPhi)
        elif exper_name in GaussianLimit_exper:
            exper = GaussianExperiment(exper_name, scattering_type, mPhi, quenching)
        elif exper_name in Poisson_exper:
            exper = PoissonExperiment(exper_name, scattering_type, mPhi)
        elif exper_name in DAMARegion_exper:
            exper = DAMAExperiment(exper_name, scattering_type, mPhi, quenching)
        elif exper_name in DAMALimit_exper:
            exper = DAMATotalRateExperiment(exper_name, scattering_type, mPhi, quenching)
        else:
            print("Error! This experiment was not implemented!")
            return
    else:
        print('Halo Independent')
        if exper_name in FoxMethod_exper and np.any(FOX_METHOD):
            print('Fox Method')
            exper = Experiment_FoxMethod(exper_name, scattering_type, mPhi)
        elif exper_name in MaximumGapLimit_exper:
            exper = MaxGapExperiment_HaloIndep(exper_name, scattering_type, mPhi)
        else:
            print("Error! This experiment was not implemented!")
            return

    # get the file name specific to the parameters used for this run
    output_file_no_extension = Output_file_name(exper_name, scattering_type, mPhi, mx, fp, fn, delta, \
        HALO_DEP, filename_tail, OUTPUT_MAIN_DIR, quenching)

    # (re-)compute the data
    if RUN_PROGRAM:
        output_file = output_file_no_extension + "_temp.dat" 
        f_handle = open(output_file, 'w')   # clear the file first
        f_handle.close()

        # calculate the upper limit, for halo-dependent or halo-independent
        if HALO_DEP:
            (mx_min, mx_max, num_steps) = mx_range
            upper_limit = exper.UpperLimit(fp, fn, delta, mx_min, mx_max, num_steps, \
                output_file)
        else:
            (vmin_min, vmin_max, vmin_step) = vmin_range
            if not np.any(FOX_METHOD):
                upper_limit = exper.UpperLimit(mx, fp, fn, delta, vmin_min, vmin_max, vmin_step, \
                    output_file)
            else:
                if FOX_METHOD[0]:
                    exper.ResponseTables(vmin_min, vmin_max, vmin_step, mx, fp, fn, delta, output_file_no_extension)
                if FOX_METHOD[1]:
                    exper.OptimalLikelihood(output_file_no_extension)
                if FOX_METHOD[2]:
                    exper.ImportOptimalLikelihood(output_file_no_extension)
                    exper.PlotOptimum()
                if FOX_METHOD[3]:
                    vminStar = 612.853861611
                    logetaStar = -25.5
                    exper.ImportOptimalLikelihood(output_file_no_extension)
                    exper.ConstrainedOptimalLikelihood(vminStar, logetaStar, output_file_no_extension, plot = True)
                if np.any(FOX_METHOD[4:]):
                    (vmin_band_min, vmin_band_max, vmin_num_steps) = vmin_FoxBand_range
                    (logeta_percent_minus, logeta_percent_plus, logeta_num_steps) = logeta_FoxBand_percent_range
                    if steepness != None:
                        (steepness_vmin, steepness_vmin_center, steepness_logeta) = steepness
                        exper.VminSamplingList(output_file_no_extension, \
                            vmin_band_min, vmin_band_max, vmin_num_steps, steepness_vmin, steepness_vmin_center)
                        exper.VminLogetaSamplingTable(output_file_no_extension, \
                            logeta_percent_minus, logeta_percent_plus, logeta_num_steps, steepness_logeta)
                    else:
                        exper.VminSamplingList(output_file_no_extension, \
                            vmin_band_min, vmin_band_max, vmin_num_steps, plot = not np.any(FOX_METHOD[5:]))
                        exper.VminLogetaSamplingTable(output_file_no_extension, \
                            logeta_percent_minus, logeta_percent_plus, logeta_num_steps, plot = not np.any(FOX_METHOD[5:]))
                if FOX_METHOD[5]:
                    print("vmin_FoxBand_range = ", vmin_band_min, " ", vmin_band_max, " ", vmin_num_steps)
                    print("logeta_FoxBand_percent_range = ", logeta_percent_minus, " ", logeta_percent_plus, " ", logeta_num_steps)                    
                    exper.LogLikelihoodList(output_file_no_extension)
                if FOX_METHOD[6]:
                    exper.ImportOptimalLikelihood(output_file_no_extension)
                    delta_logL = 1
                    interpolation_order = 1
                    exper.FoxBand(output_file_no_extension, delta_logL, interpolation_order)
                
        if HALO_DEP or not np.any(FOX_METHOD):
            print("upper_limit = ", upper_limit)
            print("diff response calls = " , exper.count_diffresponse_calls)
            print("response calls = " , exper.count_response_calls)
            output_file = output_file_no_extension + ".dat"
            print(output_file)  # write to file
            np.savetxt(output_file, upper_limit)

    # produce plot
    if MAKE_PLOT and not np.any(FOX_METHOD):
        output_file = output_file_no_extension + ".dat"
        upper_limit = np.loadtxt(output_file)
        print("upper_limit = ", upper_limit)
        Plot_Upper_Limit(exper_name, upper_limit, HALO_DEP, plot_dots, plot_close = False, \
            plot_show = False)

