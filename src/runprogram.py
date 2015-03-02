# -*- coding: utf-8 -*-
"""
Created on Sat Dec  6 15:32:34 2014

@author: Andreea
"""

from __future__ import print_function
from __future__ import division
from experiment import *
from experiment_HaloIndep import *
import matplotlib.pyplot as plt

def Plot_Upper_Limit(exper_name, upper_limit, HALO_DEP, plot_dots = True, plot_close = True, plot_show = True):
    from scipy.interpolate import interp1d

    if plot_close:
        plt.close()

    if upper_limit.size == 0:
        print("upper_limit is empty!")
    elif upper_limit.ndim == 1:
        x = [upper_limit[0]]
        y = [upper_limit[1]]
        plt.plot(x, y, "o")
    else:
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

    if HALO_DEP:
        plt.xlabel('Log10(m [GeV])')
        plt.ylabel('Log10(sigma)')
    else:
        plt.xlabel('vmin [km/s]')
        plt.ylabel('eta rho sigma / m [days^-1]')

    if plot_show:
        plt.show()


def run_program(exper_name, scattering_type, mPhi, fp, fn, delta, \
    RUN_PROGRAM, MAKE_PLOT, HALO_DEP, \
    mx = None, mx_range = None, vmin_range = None, \
    filename_tail = "", OUTPUT_MAIN_DIR = "Output/", plot_dots = True, quenching = None):

    print('name = ', exper_name)
    if HALO_DEP:
        print('Halo Dependent')
        if exper_name in MaximumGapLimit_exper:
            exper = MaxGapExperiment(exper_name, scattering_type, mPhi)
        elif exper_name in GaussianLimit_exper:
            exper = GaussianExperiment(exper_name, scattering_type, mPhi, quenching)
        elif exper_name in Poisson_Exper:
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
        if exper_name in MaximumGapLimit_exper:
            exper = MaxGapExperiment_HaloIndep(exper_name, scattering_type, mPhi)
        else:
            print("Error! This experiment was not implemented!")
            return

    output_file_no_extension = Output_file_name(exper_name, scattering_type, mPhi, mx, fp, fn, delta, \
        HALO_DEP, filename_tail, OUTPUT_MAIN_DIR, quenching)

    if RUN_PROGRAM:
        output_file = output_file_no_extension + "_temp.dat" 
        f_handle = open(output_file, 'w')   # clear the file first
        f_handle.close()

        if HALO_DEP:
            (mx_min, mx_max, num_steps) = mx_range
            upper_limit = exper.UpperLimit(fp, fn, delta, mx_min, mx_max, num_steps, \
                output_file)
        else:
            (vmin_start, vmin_end, vmin_step) = vmin_range
            upper_limit = exper.UpperLimit(mx, fp, fn, delta, vmin_start, vmin_end, vmin_step, \
                output_file)
        print("upper_limit = ", upper_limit)
        print("diff response calls = " , exper.count_diffresponse_calls)
        print("response calls = " , exper.count_response_calls)
        output_file = output_file_no_extension + ".dat"
        print(output_file)
        np.savetxt(output_file, upper_limit)

    if MAKE_PLOT:
        output_file = output_file_no_extension + ".dat"
        upper_limit = np.loadtxt(output_file)
        print("upper_limit = ", upper_limit)
        Plot_Upper_Limit(exper_name, upper_limit, HALO_DEP, plot_dots, plot_close = False, \
            plot_show = False)

