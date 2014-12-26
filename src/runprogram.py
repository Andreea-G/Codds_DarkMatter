# -*- coding: utf-8 -*-
"""
Created on Sat Dec  6 15:32:34 2014

@author: Andreea
"""

from __future__ import print_function
from __future__ import division
from experiment import *
import matplotlib.pyplot as plt

def Plot_Upper_Limit(exper_name, upper_limit, plot_dots = True, plot_close = True, plot_show = True):
    from scipy.interpolate import interp1d

    if plot_close:
        plt.close()

    if upper_limit.size == 0:
        print("upper_limit is empty!")
    elif upper_limit.ndim == 1:
        if exper_name in MaximumGapLimit_exper or exper_name in Poisson_Exper:
            x = [upper_limit[0]]
            y = [-np.log10(conversion_factor) + upper_limit[0] + upper_limit[1]]
        else:
            x = np.log10([upper_limit[0]])
            y = np.log10(upper_limit[1])
        plt.plot(x, y, "o")
    else:
        if exper_name in MaximumGapLimit_exper or exper_name in Poisson_Exper:
            x = upper_limit[:,0]
            y = -np.log10(conversion_factor) + upper_limit[:,0] + upper_limit[:,1]
        else:
            x = np.log10(upper_limit[:,0])
            y = np.log10(upper_limit[:,1])
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
    plt.xlabel('Log10(m [GeV])')
    plt.ylabel('Log10(sigma)')
    if plot_show:
        plt.show()


def run_program(exper_name, scattering_type, mPhi, fp, fn, delta, mx_min, mx_max, \
    num_steps, RUN_PROGRAM, MAKE_PLOT, filename_tail = "", plot_dots = True, quenching = None):
    print('name = ', exper_name)
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

    output_dir = OutputDirectory(OUTPUT_MAIN_DIR, scattering_type, mPhi, delta)
    if exper_name in DAMARegion_exper or exper_name in DAMALimit_exper:
        output_file_no_extension = "./" + output_dir + "pbesebTab_" + exper.name \
            + "_mxsigma"
    else:
        output_file_no_extension = "./" + output_dir + "UpperLimitSHM_" + exper.name \
            + "_mxsigma"

    if vesc != default_vesc:
        output_file_no_extension += "_vesc" \
            + str(math.trunc(round(vesc)))
    if vobs != default_vobs:
        output_file_no_extension += "_vobs" \
            + str(math.trunc(round(vobs)))
    output_file_no_extension = output_file_no_extension + FileNameTail(fp, fn, mPhi)

    output_file_no_extension += filename_tail
    if quenching != None:
        output_file_no_extension += "_q" + str(quenching)
    print(output_file_no_extension)

    if RUN_PROGRAM:
        output_file = output_file_no_extension + "_temp.dat" 
        f_handle = open(output_file, 'w')   # clear the file first
        f_handle.close()

        upper_limit = exper.UpperLimit(fp, fn, delta, mx_min, mx_max, num_steps, \
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
        Plot_Upper_Limit(exper_name, upper_limit, plot_dots, plot_close = False, \
            plot_show = False)

