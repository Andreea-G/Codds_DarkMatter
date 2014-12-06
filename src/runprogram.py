# -*- coding: utf-8 -*-
"""
Created on Sat Dec  6 15:32:34 2014

@author: Andreea
"""

from __future__ import print_function
from __future__ import division
from experiment import *
import matplotlib.pyplot as plt
     
def Plot_Upper_Limit(max_gap, plot_dots = True, plot_close = True, plot_show = True):
    from scipy.interpolate import interp1d
    
    if plot_close:
        plt.close()

    if max_gap.size == 0:
        print("max_gap is empty!")
    elif max_gap.ndim == 1:
        x = [max_gap[0]]
        y = [-np.log10(3.*SpeedOfLight**2*1e4*3600*24) + max_gap[0] + max_gap[1]]
        plt.plot(x, y, "o")
    else:
        x = max_gap[:,0]
        y = -np.log10(3.*SpeedOfLight**2*1e4*3600*24) + max_gap[:,0] + max_gap[:,1]
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


def run_program(exper_name, scattering_type, mPhi, fp, fn, delta, mx_min, mx_max, num_steps, \
    RUN_PROGRAM, MAKE_PLOT, filename_tail = "", plot_dots = True):
    exper = Experiment(exper_name, scattering_type, mPhi)
    print('name = ', exper.name)
    output_dir = OutputDirectory(OUTPUT_MAIN_DIR, scattering_type, mPhi, delta)
    output_file_no_extension = "./" + output_dir + "UpperLimitSHM_" + exper.name \
        + "_mxsigma"
    if vesc != default_vesc:
        output_file_no_extension += "_vesc" \
            + str(math.trunc(round(vesc)))
    if vobs != default_vobs:
        output_file_no_extension += "_vobs" \
            + str(math.trunc(round(vobs)))
    output_file_no_extension = output_file_no_extension + FileNameTail(fp, fn)
    output_file_no_extension += filename_tail
    print(output_file_no_extension)

    if RUN_PROGRAM:  
        output_file = output_file_no_extension + "_temp.dat" 
        f_handle = open(output_file, 'w')   # clear the file first
        f_handle.close()
        
        max_gap = exper.MaximumGapLimit(fp, fn, delta, mx_min, mx_max, num_steps, output_file)
        print("max gap = ", max_gap)
        print("diff response calls = " , exper.count_diffresponse_calls)
        print("response calls = " , exper.count_response_calls)
        output_file = output_file_no_extension + ".dat" 
        print(output_file)
        np.savetxt(output_file, max_gap)

    if MAKE_PLOT:
        output_file = output_file_no_extension + ".dat" 
        max_gap = np.loadtxt(output_file)
        print("max_gap = ", max_gap)    
        Plot_Upper_Limit(max_gap, plot_dots, plot_close = False, \
            plot_show = False)
        
