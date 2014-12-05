# -*- coding: utf-8 -*-
"""
Created on Thu Nov 20 22:52:11 2014

@author: Andreea
"""

from __future__ import print_function
from __future__ import division
from experiment import *
import profile
import matplotlib.pyplot as plt
     
def Plot_Upper_Limit(max_gap, plot_close = True, plot_show = True, plot_dots = True):
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
    RUN_PROGRAM, MAKE_PLOT):
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
    print(output_file_no_extension)

    if RUN_PROGRAM:  
        output_file += "_temp.dat" 
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
        Plot_Upper_Limit(max_gap, plot_close = False, \
            plot_show = False, plot_dots = False)
        
    
    
def main():
    implemented_exper = ["superCDMS", \
        "LUX2013zero", "LUX2013one", "LUX2013three", "LUX2013five", "LUX2013many", \
        "XENON10", "CDMSlite2013CoGeNTQ"]
    scattering_type = 'SD66'
    mPhi = 1000.
    fp = 1.
    fn = 0.
    delta = 0.

    mx_min = 10
    mx_max = 100.
    num_steps = 1

    global v0bar, vobs, vesc
    v0bar = 230 - 3 * 24.4
    vobs = v0bar + 12
    vesc = 544 - 3 * 39

    RUN_PROGRAM = False
    MAKE_PLOT = False

    exper_list = implemented_exper[0:1]
    plt.close()
    for exper_name in exper_list:
        run_program(exper_name, scattering_type, mPhi, fp, fn, delta, mx_min, mx_max, num_steps, \
            RUN_PROGRAM, MAKE_PLOT)
    plt.show()
    
if __name__ == '__main__':
    main()
#    profile.run("main()")

