# -*- coding: utf-8 -*-
"""
Created on Thu Nov 20 22:52:11 2014

@author: Andreea
"""

from __future__ import print_function
from __future__ import division
from experiment import *
import profile
     
def Plot_Upper_Limit(max_gap):
    import matplotlib.pyplot as plt
    from scipy.interpolate import interp1d
    
    x = max_gap[:,0]
    y = -np.log10(3.*SpeedOfLight**2*1e4*3600*24) + max_gap[:,0] + max_gap[:,1]
    interp = interp1d(x, y, kind = "cubic")
    plt.close
    x1 = np.linspace(x[0],x[-1],50)
    plt.plot(x, y, "o", x1, interp(x1))
    plt.xlabel('Log10(m [GeV])')
    plt.ylabel('Log10(sigma)')
    plt.show()

     
def main():
#    exper_name = "CDMSlite2013CoGeNTQ"
    exper_name = "superCDMS"
#    exper_name = "LUX2013zero"
    scattering_type = 'SD66'
    mPhi = 1000.
    fp = 1.
    fn = 0.
    delta = 0.
    
    exper = Experiment(exper_name, scattering_type, mPhi)
    print('name = ', exper.name)
    
    output_dir = OutputDirectory(OUTPUT_MAIN_DIR, scattering_type, mPhi, delta)
    output_file_no_extension = "./" + output_dir + "UpperLimitSHM_" + exper.name + "_mxsigma" \
        + FileNameTail(fp, fn)# + "_Dirac"
    print(output_file_no_extension)
    

    RUN_PROGRAM = False
    MAKE_PLOT = True

    if RUN_PROGRAM:          
        mx_min = 3.5
        mx_max = 100.
        num_steps = 30
        output_file = output_file_no_extension + "_py_temp.dat" 
        f_handle = open(output_file, 'w')   # clear the file first
        f_handle.close()
        
        max_gap = exper.MaximumGapLimit(fp, fn, delta, mx_min, mx_max, num_steps, output_file)
        print("max gap = ", max_gap)
        output_file = output_file_no_extension + "_py.dat" 
        print(output_file)
        np.savetxt(output_file, max_gap)

    if MAKE_PLOT:
        output_file = output_file_no_extension + "_py.dat" 
        max_gap = np.loadtxt(output_file)
        print("max_gap = ", max_gap)    
        Plot_Upper_Limit(max_gap)

    
    TEST_INT_RESPONSE = False
    if TEST_INT_RESPONSE:
        mx = 10.
        Eee1 = 2
        Eee2 = 
    #    ER = 5.
    #    diff_resp = exper.DifferentialResponseSHM(ER, Eee1, mx, fp, fn, delta)
    #    print("diff response = ", diff_resp)
    #    print("response = ", exper.ResponseSHM(ER, Eee1, Eee2, mx, fp, fn, delta))
        print("int response = ", exper.IntegratedResponseSHM(Eee1, Eee2, mx, fp, fn, delta))
        print("diff response calls = " , exper.count_diffresponse_calls)
        print("response calls = " , exper.count_response_calls)

    
if __name__ == '__main__':
#    main()
    profile.run("main()")

import num