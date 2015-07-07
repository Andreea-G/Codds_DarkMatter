# -*- coding: utf-8 -*-
"""
Created on Sat Dec  6 15:32:34 2014

@author: Andreea
"""

from experiment_FoxMethod_jump import *
import matplotlib.pyplot as plt
linestyles = ['-', '--', '-.', ':']


def Plot_Upper_Limit(exper_name, upper_limit, HALO_DEP, kind=None, linewidth=3,
                     plot_dots=True, plot_close=True, plot_show=True):
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
    if not hasattr(Plot_Upper_Limit, "count"):
        Plot_Upper_Limit.count = {}
    if exper_name not in Plot_Upper_Limit.count:
        Plot_Upper_Limit.count[exper_name] = -1
    Plot_Upper_Limit.count[exper_name] += 1
    linestyle = linestyles[Plot_Upper_Limit.count[exper_name] % len(linestyles)]

    from scipy.interpolate import interp1d
    if plot_close:
        plt.close()

    # make a list of the x and y coordinates of the plots, and plot them
    if upper_limit.size == 0:   # nothing to plot
        print("upper_limit is empty!")
    elif upper_limit.ndim == 1:  # only one point, so no interpolation
        x = [upper_limit[0]]
        y = [upper_limit[1]]
        plt.plot(x, y, "o")
    else:   # more than one point, so decide on the interpolation order and plot
        x = upper_limit[:, 0]
        y = upper_limit[:, 1]
        num_points = x.size
        if num_points == 2 or kind == "linear":
            interp_kind = "linear"
        elif num_points == 3 or kind == "quadratic":
            interp_kind = "quadratic"
        else:
            interp_kind = "cubic"
        interp = interp1d(x, y, kind=interp_kind)
        x1 = np.linspace(x[0], x[-1], 1000)
        if plot_dots:
            plt.plot(x, y, "o")
        plt.plot(x1, interp(x1), linestyle=linestyle,
                 linewidth=linewidth, color=Color[exper_name.split()[0]])

    # set axis labels, depending on whether it is for halo-dependent or not
    if HALO_DEP is True:
        plt.xlabel('log$_{10}(m$ [GeV]$)$')
        plt.ylabel(r'log$_{10}(\sigma)$')
    else:
        plt.xlabel('$v_{min}$ [km/s]')
        plt.ylabel(r'$\eta$ $\rho$ $\sigma / m$ $[$days$^{-1}]$')

    # show plot
    if plot_show:
        plt.show()


def run_program(exper_name, scattering_type, mPhi, fp, fn, delta,
                RUN_PROGRAM, MAKE_REGIONS, MAKE_PLOT, HALO_DEP, FOX_METHOD,
                mx=None, mx_range=None, vmin_range=None, initial_energy_bin=None,
                vmin_FoxBand_range=None, logeta_FoxBand_percent_range=None,
                steepness=None, logeta_guess=None,
                vmin_index_list=None, logeta_index_range=None,
                confidence_levels=[0.9],
                filename_tail="", OUTPUT_MAIN_DIR="Output/", extra_tail="",
                plot_dots=True, quenching=None):
    ''' Main run of the program.
        Input:
            exper_name: name of experiment
            scattering_type: 'SI' for spin-dependent, 'SDPS' for pseudo-scalar, 'SDAV' for axial-vector
            mPhi: mass of mediator
            fp and fn: couplings to proton and neutron
            delta: DM mass split
            RUN_PROGRAM: True or False, whether the data should be (re-)computed
            MAKE_PLOT: True or False, whether the data should be plotted
            HALO_DEP: True or False, whether the analysis is halo-dependent or halo-independent
            FOX_METHOD: array of True or False, whether each step of the Fox Method is to be performed
            mx: DM mass, only give for halo-independent analysis
            mx_range: (mx_min, mx_max, num_steps) = DM mass range and number or steps, only for halo-dependent
            vmin_range: (vmin_min, vmin_max, vmin_step) = vmin range and step size, only for halo-independent
            initial_energy_bin: tuple or list of 2 elements, containing the starting energy bin. Used for DAMA
                combined analysis
            vmin_FoxBand_range: (vmin_Fox_min, vmin_Fox_max, vmin_Fox_numsteps) = vmin range and number of steps,
                used for calculating the Fox band.
            logeta_FoxBand_percent_range: (logeta_Fox_percent_min, logeta_Fox_percent_max, logeta_Fox_numsteps)
                = logeta percentage range and number of steps, used for calculating the Fox band.
                The min and max logeta are calculated as a given percentage above and below the optimum value.
            steepness: (steepness_vmin, steepness_vmin_center, steepness_logeta), parameters used for nonlinear
                sampling in vminStar and logetaStar. The higher the steepnesses the more points are taken close
                to the optimum values in vminStar and logetaStar.
            filename_tail: optional tag to be added to the file name
            OUTPUT_MAIN_DIR: name of main output directory, if different from "Output/"
            extra_tail: additional tail to be added to filenames for the Fox band.
            plot_dots: True or False, whether the plot should show the data points or just the interpolation
            quenching: quenching factor, needed for experiments that can have multiple options (such as KIMS or DAMA).
    '''

    print('name = ', exper_name)
    # Select which experiment class we must use, depending on what statistical analysis
    # we need.
    if HALO_DEP:
        print('Halo Dependent')
        if exper_name in MaximumGapLimit_exper:
            class_name = MaxGapExperiment
        elif exper_name in GaussianLimit_exper:
            class_name = GaussianExperiment
        elif exper_name in Poisson_exper:
            class_name = PoissonExperiment
        elif exper_name in DAMARegion_exper:
            class_name = DAMAExperiment
        elif exper_name.split()[0] in DAMARegion_exper:
            class_name = DAMAExperimentCombined
        elif exper_name in DAMALimit_exper:
            class_name = DAMATotalRateExperiment
        else:
            print("NotImplementedError: This experiment was not implemented!")
            return
    else:
        print('Halo Independent')
        if exper_name in FoxMethod_exper and np.any(FOX_METHOD):
            print('Fox Method')
            class_name = Experiment_FoxMethod
        elif exper_name in MaximumGapLimit_exper:
            class_name = MaxGapExperiment_HaloIndep
        elif exper_name in Poisson_exper:
            class_name = PoissonExperiment_HaloIndep
        elif exper_name in GaussianLimit_exper:
            class_name = GaussianExperiment_HaloIndep
        elif exper_name in DAMARegion_exper:
            class_name = Crosses_HaloIndep
        elif exper_name.split()[0] in DAMARegion_exper:
            class_name = Crosses_HaloIndep_Combined
        else:
            print("NotImplementedError: This experiment was not implemented!")
            return

        # if delta > 0 we have to use the integration in E recoil
        if delta > 0:
            class_name.__bases__ = (Experiment_HaloIndep_ER,)

    exper = class_name(exper_name, scattering_type, mPhi, quenching)

    # get the file name specific to the parameters used for this run
    output_file_no_extension = \
        Output_file_name(exper_name, scattering_type, mPhi, mx, fp, fn, delta,
                         HALO_DEP, filename_tail, OUTPUT_MAIN_DIR, quenching)

    # (re-)compute the data
    if RUN_PROGRAM:
        output_file = output_file_no_extension + "_temp.dat"
        f_handle = open(output_file, 'w')   # clear the file first
        f_handle.close()

        # calculate the upper limit, for halo-dependent or halo-independent
        if HALO_DEP:
            (mx_min, mx_max, num_steps) = mx_range
            upper_limit = exper.UpperLimit(fp, fn, delta, mx_min, mx_max, num_steps,
                                           output_file)
        else:
            (vmin_min, vmin_max, vmin_step) = vmin_range
            if not np.any(FOX_METHOD):
                upper_limit = \
                    exper.UpperLimit(mx, fp, fn, delta, vmin_min, vmin_max, vmin_step,
                                     output_file, initial_energy_bin=initial_energy_bin)
            else:
                if FOX_METHOD.ResponseTables:
                    exper.ResponseTables(vmin_min, vmin_max, vmin_step, mx, fp, fn, delta,
                                         output_file_no_extension)
                if FOX_METHOD.OptimalLikelihood:
                    if logeta_guess is None:
                        exper.OptimalLikelihood(output_file_no_extension)
                    else:
                        exper.OptimalLikelihood(output_file_no_extension,
                                                logeta_guess=logeta_guess)
                if FOX_METHOD.ImportOptimalLikelihood:
                    exper.ImportResponseTables(output_file_no_extension, plot=True)
                    exper.ImportOptimalLikelihood(output_file_no_extension, plot=True)
                    exper.PlotOptimum()
                if FOX_METHOD.ConstrainedOptimalLikelihood:
                    # Tests for delta = 0:
                    (vminStar, logetaStar) = (500, -25)
                    # Tests for delta = -50:
#                    (vminStar, logetaStar) = (185.572266287, -19.16840262)
                    exper.ImportOptimalLikelihood(output_file_no_extension)
                    exper.ConstrainedOptimalLikelihood(vminStar, logetaStar, plot=True)
                if np.any(FOX_METHOD[4:]):
                    if FOX_METHOD._fields[4] != 'VminLogetaSamplingTable':
                        raise AttributeError("FOX_METHOD's attribute is not as expected.")
                    (vmin_band_min, vmin_band_max, vmin_num_steps) = vmin_FoxBand_range
                    (logeta_percent_minus, logeta_percent_plus, logeta_num_steps) = \
                        logeta_FoxBand_percent_range
                    if steepness is not None:
                        (steepness_vmin, steepness_vmin_center, steepness_logeta) = \
                            steepness
                        print("Steepness: ", steepness_vmin, ", ",
                              steepness_vmin_center, ", ", steepness_logeta)
                        exper.VminSamplingList(output_file_no_extension,
                                               vmin_band_min, vmin_band_max, vmin_num_steps,
                                               steepness_vmin, steepness_vmin_center,
                                               plot=not np.any(FOX_METHOD[5:]))
                        exper.VminLogetaSamplingTable(output_file_no_extension,
                                                      logeta_percent_minus, logeta_percent_plus, logeta_num_steps,
                                                      steepness_logeta,
                                                      plot=not np.any(FOX_METHOD[5:]))
                    else:
                        print("Steepness: Default")
                        exper.VminSamplingList(output_file_no_extension,
                                               vmin_band_min, vmin_band_max, vmin_num_steps,
                                               plot=not np.any(FOX_METHOD[5:]))
                        exper.VminLogetaSamplingTable(output_file_no_extension,
                                                      logeta_percent_minus, logeta_percent_plus, logeta_num_steps,
                                                      plot=not np.any(FOX_METHOD[5:]))
                if FOX_METHOD.LogLikelihoodList:
                    print("vmin_FoxBand_range = ", vmin_band_min, " ",
                          vmin_band_max, " ", vmin_num_steps)
                    print("logeta_FoxBand_percent_range = ", logeta_percent_minus, " ",
                          logeta_percent_plus, " ", logeta_num_steps)
                    exper.LogLikelihoodList(output_file_no_extension, extra_tail=extra_tail,
                                            vmin_index_list=vmin_index_list,
                                            logeta_index_range=logeta_index_range)
                if FOX_METHOD.FoxBand:
                    exper.ImportOptimalLikelihood(output_file_no_extension)
                    interpolation_order = 2
                    delta_logL = [chi_squared1(c) for c in confidence_levels]
                    for d_logL in delta_logL:
                        multiplot = (d_logL == delta_logL[0]) and MAKE_PLOT
                        exper.FoxBand(output_file_no_extension, d_logL,
                                      interpolation_order, extra_tail=extra_tail,
                                      multiplot=multiplot)

        if HALO_DEP or not np.any(FOX_METHOD):
            print("upper_limit = ", upper_limit)
            print("diff response calls = ", exper.count_diffresponse_calls)
            print("response calls = ", exper.count_response_calls)
            output_file = output_file_no_extension + ".dat"
            print(output_file)  # write to file
            np.savetxt(output_file, upper_limit)

    # make regions
    if MAKE_REGIONS and exper_name.split()[0] in DAMARegion_exper:
        output_file = output_file_no_extension + ".dat"
        for CL in confidence_levels:
            output_file_regions = output_file_no_extension + \
                                  "_" + str(round(sigma_dev(CL), 2)) + "sigma"
            output_file_lower = output_file_regions + "_lower_limit.dat"
            output_file_upper = output_file_regions + "_upper_limit.dat"
            exper.Region(delta, CL, output_file, output_file_lower, output_file_upper)

    # produce plot
    if MAKE_PLOT and not np.any(FOX_METHOD[:-1]):
        if exper_name.split()[0] in DAMARegion_exper:
            for CL in confidence_levels:
                if hasattr(Plot_Upper_Limit, 'count'):
                    Plot_Upper_Limit.count[exper_name] = -1
                output_file_regions = output_file_no_extension + \
                    "_" + str(round(sigma_dev(CL), 2)) + "sigma"
                output_file_lower = output_file_regions + "_lower_limit.dat"
                output_file_upper = output_file_regions + "_upper_limit.dat"
                lower_limit = np.loadtxt(output_file_lower)
                upper_limit = np.loadtxt(output_file_upper)
                Plot_Upper_Limit(exper_name, lower_limit, HALO_DEP,
                                 plot_dots=plot_dots, plot_close=False, plot_show=False)
                Plot_Upper_Limit.count[exper_name] -= 1
                Plot_Upper_Limit(exper_name, upper_limit, HALO_DEP,
                                 plot_dots=plot_dots, plot_close=False, plot_show=False)
        else:
            output_file = output_file_no_extension + ".dat"
            upper_limit = np.loadtxt(output_file)
            print("upper_limit = ", upper_limit)
            Plot_Upper_Limit(exper_name, upper_limit, HALO_DEP,
                             plot_dots=plot_dots, plot_close=False, plot_show=False)

    # make band plot
    if FOX_METHOD.FoxBandPlot and exper_name == "CDMSSi2012":
        output_file = output_file_no_extension + ".dat"
        exper.ImportOptimalLikelihood(output_file_no_extension)
        interp_kind = 'cubic'

#        exper.PlotSamplingTable(output_file_no_extension,
#                                plot_close=False, plot_show=False, plot_optimum=False)
        delta_logL = [chi_squared1(c) for c in confidence_levels]
        print("delta_logL =", delta_logL)
        for d_logL in delta_logL:
            exper.ImportFoxBand(output_file_no_extension, d_logL, extra_tail=extra_tail)
            Plot_Upper_Limit(exper_name, exper.vmin_logeta_band_low, HALO_DEP,
                             kind=interp_kind,
                             plot_dots=plot_dots, plot_close=False, plot_show=False)
            Plot_Upper_Limit(exper_name, exper.vmin_logeta_band_up, HALO_DEP,
                             kind=interp_kind,
                             plot_dots=plot_dots, plot_close=False, plot_show=False)
        exper.PlotOptimum(ylim_percentage=(1.2, 0.8), plot_close=False, plot_show=False)
