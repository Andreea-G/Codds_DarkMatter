# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 00:47:37 2015

@author: Andreea
"""

#TODO! This only works for CDMSSi! 

from experiment_HaloIndep import *
from interp import interp1d
#from scipy.interpolate import interp1d
from scipy.optimize import minimize, brentq, fsolve
import matplotlib.pyplot as plt
import os   # for speaking

DEBUG = T

class Experiment_FoxMethod(Experiment_HaloIndep):
    def __init__(self, expername, scattering_type, mPhi = mPhiRef):
        Experiment_HaloIndep.__init__(self, expername, scattering_type, mPhi)
        module = import_file(INPUT_DIR + expername + ".py")
        self.ERecoilList = module.ERecoilList
        self.mu_BKG_i = module.mu_BKG_i
        self.NBKG = module.NBKG

    def VMinSortedList(self, mx, fp, fn, delta):
        ''' The list of vmin corresponsing to measured recoil energies, sorted in increasing order. 
            Will be useful as starting guesses.
        '''
        self.vmin_sorted_list = np.sort(VMin(self.ERecoilList, self.mT[0], mx, delta))
        return

    def ResponseTables(self, vmin_min, vmin_max, vmin_step, mx, fp, fn, delta, output_file_tail):
        ''' Computes response tables:
            - self.diff_response_tab is a table of [vmin, DifferentialResponse(Eee_i)] pairs for each vmin in the range 
        [vminmin, vminmax], corresponding to measured recoil energies Eee_i. 
        It is a 3D matrix where
                axis = 0 has dimension self.ERecoilList.size()
                axis = 1 has dimension vmin_list.size() + 1 (where +1 is because we prepend zeros for vmin = 0)
                axis = 2 has dimension 2 for the pairs of [vmin, diff_response]
            - self.response_tab is a table of [vmin, Response] pairs for each vmin in the range 
        [vminmin, vminmax], corresponding to DifferentialResponse integrated over the full energy range.
        It is a 2D matrix where 
                axis = 1 has dimension vmin_list.size() + 1 (where +1 is because we prepend zeros for vmin = 0)
                axis = 2 has dimension 2 for the pairs of [vmin, diff_response]
            Input:
                vmin range and vmin step size
                mx, fp, fn, delta
                output_file_tail since the results for self.vmin_sorted_list, self.diff_response_tab and 
        self.response_tab are each written to files.
        '''
        self.VMinSortedList(mx, fp, fn, delta)
        # write it to file
        file =  output_file_tail + "_VMinSortedList.dat"
        print(file)  # write to file
        np.savetxt(file, self.vmin_sorted_list)

        if delta == 0:
            branches = [1]
        else:
            branches = [1, -1]
        self.vmin_linspace = np.linspace(vmin_min, vmin_max, (vmin_max - vmin_min)/vmin_step + 1)
        self.diff_response_tab = np.zeros((self.ERecoilList.size, 1))  
        self.response_tab = np.zeros(1)
        for vmin in self.vmin_linspace:
            diff_resp_list = np.zeros((1,3))
            resp = 0
            for sign in branches:
                (qER, const_factor) = self.ConstFactor(vmin, mx, fp, fn, delta, sign)
                diff_resp_list += np.array([self.DifferentialResponse(Eee, qER, const_factor) for Eee in self.ERecoilList])
                resp += integrate.quad(self.DifferentialResponse, self.Ethreshold, self.Emaximum, \
                    args=(qER, const_factor), epsrel = PRECISSION, epsabs = 0)[0]
            self.diff_response_tab = np.append(self.diff_response_tab, diff_resp_list.transpose(), axis = 1)
            self.response_tab = np.append(self.response_tab, [resp], axis = 0)
        self.vmin_linspace = np.insert(self.vmin_linspace, 0., 0)
        # write to file
        file = output_file_tail + "_VminLinspace.dat"
        print(file)  # write to file
        np.savetxt(file, self.vmin_linspace)
        file = output_file_tail + "_DiffRespTable.dat"
        print(file)  # write to file
        np.savetxt(file, self.diff_response_tab)
        file = output_file_tail + "_RespTable.dat"
        print(file)  # write to file
        np.savetxt(file, self.response_tab)
        return

    def PlotDifferentialResponse(self):
        plt.close()
        for i in range(self.ERecoilList.size):
            plt.plot(self.vmin_linspace, np.array([self.diff_response_interp[i](v) for v in self.vmin_linspace]))
        plt.show()

    def PlotResponse(self):
        plt.close()
        plt.plot(self.vmin_linspace, np.array([self.response_interp(v) for v in self.vmin_linspace]))
        plt.show()
        
    def ImportResponseTables(self, output_file_tail, plot = True):
        file =  output_file_tail + "_VMinSortedList.dat"
        with open(file, 'r') as f_handle:
            self.vmin_sorted_list = np.loadtxt(f_handle)
        file = output_file_tail + "_VminLinspace.dat"
        with open(file, 'r') as f_handle:
            self.vmin_linspace = np.loadtxt(f_handle)
        file = output_file_tail + "_DiffRespTable.dat"
        with open(file, 'r') as f_handle:
            self.diff_response_tab = np.loadtxt(f_handle)
        file = output_file_tail + "_RespTable.dat"
        with open(file, 'r') as f_handle:
            self.response_tab = np.loadtxt(f_handle)
        self.diff_response_interp = np.array([interp1d(self.vmin_linspace, dr) for dr in self.diff_response_tab])
        self.response_interp = interp1d(self.vmin_linspace, self.response_tab)
        if plot:
            self.PlotDifferentialResponse()
            self.PlotResponse()
        return

    def VminIntegratedResponseTable(self, vmin_list):
        return np.array([[integrate.quad(self.diff_response_interp[i], vmin_list[a], vmin_list[a + 1], epsrel = PRECISSION, epsabs = 0)[0] \
            for a in range(vmin_list.size - 1)] for i in range(self.ERecoilList.size)])

    def IntegratedResponseTable(self, vmin_list):
        return np.array([integrate.quad(self.response_interp, vmin_list[a], vmin_list[a + 1], epsrel = PRECISSION, epsabs = 0)[0] \
            for a in range(vmin_list.size - 1)])

    def MinusLogLikelihood(self, vars_list, vminStar = None, logetaStar = None, vminStar_index = None):
#        print("vars_list = ", vars_list)
#        print("Star vars: ", vminStar, " ", logetaStar)
        if vminStar == None:
            vmin_list_w0 = vars_list[: vars_list.size/2]
            logeta_list = vars_list[vars_list.size/2 :]
        else:
            vmin_list_w0 = np.insert(vars_list[: vars_list.size/2], vminStar_index, vminStar)
            logeta_list = np.insert(vars_list[vars_list.size/2 :], vminStar_index, logetaStar)
        vmin_list_w0 = np.insert(vmin_list_w0, 0, 0)
        vmin_resp_integr = self.VminIntegratedResponseTable(vmin_list_w0)
        resp_integr = self.IntegratedResponseTable(vmin_list_w0)
        if np.any(vmin_resp_integr < 0):
            print("@@@@@@@@@@@@@ vmin_resp_integr < 0 @@@@@@@@@@@@@")
            print("vmin_list_w0 = ", vmin_list_w0)
            print("vmin_resp_integr = ", vmin_resp_integr)
            print("@@@@@@@@@@@@@")
            print("vmin_linspace = ", self.vmin_linspace)
            print("diff resp = ", np.array([self.diff_response_interp[0](v) for v in self.vmin_linspace]))
            self.PlotDifferentialResponse()
            raise ValueError
        if np.any(resp_integr < 0):
            print("@@@@@@@@@@@@@ resp_integr < 0 @@@@@@@@@@@@@")
            print("vmin_list_w0 = ", vmin_list_w0)
            print("resp_integr = ", resp_integr)
            print("@@@@@@@@@@@@@")
            self.PlotResponse()
            raise ValueError
        mu_i = self.Exposure * np.dot(vmin_resp_integr, 10**logeta_list)
        Nsignal = self.Exposure * np.dot(10**logeta_list, resp_integr)
        result = self.NBKG + Nsignal - np.log(self.mu_BKG_i + mu_i).sum()
#        print("result = ", result)
        return result

    def _MinusLogLikelihood(self, vars_list, constr_func = None, vminStar = None, logetaStar = None, vminStar_index = None):
        if DEBUG:
            print("*** vars_list = ", vars_list)
        constraints = constr_func(vars_list)
        if DEBUG:
            print("constraints = ", constraints)
        constr_not_valid = constraints < 0
        if DEBUG:
            print("constr_not_valid = ", constr_not_valid)
        if np.any(constr_not_valid):
            constr_list = constraints[constr_not_valid]
            if DEBUG: 
                print("Constraints not valid!!********************************************************************************")
                print("constr sum = ", -constr_list.sum())
            return -constr_list.sum() * 10**6
        return self.MinusLogLikelihood(vars_list, vminStar = vminStar, logetaStar = logetaStar, vminStar_index = vminStar_index)

    def OptimalLikelihood(self, output_file_tail, logeta_guess = -23.85):
        self.ImportResponseTables(output_file_tail, plot = False)
        vars_guess = np.append(self.vmin_sorted_list, logeta_guess * np.ones(self.vmin_sorted_list.size))
        print("vars_guess = ", vars_guess)
        vmin_max = self.vmin_linspace[-1]
        def constr_func(x, vmin_max = vmin_max): 
            #  0 -  8: bounds: 3 * (x.size/2) constraints = 9 for x.size/2 = 3
            #  9 - 12: sorted array: 2 * (x.size/2 - 1) constraints = 4 for x.size/2 = 3
            constraints = np.concatenate([x[:x.size/2], vmin_max - x[:x.size/2], -x[x.size/2:], \
                np.diff(x[:x.size/2]), np.diff(-x[x.size/2:])])
            is_not_close = np.logical_not(np.isclose(constraints, np.zeros_like(constraints), atol = 1e-5))
            is_not_close[:3 * (x.size/2)] = T
            res = np.where(is_not_close, constraints, np.abs(constraints))
            if DEBUG:
                print("***constr = ", res)
                print("tf = ", res > 0)
            return res
        constr = ({'type': 'ineq', 'fun': constr_func})
#        optimum_log_likelihood = minimize(self.MinusLogLikelihood, vars_guess, constraints = constr)
        optimum_log_likelihood = minimize(self._MinusLogLikelihood, vars_guess, args = (constr_func,), constraints = constr)
        print(optimum_log_likelihood)
        print("vars_guess = ", repr(vars_guess))
        file = output_file_tail + "_GloballyOptimalLikelihood.dat"
        print(file)  # write to file
        np.savetxt(file, np.append([optimum_log_likelihood.fun], optimum_log_likelihood.x))
        return 

    def ImportOptimalLikelihood(self, output_file_tail):
        self.ImportResponseTables(output_file_tail, plot = False)
        file = output_file_tail + "_GloballyOptimalLikelihood.dat"
        with open(file, 'r') as f_handle:
            optimal_result = np.loadtxt(f_handle)
        self.optimal_logL = optimal_result[0]
        self.optimal_vmin = optimal_result[1 : optimal_result.size/2 + 1]
        self.optimal_logeta = optimal_result[optimal_result.size/2 + 1 :]
        return

    def PlotStepFunction(self, vmin_list, logeta_list, xlim_percentage = (0., 1.1), ylim_percentage = (1.01, 0.99), plot_close = True, plot_show = True):
        if plot_close: 
            plt.close()
        print(vmin_list)
        print(logeta_list)
        x = np.insert(vmin_list, 0, 0)
        y = np.insert(logeta_list, 0, logeta_list[0])
        plt.step(x, y)
        plt.plot(x, y, 'o')
        plt.xlim([vmin_list[0] * xlim_percentage[0], vmin_list[-1] * xlim_percentage[1]])
        plt.ylim([max(logeta_list[-1] * ylim_percentage[0], -60), max(logeta_list[0] * ylim_percentage[1], -35)])
        if plot_show:
            plt.show()
        return

    def PlotOptimum(self, xlim_percentage = (0., 1.1), ylim_percentage = (1.01, 0.99), plot_close = True, plot_show = True):
        self.PlotStepFunction(self.optimal_vmin, self.optimal_logeta,\
            xlim_percentage = xlim_percentage, ylim_percentage = ylim_percentage,
            plot_close = plot_close, plot_show = plot_show)
        return

    def PlotConstrainedOptimum(self, vminStar, logetaStar, vminStar_index, \
        xlim_percentage = (0., 1.1), ylim_percentage = (1.01, 0.99), plot_close = True, plot_show = True):
        self.PlotStepFunction(self.optimal_vmin, self.optimal_logeta, plot_close = plot_close, plot_show = False)
        self.PlotStepFunction(np.insert(self.constr_optimal_vmin, vminStar_index, vminStar), \
            np.insert(self.constr_optimal_logeta, vminStar_index, logetaStar), \
            xlim_percentage = xlim_percentage, ylim_percentage = ylim_percentage,
            plot_close = False, plot_show = False)
        plt.plot(vminStar, logetaStar, '*')
        if plot_show:
            plt.show()
        return
        
    def ConstrainedOptimalLikelihood(self, vminStar, logetaStar, logeta_guess = -25., plot = False):
        vminStar_index = 0
        while vminStar_index < self.optimal_vmin.size and vminStar > self.optimal_vmin[vminStar_index]:
            vminStar_index += 1
        vmin_max = self.vmin_linspace[-1]
        if DEBUG:
            print("vminStar_index =", vminStar_index)
        logeta_guess = self.optimal_logeta
        logeta_guess = np.concatenate([np.maximum(logeta_guess[:vminStar_index], np.ones(vminStar_index)*logetaStar), \
            np.minimum(logeta_guess[vminStar_index:], np.ones(logeta_guess.size - vminStar_index)*logetaStar)])
        vars_guess = np.append(self.optimal_vmin, logeta_guess)
                    
        def constr_func(x, vminStar = vminStar, logetaStar = logetaStar, vmin_max = vmin_max, vminStar_index = vminStar_index): 
            #  0 -  8: bounds: 3 * (x.size/2) constraints = 9 for x.size/2 = 3
            #  9 - 12: sorted array: 2 * (x.size/2 - 1) constraints = 4 for x.size/2 = 3
            # 13 - 15: vminStar_index: x.size/2 constraints = 3 for x.size/2 = 3
            # 16 - 18: vminStar and logetaStar: x.size/2 constraints = 3 for x.size/2 = 3
            constraints = np.concatenate([x[:x.size/2], vmin_max - x[:x.size/2], -x[x.size/2:], \
                np.diff(x[:x.size/2]), np.diff(-x[x.size/2:]), \
                vminStar - x[:vminStar_index], x[vminStar_index : x.size/2] - vminStar, \
                (x[:x.size/2] - vminStar) * (-x[x.size/2:] + logetaStar)])
            is_not_close = np.logical_not(np.isclose(constraints, np.zeros_like(constraints), atol = 1e-5))
            is_not_close[:3 * (x.size/2)] = T
            res = np.where(is_not_close, constraints, np.abs(constraints))
            if DEBUG:
                print("***constr = ", res)
                print("tf = ", res > 0)
            return res
        
        constr = ({'type': 'ineq', 'fun': constr_func})

        constr_optimum_log_likelihood = minimize(self._MinusLogLikelihood, vars_guess, \
            args = (constr_func, vminStar, logetaStar, vminStar_index), constraints = constr)
        if DEBUG:
            print(constr_optimum_log_likelihood)
            constraints = constr_func(constr_optimum_log_likelihood.x)
            print("constraints = ", repr(constraints))
            is_not_close = np.logical_not(np.isclose(constraints, np.zeros_like(constraints)))
            constr_valid = np.logical_not(np.logical_and(constraints < 0, is_not_close))
            print("constr_valid = ", repr(constr_valid))
            print("MinusLogLikelihood = ", self.MinusLogLikelihood(constr_optimum_log_likelihood.x, \
                vminStar = vminStar, logetaStar = logetaStar, vminStar_index = vminStar_index))
            print("vars_guess = ", repr(vars_guess))
            print("optimum_logL = ", self.optimal_logL)
            print("vminStar = ", vminStar)
            print("logetaStar = ", logetaStar)

        self.constr_optimal_logl = constr_optimum_log_likelihood.fun
        vars_result = constr_optimum_log_likelihood.x

        self.constr_optimal_vmin = vars_result[: vars_result.size/2]
        self.constr_optimal_logeta = vars_result[vars_result.size/2:]
#        print(self.constr_optimal_vmin)
#        print(self.constr_optimal_logeta)
        if plot:
            os.system("say 'Finished plot'")
            self.PlotConstrainedOptimum(vminStar, logetaStar, vminStar_index, xlim_percentage = (0., 1.1), ylim_percentage = (1.2, 0.8))
        return self.constr_optimal_logl

            
    def VminSamplingList(self, output_file_tail, vmin_min, vmin_max, vmin_num_steps, \
        steepness_vmin = 1.5, steepness_vmin_center = 2.5, plot = False):
        self.ImportOptimalLikelihood(output_file_tail)
        xmin = vmin_min
        xmax = vmin_max
        # TODO! This +4 is to compensate for a loss of ~4 points (not always 4 though), 
        # and it's due to taking floor later on.
        # Fid a better way to deal with this.
        x_num_steps = vmin_num_steps #+ 4  
        s = steepness_vmin
        sc = steepness_vmin_center
        
        x_lin = np.linspace(xmin, xmax, 1000)
        x0_list = self.optimal_vmin
        numx0 = x0_list.size
        
        print("x0 = ", x0_list)
        UnitStep = lambda x: (np.sign(x) + 1) / 2
        g1 = lambda x, x0, s0, xmin = xmin: np.log10(UnitStep(x - x0) + \
            UnitStep(x0 - x) * (x0 - xmin)/ (x + 10**s0 * (-x + x0) - xmin))
        g2 = lambda x, x0, s0, xmax = xmax: np.log10(UnitStep(x0 - x) + \
            UnitStep(x - x0) * (x + 10**s0 * (-x + x0) - xmax)/ (x0 - xmax))        
#        print("g2 = ", g2(x, x0_list, sc))
        g = lambda x, x0, s1, s2: g1(x, x0, s1) + g2(x, x0, s2)
#        print("g = ", g(x, x0_list, sc, sc))
        
        s_list = np.array([[s, sc]] + [[sc,sc]] * (numx0 - 2) + [[sc, s]])
        g_total = lambda x, sign = 1, x0 = x0_list, s_list = s_list: \
            np.array([sign * g(x, x0_list[i], s_list[i,0], s_list[i,1]) \
            for i in range(x0_list.size)]).prod(axis = 0)
        g_lin = g_total(x_lin)
        
        xT_guess = (x0_list[:-1] + x0_list[1:]) / 2
#        print("xT_guess = ", xT_guess)
        bounds = np.array([(x0_list[i], x0_list[i + 1]) for i in range(x0_list.size - 1)])
#        print("bounds = ", bounds)
#        print("gtot = ", g_total(xT_guess))
        x_turns_max = np.array([minimize(g_total, np.array(xT_guess[i]), \
            args = (-1,), bounds = [bounds[i]]).x for i in range(0, xT_guess.size, 2)])
        x_turns_min = np.array([minimize(g_total, np.array(xT_guess[i]), \
            bounds = [bounds[i]]).x for i in range(1, xT_guess.size, 2)])
        x_turns = np.sort(np.append(x_turns_max, x_turns_min))
        x_turns = np.append(np.insert(x_turns, 0, xmin), [xmax])
        y_turns = g_total(x_turns)
        print("x_turns = ", x_turns)
        print("y_turns = ", y_turns)
        
        g_inverse = lambda y, x1, x2: brentq(lambda x: g_total(x) - y, x1, x2)
        g_inverse_list = lambda y_list, x1, x2: np.array([g_inverse(y, x1, x2) for  y in y_list])
        y_diff = np.diff(y_turns)
        y_diff_sum = np.abs(y_diff).sum()
        print("y_diff = ", y_diff)
        num_steps = np.array([max(1, np.floor(x_num_steps * np.abs(yd)/y_diff_sum)) for yd in y_diff])
        print("num_steps = ", num_steps)
        y_list = np.array([np.linspace(y_turns[i], y_turns[i+1], num_steps[i]) for i in range(num_steps.size)])
        x_list = np.array([g_inverse_list(y_list[i], x_turns[i], x_turns[i+1]) for i in range(y_list.size)])
        x_list = np.concatenate(x_list)
        y_list = np.concatenate(y_list)
        x_list = x_list[np.array([x_list[i] != x_list[i+1] for i in range(x_list.size - 1)] + [True])]
        y_list = y_list[np.array([y_list[i] != y_list[i+1] for i in range(y_list.size - 1)] + [True])]
        self.vmin_sampling_list = x_list
        
        if plot:
            plt.close()
            plt.plot(x_lin, g_lin)
            plt.plot(x_turns, y_turns, 'o')
            plt.plot(x_list, y_list, '*')
            plt.xlim([xmin, xmax])
            plt.ylim([-s * sc**(numx0 - 1), s * sc**(numx0 - 1)])
            plt.show()
        
        return

    def OptimumStepFunction(self, vmin):
        index = 0
        while index < self.optimal_vmin.size and vmin > self.optimal_vmin[index]:
            index += 1
        if index == self.optimal_vmin.size:
            return self.optimal_logeta[-1]*10
        return self.optimal_logeta[index]

    def VminLogetaSamplingTable(self, output_file_tail, \
        logeta_percent_minus, logeta_percent_plus, logeta_num_steps, \
        steepness_logeta = 1, plot = False):
        print(self.optimal_vmin)
        print(self.optimal_logeta)
        
        logeta_num_steps_minus = logeta_num_steps * \
            logeta_percent_minus / (logeta_percent_minus + logeta_percent_plus)
        logeta_num_steps_plus = logeta_num_steps * \
            logeta_percent_plus / (logeta_percent_minus + logeta_percent_plus)
        
        s = steepness_logeta
        
        f = lambda x, xm, i, s0 = s: \
            (xm - x) / (10**s0 - 1) * 10**i + (10**s0 * x - xm) / (10**s0 - 1)
        self.vmin_logeta_sampling_table = []
        vmin_last_step = self.optimal_vmin[-1]
        for vmin in self.vmin_sampling_list:
            logeta_opt = self.OptimumStepFunction(min(vmin, vmin_last_step))
            logeta_min = logeta_opt * (1 + logeta_percent_minus)
            logeta_max = logeta_opt * (1 - logeta_percent_plus)
            logeta_list_minus = [[vmin, f(logeta_opt, logeta_min, i)] for i in np.linspace(s, 0, logeta_num_steps_minus)]
            logeta_list_plus = [[vmin, f(logeta_opt, logeta_max, i)] for i in np.linspace(s / logeta_num_steps_plus, s, logeta_num_steps_plus)]
            self.vmin_logeta_sampling_table += [logeta_list_minus + logeta_list_plus]
        
        self.vmin_logeta_sampling_table = np.array(self.vmin_logeta_sampling_table)
        
        if plot:
            plt.close()
            for tab in self.vmin_logeta_sampling_table:
                plt.plot(tab[:,0], tab[:,1], 'o')
            self.PlotOptimum(xlim_percentage = (0.9, 1.1), ylim_percentage = (1.2, 0.8), plot_close = False, plot_show = True)
        return
    
    def ImportSamplingTable(self, output_file_tail):
        self.ImportOptimalLikelihood(output_file_tail)
        file = output_file_tail + "_vminLogetaSamplingTable.dat"
        with open(file, 'r') as f_handle:
            self.vmin_logeta_sampling_table = np.loadtxt(f_handle)
        return
    
    def LogLikelihoodList(self, output_file_tail):
        ''' Gives a list of the form [[logetaStar_i0, logL_i0], [logetaStar_i1, logL_i1], ...] needed for 1D interpolation, 
        where i is the index corresponding to vminStar_i.
            Imports:
            - vmin_logeta_sampling_table list of {vminStar, logetaStar} at which the constrained optimal likelihood will be computed. 
            Produces a table of the form:
            - vmin_logeta_sampling_table[i] is a list of [[vminStar_i, logetaStar_0], [vminStar_i, logetaStar_1], ...] corresponding to vminStar_i
        '''
        self.logL_list_table = []
        temp_file = output_file_tail + "_LogetaStarLogLikelihoodList.dat"
        f_handle = open(temp_file, 'w')   # clear the file first
        f_handle.close()
        
        for index in range(0, self.vmin_logeta_sampling_table.shape[0]):
            print("index = ", index)
            vminStar = self.vmin_logeta_sampling_table[index, 0, 0]
            print("vminStar = ", vminStar)
            logetaStar_list = self.vmin_logeta_sampling_table[index, :, 1]
#            print("logetaStar_list = ", logetaStar_list)
            try:
                table = np.array([[logetaStar, self.ConstrainedOptimalLikelihood(vminStar, logetaStar)] \
                    for logetaStar in logetaStar_list])
            except:
                os.system("say 'Finished program'")
            print("table = ", table)
            self.logL_list_table += [table]
            with open(temp_file,'ab') as f_handle:
                np.savetxt(f_handle, table)
#            print("self.logL_list_table = ", self.logL_list_table)
        
        self.logL_list_table = np.concatenate(self.logL_list_table)
        
        file = output_file_tail + "_LogetaStarLogLikelihoodList.dat"
        print(file)  # write to file
        np.savetxt(file, self.logL_list_table)
    
    def _logL_interp(vars_list, constraints):
        constr_not_valid = constraints(vars_list)[:-1] < 0
        if np.any(constr_not_valid):
            constr_list = constraints(vars_list)[constr_not_valid]
            return -constr_list.sum() * 10**2
        return logL_interp(vars_list)


    def FoxBand(self, output_file_tail, delta_logL, interpolation_order, multiplot = True, plot = True):
        file = output_file_tail + "_LogetaStarLogLikelihoodList.dat"
        with open(file, 'r') as f_handle:
            self.logL_list_table = np.loadtxt(f_handle)
        
        n = self.vmin_sampling_list.size
        shape = self.logL_list_table.shape[0]
        self.logL_list_table = self.logL_list_table.reshape(n, shape/n, self.logL_list_table.shape[1])
#        print("self.logL_list_table = ", self.logL_list_table)
        print("self.vmin_sampling_list = ", self.vmin_sampling_list)
        print("self.logL_list_table = ", self.logL_list_table)
        
        self.vmin_logeta_band_low = []
        self.vmin_logeta_band_up = []
        vmin_last_step = self.optimal_vmin[-1]
        if multiplot:
            plt.close()
        for index in range(self.logL_list_table.shape[0]):
            print("vmin = ", self.vmin_sampling_list[index])
            logeta_optim = self.OptimumStepFunction(min(self.vmin_sampling_list[index], vmin_last_step))
            x = self.logL_list_table[index, :, 0]   # this is logeta
            y = self.logL_list_table[index, :, 1]   # this is logL
            logL_interp = interp1d(x, y)
            
            def _logL_interp(vars_list, constraints):
                constr_not_valid = constraints(vars_list)[:-1] < 0
                if np.any(constr_not_valid):
                    constr_list = constraints(vars_list)[constr_not_valid]
                    return -constr_list.sum() * 10**2
                return logL_interp(vars_list)
                        
            print(self.optimal_logL - delta_logL)
            print(np.array([self.logL_list_table[index, 0, 0]]), " ", self.logL_list_table[index, -1, 0])
            print(logeta_optim)
            constr_func = lambda logeta, logeta_min = np.array([self.logL_list_table[index, 0, 0]]), \
                logeta_max = np.array([self.logL_list_table[index, -1, 0]]): \
                np.concatenate([logeta - logeta_min, logeta_max - logeta])
            constr = ({'type': 'ineq', 'fun': constr_func})
            logeta_minimLogL = minimize(_logL_interp, np.array([logeta_optim]), args = (constr_func,), constraints = constr).x
            print(logeta_minimLogL)
            
            print("x = ", x)
            print("y = ", y)
            if multiplot:
                plt.close()
                plt.plot(x, y, 'o-')
                plt.plot(x, (self.optimal_logL + delta_logL) * np.ones_like(y))
                plt.xlim(x[0], x[-1])
                plt.ylim(-5, 20)
                plt.show()

            try:
                if y[0] > self.optimal_logL + delta_logL and logeta_minimLogL < self.optimal_logL + delta_logL:
                    self.vmin_logeta_band_low += [[self.vmin_sampling_list[index], \
                        brentq(lambda logeta: logL_interp(logeta) - self.optimal_logL - delta_logL, \
                            self.logL_list_table[index, 0, 0], logeta_minimLogL)]]
                if y[-1] > self.optimal_logL + delta_logL and logeta_minimLogL < self.optimal_logL + delta_logL:
                    print("a, b: ", logeta_optim, " ", self.logL_list_table[index, -1, 0])
                    print("f: ", logL_interp(logeta_optim) - self.optimal_logL - delta_logL, " ", \
                        logL_interp(self.logL_list_table[index, -1, 0]) - self.optimal_logL - delta_logL)
                    self.vmin_logeta_band_up += [[self.vmin_sampling_list[index], \
                        brentq(lambda logeta: logL_interp(logeta) - self.optimal_logL - delta_logL, \
                            logeta_minimLogL, self.logL_list_table[index, -1, 0])]]
#            except ValueError:
#                plt.close()
#                plt.plot(x, y, 'o-')
#                plt.plot(x, (self.optimal_logL + delta_logL) * np.ones_like(y))
#                plt.xlim(x[0], x[-1])
#                plt.ylim([3,10])
#
#                plt.plot(x, y, 'o')
#                plt.plot(logeta_optim, logL_interp(logeta_optim), '*')
#                plt.plot(logeta_optim, self.optimal_logL, '*')
#                print("ValueError")
#                plt.show()
##                pass
            finally:
                None
        if multiplot:
            plt.show()
        
        self.vmin_logeta_band_low = np.array(self.vmin_logeta_band_low)
        self.vmin_logeta_band_up = np.array(self.vmin_logeta_band_up)
        
        print("lower band: ", self.vmin_logeta_band_low)
        print("upper band: ", self.vmin_logeta_band_up)        
        
        plt.close()
        plt.plot(self.vmin_logeta_band_low[:,0], self.vmin_logeta_band_low[:,1], 'o-')
        plt.plot(self.vmin_logeta_band_up[:,0], self.vmin_logeta_band_up[:, 1], 'o-')
        self.PlotOptimum(ylim_percentage = (1.2, 0.8), plot_close = F, plot_show = T)
        
        file = output_file_tail + "_FoxBand_low.dat"
        print(file)  # write to file
        np.savetxt(file, self.vmin_logeta_band_low)
        file = output_file_tail + "_FoxBand_up.dat"
        print(file)  # write to file
        np.savetxt(file, self.vmin_logeta_band_up)
        
        return
