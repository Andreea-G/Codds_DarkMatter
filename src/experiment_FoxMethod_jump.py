# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 00:47:37 2015

@author: Andreea
"""

#TODO! This only works for CDMSSi! 

from experiment_HaloIndep import *
import interp_uniform as unif
#import interp as unif
from interp import interp1d
#from scipy.interpolate import interp1d
from scipy.optimize import brentq, minimize
from basinhopping import *
import matplotlib.pyplot as plt
import os   # for speaking
#from multiprocessing import Pool
import parallel_map as par

DEBUG = F
DEBUG_FULL = F
USE_BASINHOPPING = T
ALLOW_MOVE = T
#METHOD = 'COBYLA'
METHOD = 'SLSQP'

class ConstraintsFunction(object):
    def __init__(self, *args): #TODO don't know why I can't put vmin_max = vmin_max in the arguments
        self.vminStar = args[0]
        self.logetaStar = args[1]
        self.vminStar_index = args[2]
        self.vmin_max = 2000
    def __call__(self, x):
        '''  0 -  8: bounds: 3 * (x.size/2) constraints = 9 for x.size/2 = 3
             9 - 12: sorted array: 2 * (x.size/2 - 1) constraints = 4 for x.size/2 = 3
            13 - 15: vminStar_index: x.size/2 constraints = 3 for x.size/2 = 3
            16 - 18: vminStar and logetaStar: x.size/2 constraints = 3 for x.size/2 = 3 '''
        constraints = np.concatenate([x[:x.size/2], self.vmin_max - x[:x.size/2], -x[x.size/2:], \
            np.diff(x[:x.size/2]), np.diff(-x[x.size/2:]), \
            self.vminStar - x[:self.vminStar_index], x[self.vminStar_index : x.size/2] - self.vminStar, \
            (x[:x.size/2] - self.vminStar) * (-x[x.size/2:] + self.logetaStar)])
        is_not_close = np.logical_not(np.isclose(constraints, np.zeros_like(constraints), atol = 1e-5))
        is_not_close[:3 * (x.size/2)] = True
        constr = np.where(is_not_close, constraints, np.abs(constraints))
#                print("~~~ constraints = ", constraints)
#                print("is_not_close = ", is_not_close)
#                print("constr = ", constr)
#                print("vminStar = ", vminStar)
#                print("logetaStar = ", logetaStar)
#                print("vminStar_index = ", vminStar_index)
#                print("~~~~~")
        if np.any(np.isnan(constr)):
            raise ValueError
        return constr

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
        ''' Computes response tables
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
        file =  output_file_tail + "_VminSortedList.dat"
        print(file)  
        np.savetxt(file, self.vmin_sorted_list)

        if delta == 0:
            branches = [1]
        else:
            branches = [1, -1]
        self.vmin_linspace = np.linspace(vmin_min, vmin_max, (vmin_max - vmin_min)/vmin_step + 1)
        self.diff_response_tab = np.zeros((self.ERecoilList.size, 1))
        self.response_tab = np.zeros(1)
        self.curly_H_tab = np.zeros((self.ERecoilList.size, 1))
        self.xi_tab = np.zeros(1)
        xi = 0
        vmin_prev = 0
        for vmin in self.vmin_linspace:
            print("vmin = ", vmin)
            diff_resp_list = np.zeros((1,3))
            resp = 0
            curly_H = np.zeros((1,3))
            for sign in branches:
                (ER, qER, const_factor) = self.ConstFactor(vmin, mx, fp, fn, delta, sign)
                v_delta = 0 #TODO! generalize to endothermic
                diff_resp_list += np.array([self.DifferentialResponse(Eee, qER, const_factor) for Eee in self.ERecoilList])
                resp += integrate.quad(self.DifferentialResponse, self.Ethreshold, self.Emaximum, \
                    args=(qER, const_factor), epsrel = PRECISSION, epsabs = 0)[0]
                curly_H += np.array([[integrate.quad(self.DifferentialResponse_Full, v_delta, vmin, \
                    args=(Eee, mx, fp, fn, delta, sign), epsrel = PRECISSION, epsabs = 0)[0] for Eee in self.ERecoilList]])
            xi += self.Exposure * self.IntegratedResponse_Other(vmin_prev, vmin, self.Ethreshold, self.Emaximum, mx, fp, fn, delta)
            vmin_prev = vmin
            self.diff_response_tab = np.append(self.diff_response_tab, diff_resp_list.transpose(), axis = 1)
            self.response_tab = np.append(self.response_tab, [resp], axis = 0)
            self.curly_H_tab = np.append(self.curly_H_tab, curly_H.transpose(), axis = 1)
            self.xi_tab = np.append(self.xi_tab, [xi], axis = 0)
        self.vmin_linspace = np.insert(self.vmin_linspace, 0., 0)
        file = output_file_tail + "_VminLinspace.dat"
        print(file)
        np.savetxt(file, self.vmin_linspace)
        file = output_file_tail + "_DiffRespTable.dat"
        print(file)
        np.savetxt(file, self.diff_response_tab)
        file = output_file_tail + "_RespTable.dat"
        print(file)
        np.savetxt(file, self.response_tab)
        file = output_file_tail + "_CurlyHTable.dat"
        print(file)
        np.savetxt(file, self.curly_H_tab)
        file = output_file_tail + "_XiTable.dat"
        print(file)
        np.savetxt(file, self.xi_tab)
        os.system("say Finished response tables.")
        return
    
    def PlotTable(self, interpolation, dimension = 0, xlim = None, ylim = None, 
                  title = None, plot_close = True, plot_show = True, show_zero_axis = False):
        if plot_close:
            plt.close()
        if dimension == 0:   # only one function
            plt.plot(self.vmin_linspace, np.array([interpolation(v) for v in self.vmin_linspace]))
        elif dimension == 1: # list of interpolated functions for each energy in self.ERecoilList
            for i in range(self.ERecoilList.size):
                plt.plot(self.vmin_linspace, np.array([interpolation[i](v) for v in self.vmin_linspace]))
        else:
            print("Wrong dimension")
            raise TypeError
        if show_zero_axis:
            plt.plot(self.vmin_linspace, np.zeros(self.vmin_linspace.size))
        if xlim != None:
            plt.xlim(xlim)
        if ylim != None:
            plt.ylim(ylim)
        if title != None:
            plt.title(title)
        if plot_show:
            plt.show()
    
    def ImportResponseTables(self, output_file_tail, plot = True):
        file =  output_file_tail + "_VminSortedList.dat"
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
        file = output_file_tail + "_CurlyHTable.dat"
        with open(file, 'r') as f_handle:
            self.curly_H_tab = np.loadtxt(f_handle)
        file = output_file_tail + "_XiTable.dat"
        with open(file, 'r') as f_handle:
            self.xi_tab = np.loadtxt(f_handle)
        self.diff_response_interp = np.array([unif.interp1d(self.vmin_linspace, dr) for dr in self.diff_response_tab])
        self.response_interp = unif.interp1d(self.vmin_linspace, self.response_tab)
        self.curly_H_interp = np.array([unif.interp1d(self.vmin_linspace, h) for h in self.curly_H_tab])
        
        if plot:
            self.PlotTable(self.diff_response_interp, dimension = 1)
            self.PlotTable(self.response_interp, dimension = 0)
            self.PlotTable(self.curly_H_interp, dimension = 1, title = 'Curly H')
        return

    def VminIntegratedResponseTable(self, vmin_list):
        return np.array([[integrate.quad(self.diff_response_interp[i], vmin_list[a], vmin_list[a + 1], epsrel = PRECISSION, epsabs = 0)[0] \
            for a in range(vmin_list.size - 1)] for i in range(self.ERecoilList.size)])

    def IntegratedResponseTable(self, vmin_list):
        return np.array([integrate.quad(self.response_interp, vmin_list[a], vmin_list[a + 1], epsrel = PRECISSION, epsabs = 0)[0] \
            for a in range(vmin_list.size - 1)])

    def MinusLogLikelihood(self, vars_list, vminStar = None, logetaStar = None, vminStar_index = None):
        if vminStar == None:
            vmin_list_w0 = vars_list[: vars_list.size/2]
            logeta_list = vars_list[vars_list.size/2 :]
        else:
            vmin_list_w0 = np.insert(vars_list[: vars_list.size/2], vminStar_index, vminStar)
            logeta_list = np.insert(vars_list[vars_list.size/2 :], vminStar_index, logetaStar)
        vmin_list_w0 = np.insert(vmin_list_w0, 0, 0)
        vmin_resp_integr = self.VminIntegratedResponseTable(vmin_list_w0)
        resp_integr = self.IntegratedResponseTable(vmin_list_w0)
        mu_i = self.Exposure * np.dot(vmin_resp_integr, 10**logeta_list)
        Nsignal = self.Exposure * np.dot(10**logeta_list, resp_integr)
        if vminStar == None:
            self.gamma_i = (self.mu_BKG_i + mu_i) / self.Exposure
        result = self.NBKG + Nsignal - np.log(self.mu_BKG_i + mu_i).sum()
        if np.any(self.mu_BKG_i + mu_i < 0):
            raise ValueError
        return result

    def _MinusLogLikelihood(self, vars_list, constr_func = None, vminStar = None, logetaStar = None, vminStar_index = None):
        constraints = constr_func(vars_list)
        constr_not_valid = constraints < 0
        if DEBUG:
            print("*** vars_list = ", repr(vars_list))
        if DEBUG_FULL:
            print("vminStar = ", vminStar)
            print("logetaStar = ", logetaStar)
            print("constraints = ", repr(constraints))
            print("constr_not_valid = ", repr(constr_not_valid))
        try:
            return self.MinusLogLikelihood(vars_list, vminStar = vminStar, logetaStar = logetaStar, vminStar_index = vminStar_index)
        except:
            if np.any(constr_not_valid):
                constr_list = constraints[constr_not_valid]
                if DEBUG_FULL: 
                    print("Constraints not valid!!********************************************************************************")
                    print("constr sum = ", -constr_list.sum())
                return min(max(-constr_list.sum(), 0.001) * 1e6, 1e6)
            else:
                print("What just happened???????????????????????????????????????????????????????????????????????????????")
                return self.MinusLogLikelihood(vars_list, vminStar = vminStar, logetaStar = logetaStar, vminStar_index = vminStar_index)

    def OptimalLikelihood(self, output_file_tail, logeta_guess = -24):
        self.ImportResponseTables(output_file_tail, plot = False)
        vars_guess = np.append(self.vmin_sorted_list, logeta_guess * np.ones(self.vmin_sorted_list.size))
        print("vars_guess = ", vars_guess)
        vmin_max = self.vmin_linspace[-1]
        def constr_func(x, vmin_max = vmin_max): 
            ''' 0 -  8: bounds: 3 * (x.size/2) constraints = 9 for x.size/2 = 3
                9 - 12: sorted array: 2 * (x.size/2 - 1) constraints = 4 for x.size/2 = 3 '''
            constraints = np.concatenate([x[:x.size/2], vmin_max - x[:x.size/2], -x[x.size/2:], \
                np.diff(x[:x.size/2]), np.diff(-x[x.size/2:])])
            is_not_close = np.logical_not(np.isclose(constraints, np.zeros_like(constraints), atol = 1e-5))
            is_not_close[:3 * (x.size/2)] = T
            constr = np.where(is_not_close, constraints, np.abs(constraints))
            if DEBUG:
                print("***constr = ", repr(constr))
                print("tf = ", repr(constr < 0))
            return constr
        constr = ({'type': 'ineq', 'fun': constr_func})
        
        np.random.seed(0)
        if USE_BASINHOPPING:
            minimizer_kwargs = {"constraints": constr, "args": (constr_func,)}
            optimum_log_likelihood = basinhopping(self._MinusLogLikelihood, vars_guess, \
                minimizer_kwargs = minimizer_kwargs, niter = 10, stepsize = 0.1)
        else:
            optimum_log_likelihood = minimize(self._MinusLogLikelihood, vars_guess, args = (constr_func,), constraints = constr)

        print(optimum_log_likelihood)
        print("MinusLogLikelihood = ", self.MinusLogLikelihood(optimum_log_likelihood.x))
        print("vars_guess = ", repr(vars_guess))
        file = output_file_tail + "_GloballyOptimalLikelihood.dat"
        print(file)  
        np.savetxt(file, np.append([optimum_log_likelihood.fun], optimum_log_likelihood.x))
        os.system("say 'Finished finding optimum'")
        return 

    def ImportOptimalLikelihood(self, output_file_tail, plot = False):
        self.ImportResponseTables(output_file_tail, plot = False)
        file = output_file_tail + "_GloballyOptimalLikelihood.dat"
        with open(file, 'r') as f_handle:
            optimal_result = np.loadtxt(f_handle)
        self.optimal_logL = optimal_result[0]
        self.optimal_vmin = optimal_result[1 : optimal_result.size/2 + 1]
        self.optimal_logeta = optimal_result[optimal_result.size/2 + 1 :]
        print("optimal result = ", optimal_result)

        if plot:
            self.MinusLogLikelihood(optimal_result[1:]) # to get self.gamma_i
            self.xi_interp = unif.interp1d(self.vmin_linspace, self.xi_tab)
            self.h_sum_tab = np.sum([self.curly_H_tab[i] / self.gamma_i[i] for i in range(self.optimal_vmin.size)], axis = 0)
            self.q_tab = 2 * (self.xi_tab - self.h_sum_tab)
            self.h_sum_interp = unif.interp1d(self.vmin_linspace, self.h_sum_tab)
            self.q_interp = unif.interp1d(self.vmin_linspace, self.q_tab)
            self.PlotTable(self.xi_interp, dimension = 0, plot_show = False)
            self.PlotTable(self.h_sum_interp, dimension = 0, xlim = [0, 2000], ylim = [-2e24, 2e24], title = 'Xi, H_sum', plot_close = False)
            self.PlotTable(self.q_interp, dimension = 0, xlim = [0, 2000], ylim = [-2e24, 2e24], title = 'q', show_zero_axis = True)
        return

    def PlotStepFunction(self, vmin_list, logeta_list, xlim_percentage = (0., 1.1), ylim_percentage = (1.01, 0.99), \
        plot_close = True, plot_show = True, mark = 'o', color = None):
        if plot_close: 
            plt.close()
        print(vmin_list)
        print(logeta_list)
        x = np.insert(vmin_list, 0, 0)
        y = np.insert(logeta_list, 0, logeta_list[0])
        plt.step(x, y)
        if color != None:
            plt.plot(x, y, mark, color = color)
        else:
            plt.plot(x, y, mark)
#        plt.xlim([vmin_list[0] * xlim_percentage[0], vmin_list[-1] * xlim_percentage[1]])
        plt.xlim([0, 1000])
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
            plot_close = False, plot_show = False, mark = 'x', color = 'k')
        plt.plot(vminStar, logetaStar, '*')
        if plot_show:
            plt.show()
        return
        
    def _ConstrainedOptimalLikelihood(self, vminStar, logetaStar, vminStar_index, plot = False):
#        vmin_max = self.vmin_linspace[-1]
        if DEBUG:
            print("vminStar_index =", vminStar_index)
        vmin_guess = np.concatenate([np.minimum(self.optimal_vmin[:vminStar_index], np.ones(vminStar_index) * vminStar), \
            np.maximum(self.optimal_vmin[vminStar_index:], np.ones(self.optimal_vmin.size - vminStar_index) * vminStar)])
        vmin_guess = np.append(np.array([self.optimal_vmin[ind] if self.optimal_vmin[ind] < vminStar \
                                        else vminStar * (1 - 0.001*(vminStar_index - ind)) \
                                        for ind in range(vminStar_index)]), \
                               np.array([self.optimal_vmin[ind] if self.optimal_vmin[ind] > vminStar \
                                       else vminStar * (1 + 0.001*(ind - vminStar_index - 1)) \
                                       for ind in range(vminStar_index, self.optimal_vmin.size)]))
        logeta_guess = self.optimal_logeta
        logeta_guess = np.concatenate([np.maximum(logeta_guess[:vminStar_index], np.ones(vminStar_index)*logetaStar), \
            np.minimum(logeta_guess[vminStar_index:], np.ones(logeta_guess.size - vminStar_index)*logetaStar)])
        vars_guess = np.append(vmin_guess, logeta_guess)

        constr_func = ConstraintsFunction(vminStar, logetaStar, vminStar_index)
        constr = ({'type': 'ineq', 'fun': constr_func})
        args = (constr_func, vminStar, logetaStar, vminStar_index)
        
        sol_not_found = True
        attempts = 3
        np.random.seed(1)
        random_variation = 1e-5

        if USE_BASINHOPPING:
            class TakeStep(object):
                def __init__(self, stepsize = 0.1):
                    pass
                    self.stepsize = stepsize
                def __call__(self, x):
                    x[:x.size/2] += np.random.uniform(-5. * self.stepsize, 5. * self.stepsize, x[x.size/2:].shape)
                    x[x.size/2:] += np.random.uniform(-self.stepsize, self.stepsize, x[x.size/2:].shape)
                    return x
            take_step = TakeStep()

            class AdaptiveKwargs(object):
                def __init__(self, kwargs, random_variation = random_variation):
                    self.kwargs = kwargs
                    self.random_variation = random_variation
                def __call__(self):
                    new_kwargs = {}
                    constr_func_args = (self.kwargs['args'][1] * (1 + self.random_variation * np.random.uniform(-1, 1)), \
                                        self.kwargs['args'][2] * (1 + self.random_variation * np.random.uniform(-1, 1)), \
                                        self.kwargs['args'][3])
                    constr_func = ConstraintsFunction(*constr_func_args)
                    new_kwargs['args'] = (constr_func,) + constr_func_args
                    new_kwargs['constraints'] = ({'type': 'ineq', 'fun': constr_func})
                    if 'method' in self.kwargs:
                        new_kwargs['method'] = self.kwargs['method']
                    return new_kwargs
    
            minimizer_kwargs = {"constraints": constr, "args": args, "method": METHOD}
            adapt_kwargs = AdaptiveKwargs(minimizer_kwargs, random_variation)
        
        while sol_not_found and attempts > 0:
            try:
                if USE_BASINHOPPING:
                    constr_optimum_log_likelihood = basinhopping(self._MinusLogLikelihood, vars_guess, \
                        minimizer_kwargs = minimizer_kwargs, niter = 3, take_step = take_step, adapt_kwargs = adapt_kwargs, stepsize = 0.2)
                else:
                    constr_optimum_log_likelihood = minimize(self._MinusLogLikelihood, vars_guess, \
                        args = args, constraints = constr, method = METHOD)
                constraints = constr_func(constr_optimum_log_likelihood.x)
                is_not_close = np.logical_not(np.isclose(constraints, np.zeros_like(constraints)))
                constr_not_valid = np.logical_and(constraints < 0, is_not_close)
                sol_not_found = np.any(constr_not_valid)
            except ValueError:
                sol_not_found = True
                pass
            
            attempts -= 1
            args = (constr_func, vminStar * (1 + random_variation * np.random.uniform(-1, 1)), \
                    logetaStar * (1 + random_variation * np.random.uniform(-1, 1)), vminStar_index)
            if USE_BASINHOPPING:
                minimizer_kwargs = {"constraints": constr, "args": args}
                
            if DEBUG and sol_not_found:
                print(attempts, " attempts left! ####################################################################################################")
#                os.system("say Error"); os.system("say " + str(attempts) + " attempts left")
                print("sol_not_found = ", sol_not_found)
        if sol_not_found:
            if DEBUG:
                print("ValueError: sol not found")
            raise ValueError
        
        if DEBUG:
            print(constr_optimum_log_likelihood)
            print("kwargs = ", constr_optimum_log_likelihood.minimizer.kwargs)
            print("args = ", constr_optimum_log_likelihood.minimizer.kwargs['args'])
            print("optimum_logL = ", self.optimal_logL)
            print("constraints = ", repr(constraints))
            print("constr_not_valid = ", repr(constr_not_valid))
            print("vars_guess = ", repr(vars_guess))
            print("optimum_logL = ", self.optimal_logL)
            print("vminStar_index = ", vminStar_index)
            
        return constr_optimum_log_likelihood

    def ConstrainedOptimalLikelihood(self, vminStar, logetaStar, plot = False):
        vminStar_index = 0
        while vminStar_index < self.optimal_vmin.size and vminStar > self.optimal_vmin[vminStar_index]:
            vminStar_index += 1
        
        try:
            constr_optimum_log_likelihood = self._ConstrainedOptimalLikelihood(vminStar, logetaStar, vminStar_index, plot = plot)
        except ValueError:
            optim_logL = 10**6
            pass
        else:
            optim_logL = constr_optimum_log_likelihood.fun
            original_optimum = constr_optimum_log_likelihood
        
        vminStar_index_original = vminStar_index
        index = vminStar_index
#        move_left = ALLOW_MOVE and vminStar_index > 0
        while ALLOW_MOVE and index > 0:
            try:
                index -= 1
                new_optimum = self._ConstrainedOptimalLikelihood(vminStar, logetaStar, index, plot = plot)
            except ValueError:
                pass
#                if vminStar_index <= 0:
#                    move_left = False
            else:
                if new_optimum.fun < optim_logL:
                    os.system("say Moved left")
                    print("#################################################################################################################################################")
                    print("Moved left, index is now ", index)
                    vminStar_index = index
                    constr_optimum_log_likelihood = new_optimum
                    optim_logL = constr_optimum_log_likelihood.fun
#                if vminStar_index <= 0:
#                    move_left = False
        index = vminStar_index_original
#        move_right = ALLOW_MOVE and vminStar_index < self.optimal_vmin.size
        while ALLOW_MOVE and index < self.optimal_vmin.size:
            try:
                index += 1
                new_optimum = self._ConstrainedOptimalLikelihood(vminStar, logetaStar, index, plot = plot)
            except ValueError:
                pass
#                if vminStar_index >= self.optimal_vmin.size:
#                    move_right = False
            else:
                if new_optimum.fun < optim_logL:
                    os.system("say Moved right")
                    print("#################################################################################################################################################")
                    print("Moved right, index is now ", index)
                    vminStar_index = index
                    constr_optimum_log_likelihood = new_optimum
                    optim_logL = constr_optimum_log_likelihood.fun
#                if vminStar_index >= self.optimal_vmin.size:
#                    move_right = False
        if optim_logL == 10**6:
#            os.system("say No solution was found")
            raise ValueError
        
        self.constr_optimal_logl = constr_optimum_log_likelihood.fun
        vars_result = constr_optimum_log_likelihood.x

        self.constr_optimal_vmin = vars_result[: vars_result.size/2]
        self.constr_optimal_logeta = vars_result[vars_result.size/2:]
#        print(self.constr_optimal_vmin)
#        print(self.constr_optimal_logeta)
        
        if plot:
            print("vminStar = ", vminStar)
            print("logetaStar = ", logetaStar)
            print("vminStar_index = ", vminStar_index)
            try:
                print("original: ", original_optimum)
                print("new: ", constr_optimum_log_likelihood)
                constr_func = ConstraintsFunction(vminStar, logetaStar, vminStar_index)
                constraints = constr_func(constr_optimum_log_likelihood.x)
                is_not_close = np.logical_not(np.isclose(constraints, np.zeros_like(constraints)))
                constr_not_valid = np.logical_and(constraints < 0, is_not_close)
                sol_not_found = np.any(constr_not_valid)
                print("constraints = ", constraints)
                print("is_not_close = ", is_not_close)
                print("constr_not_valid = ", constr_not_valid)
                print("sol_not_found = ", sol_not_found)
            except:
                print("Error")
                pass
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
        # Find a better way to deal with this.
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
        g = lambda x, x0, s1, s2: g1(x, x0, s1) + g2(x, x0, s2)
        
        s_list = np.array([[s, sc]] + [[sc,sc]] * (numx0 - 2) + [[sc, s]])
        g_total = lambda x, sign = 1, x0 = x0_list, s_list = s_list: \
            np.array([sign * g(x, x0_list[i], s_list[i,0], s_list[i,1]) \
            for i in range(x0_list.size)]).prod(axis = 0)
        g_lin = g_total(x_lin)
        
        xT_guess = (x0_list[:-1] + x0_list[1:]) / 2
        bounds = np.array([(x0_list[i], x0_list[i + 1]) for i in range(x0_list.size - 1)])
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
            plt.ylim([min(-s * sc**(numx0 - 1), np.min(y_turns)), max(s * sc**(numx0 - 1), np.max(y_turns))])
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
            print("sampling_size = ", self.vmin_logeta_sampling_table.shape)
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
    
    def GetLikelihoodTable(self, kwargs):
        ''' Prints to file lists of the form [logetaStar_ij, logL_ij] needed for 1D interpolation, 
        where i is the index corresponding to vminStar_i and j is the index for each logetaStar. Each file corresponds to a different index i.
            Here only one file is written for a specific vminStar.
            Input:
            - kwargs containing a dictionary of the form:
            {'index': index, 'output_file_tail': output_file_tail}
            where index is the index of vminStar and output_file_tail is part of the file name where the table will be printed.
            Output:
            - none
        '''
        index = kwargs['index']
        output_file_tail = kwargs['output_file_tail']
        print("index = ", index)
        print("output_file_tail = ", output_file_tail)
        vminStar = self.vmin_logeta_sampling_table[index, 0, 0]
        logetaStar_list = self.vmin_logeta_sampling_table[index, :, 1][24:25]
        print("vminStar = ", vminStar)
        table = np.empty((0,2))
        for logetaStar in logetaStar_list:
            try:
                constr_opt = self.ConstrainedOptimalLikelihood(vminStar, logetaStar, plot = True)
            except:
                print("error")
                os.system("say Error")
                pass
            else:
                print("logetaStar = ", logetaStar, "; constr_opt = ", constr_opt)
                table = np.append(table, [[logetaStar, constr_opt]], axis = 0)
#                table = np.append(table, [logetaStar])
        print("vminStar = ", vminStar, "; table = ", table)
        if True:
            temp_file = output_file_tail + '_' + str(index) + '_LogetaStarLogLikelihoodList_junk.dat'
            print(temp_file)
            np.savetxt(temp_file, table)
        return
    
    def LogLikelihoodList(self, output_file_tail, processes = None):
        ''' Loops thorugh the list of all vminStar and calls GetLikelihoodTable, which will print the likelihood tables to files.
            Input:
            - output_file_tail: part of the files name
            - processes: number of processes for parallel programming
            Output:
            - none
        '''
        index_list = range(24, 25)#, self.vmin_logeta_sampling_table.shape[0])
        kwargs = ({'index': index, \
            'output_file_tail': output_file_tail} for index in index_list)
#        pool = Pool()
#        pool.map(GetTable, kwargs)
        par.parmap(self.GetLikelihoodTable, kwargs, processes)
        return
    
    def _logL_interp(vars_list, constraints):
        constr_not_valid = constraints(vars_list)[:-1] < 0
        if np.any(constr_not_valid):
            constr_list = constraints(vars_list)[constr_not_valid]
            return -constr_list.sum() * 10**2
        return logL_interp(vars_list)

    def FoxBand(self, output_file_tail, delta_logL, interpolation_order, multiplot = True, plot = True):
        print("self.vmin_sampling_list = ", self.vmin_sampling_list)
        self.vmin_logeta_band_low = []
        self.vmin_logeta_band_up = []
        vmin_last_step = self.optimal_vmin[-1]
        if multiplot:
            plt.close()
        for index in range(self.vmin_sampling_list.size):
            print("index = ", index)
            print("vmin = ", self.vmin_sampling_list[index])
            logeta_optim = self.OptimumStepFunction(min(self.vmin_sampling_list[index], vmin_last_step))
            file = output_file_tail + "_" + str(index) + "_LogetaStarLogLikelihoodList.dat"
            try:
                with open(file, 'r') as f_handle:
                    table = np.loadtxt(f_handle)
            except:
                continue
            x = table[:, 0]   # this is logeta
            y = table[:, 1]   # this is logL
            logL_interp = interp1d(x, y)
            
            def _logL_interp(vars_list, constraints):
                constr_not_valid = constraints(vars_list)[:-1] < 0
                if np.any(constr_not_valid):
                    constr_list = constraints(vars_list)[constr_not_valid]
                    return -constr_list.sum() * 1e2
                return logL_interp(vars_list)
                        
            print(self.optimal_logL - delta_logL)
            print(np.array([table[0, 0]]), " ", table[-1, 0])
            print(logeta_optim)
            constr_func = lambda logeta, logeta_min = np.array([table[0, 0]]), \
                logeta_max = np.array([table[-1, 0]]): \
                np.concatenate([logeta - logeta_min, logeta_max - logeta])
            constr = ({'type': 'ineq', 'fun': constr_func})
            try:
                logeta_minimLogL = minimize(_logL_interp, np.array([logeta_optim]), args = (constr_func,), constraints = constr).x
            except ValueError:
                print("ValueError at logeta_minimLogL")
                logeta_minimLogL = self.optimal_logL
                pass
            print(logeta_minimLogL)
            
            print("x = ", x)
            print("y = ", y)
            if multiplot:
                plt.close()
                plt.plot(x, y, 'o-')
                plt.plot(x, (self.optimal_logL + 1) * np.ones_like(y))
                plt.plot(x, (self.optimal_logL + 4) * np.ones_like(y))
                plt.title("index = " + str(index) + ", v_min = " + str(self.vmin_sampling_list[index]) + "km/s")
                plt.xlim(x[0], x[-1])
                plt.ylim(-5, 20)
                plt.show()

            try:
                if y[0] > self.optimal_logL + delta_logL and logeta_minimLogL < self.optimal_logL + delta_logL:
                    self.vmin_logeta_band_low += [[self.vmin_sampling_list[index], \
                        brentq(lambda logeta: logL_interp(logeta) - self.optimal_logL - delta_logL, \
                            table[0, 0], logeta_minimLogL)]]
                if y[-1] > self.optimal_logL + delta_logL and logeta_minimLogL < self.optimal_logL + delta_logL:
                    print("a, b: ", logeta_optim, " ", table[-1, 0])
                    print("f: ", logL_interp(logeta_optim) - self.optimal_logL - delta_logL, " ", \
                        logL_interp(table[-1, 0]) - self.optimal_logL - delta_logL)
                    self.vmin_logeta_band_up += [[self.vmin_sampling_list[index], \
                        brentq(lambda logeta: logL_interp(logeta) - self.optimal_logL - delta_logL, \
                            logeta_minimLogL, table[-1, 0])]]
            except ValueError:
                plt.close()
                plt.plot(x, (self.optimal_logL + 1) * np.ones_like(y))
                plt.plot(x, (self.optimal_logL + 4) * np.ones_like(y))
                plt.title("v_min = " + str(self.vmin_sampling_list[index]) + "km/s")
                plt.xlim(x[0], x[-1])
                plt.ylim([-5,20])
                plt.plot(x, y, 'o-', color="r")
                plt.plot(logeta_optim, logL_interp(logeta_optim), '*')
                plt.plot(logeta_optim, self.optimal_logL, '*')
                print("ValueError")
                plt.show()
#                raise
                pass
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
        print(file)  
        np.savetxt(file, self.vmin_logeta_band_low)
        file = output_file_tail + "_FoxBand_up.dat"
        print(file)  
        np.savetxt(file, self.vmin_logeta_band_up)
        
        return

    def ImportFoxBand(self, output_file_tail):
        file = output_file_tail + "_FoxBand_low.dat"
        with open(file, 'r') as f_handle:
            self.vmin_logeta_band_low = np.loadtxt(f_handle)
        file = output_file_tail + "_FoxBand_up.dat"
        with open(file, 'r') as f_handle:
            self.vmin_logeta_band_up = np.loadtxt(f_handle)
        plt.close()
        plt.plot(self.vmin_logeta_band_low[:,0], self.vmin_logeta_band_low[:,1], 'o-')
        plt.plot(self.vmin_logeta_band_up[:,0], self.vmin_logeta_band_up[:, 1], 'o-')
        self.PlotOptimum(ylim_percentage = (1.2, 0.8), plot_close = F, plot_show = T)
        return

