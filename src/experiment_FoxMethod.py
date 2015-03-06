# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 00:47:37 2015

@author: Andreea
"""

#TODO! This only works for CDMSSi! 

from experiment_HaloIndep import *
from interp import interp1d
#from scipy.interpolate import interp1d
from scipy.optimize import minimize, brentq
import matplotlib.pyplot as plt

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

    def ImportResponseTables(self, output_file_tail):
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
#        plt.close()
#        plt.plot(self.vmin_linspace, self.diff_response_tab[0])
#        plt.plot(self.vmin_linspace, self.diff_response_tab[1])
#        plt.plot(self.vmin_linspace, self.diff_response_tab[2])
#        plt.plot(self.vmin_linspace, self.response_tab)
#        plt.show()
        self.diff_response_interp = np.array([interp1d(self.vmin_linspace, dr) for dr in self.diff_response_tab])
        self.response_interp = interp1d(self.vmin_linspace, self.response_tab)
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
        mu_i = self.Exposure * np.dot(self.VminIntegratedResponseTable(vmin_list_w0), 10**logeta_list)
        Nsignal = self.Exposure * np.dot(10**logeta_list, self.IntegratedResponseTable(vmin_list_w0))
        result = self.NBKG + Nsignal - np.log10(self.mu_BKG_i + mu_i).sum()
#        print("result = ", result)
        return result

    def OptimalLikelihood(self, output_file_tail, logeta_guess = -25.):
        self.ImportResponseTables(output_file_tail)
        vars_guess = np.append(self.vmin_sorted_list, logeta_guess * np.ones(self.vmin_sorted_list.size))
        print("vars_guess = ", vars_guess)
        bounds = np.array([(0, 1000)] * self.vmin_sorted_list.size + [(None, 0)] * self.vmin_sorted_list.size)
        constr_func = lambda x:  np.append(np.diff(x[:x.size/2]), np.diff(-x[x.size/2:]))
        constr = ({'type': 'ineq', 'fun': constr_func})
        optimum_log_likelihood = minimize(self.MinusLogLikelihood, vars_guess, bounds = bounds, constraints = constr)
        print(optimum_log_likelihood)
        file = output_file_tail + "_GloballyOptimalLikelihood.dat"
        print(file)  # write to file
        np.savetxt(file, np.append([optimum_log_likelihood.fun], optimum_log_likelihood.x))
        return 

    def ImportOptimalLikelihood(self, output_file_tail):
        self.ImportResponseTables(output_file_tail)
        file = output_file_tail + "_GloballyOptimalLikelihood.dat"
        with open(file, 'r') as f_handle:
            optimal_result = np.loadtxt(f_handle)
        self.optimal_logl = optimal_result[0]
        self.optimal_vmin = optimal_result[1 : optimal_result.size/2 + 1]
        self.optimal_logeta = optimal_result[optimal_result.size/2 + 1 :]
        return

    def PlotStepFunction(self, vmin_list, logeta_list, plt_close = True, plt_show = True):
        if plt_close: 
            plt.close()
        print(vmin_list)
        print(logeta_list)
        plt.step(np.insert(vmin_list, 0, 0),np.insert(logeta_list, 0, logeta_list[0]))
        plt.xlim([vmin_list[0] * 0.9, vmin_list[-1] * 1.1])
        plt.ylim([logeta_list[-1] * 1.02, logeta_list[0] * 0.99])
        plt.xlim([0, vmin_list[-1] * 1.1])
        plt.ylim([logeta_list[-1] * 1.02, logeta_list[0] * 0.99])
        if plt_show:
            plt.show()
        return

    def PlotOptimum(self):
        self.PlotStepFunction(self.optimal_vmin, self.optimal_logeta)
        return

    def PlotConstrainedOptimum(self, vminStar, logetaStar, vminStar_index):
        self.PlotStepFunction(self.optimal_vmin, self.optimal_logeta, plt_show = False)
        self.PlotStepFunction(np.insert(self.constr_optimal_vmin, vminStar_index, vminStar), \
            np.insert(self.constr_optimal_logeta, vminStar_index, logetaStar), \
            plt_close = False)
        return

    def ConstrainedOptimalLikelihood(self, vminStar, logetaStar, logeta_guess = -25., plot = False):
        vars_guess = np.append(self.optimal_vmin, self.optimal_logeta)
        print("vars_guess = ", vars_guess)
        size = self.optimal_vmin.size
        
        vminStar_index = 0
        while vminStar > self.optimal_vmin[vminStar_index]:
            vminStar_index += 1
        print("index = ", vminStar_index)
        
        bounds = np.array([(0, 1000)] * size + [(None, 0)] * size)
        
        constr_func = lambda x, vminStar = vminStar, logetaStar = logetaStar : \
            np.append(np.append(np.append(np.append(np.append(np.diff(x[:x.size/2]), np.diff(-x[x.size/2:])), \
            vminStar - x[:vminStar_index]), x[vminStar_index : x.size/2] - vminStar), \
            x[x.size/2 : x.size/2 + vminStar_index] - logetaStar), logetaStar - x[x.size/2 + vminStar_index :])
        constr = ({'type': 'ineq', 'fun': constr_func})
        constr_optimum_log_likelihood = minimize(self.MinusLogLikelihood, vars_guess, \
            args = (vminStar, logetaStar, vminStar_index), bounds = bounds, constraints = constr)
        print(constr_optimum_log_likelihood)
        
        self.constr_optimal_logl = constr_optimum_log_likelihood.fun
        vars_result = constr_optimum_log_likelihood.x
        print("vars_result = ", vars_result)
        self.constr_optimal_vmin = vars_result[: vars_result.size/2]
        self.constr_optimal_logeta = vars_result[vars_result.size/2:]
        print(self.constr_optimal_vmin)
        print(self.constr_optimal_logeta)
        if plot:
            self.PlotConstrainedOptimum(vminStar, logetaStar, vminStar_index)
        return self.constr_optimal_logl

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
        self.ImportSamplingTable(output_file_tail)
        self.logL_list_table = np.empty((0, 2))
        for index in range(vmin_logeta_sampling_table.size):
            vminStar = vmin_logeta_sampling_table[index, 0, 0]
            logetaStar_list = self.vmin_logeta_sampling_table[index, All, 1]
            table = np.array([[logetaStar, ConstrainedOptimalLikelihood(vminStar, logetaStar)] \
                for logetaStar in logetaStar_list])
            self.logL_list_table = np.append(self.constr_optimal_logl, table, axis = 0)
        
        file = output_file_tail + "_LogetaStarLogLikelihoodList.dat"
        print(file)  # write to file
        np.savetxt(file, self.logL_list_table)
        return
            
    def VminSamplingList(self, output_file_tail, vmin_min, vmin_max, vmin_num_steps, \
        steepnessVmin = 1.5, steepnessVminCenter = 2.5, plot = True):
        self.ImportOptimalLikelihood(output_file_tail)
        xmin = vmin_min
        xmax = vmin_max
        x_num_steps = vmin_num_steps
        s = steepnessVmin
        sc = steepnessVminCenter
        
        x_lin = np.linspace(xmin, xmax, 1000)
        x0_list = self.optimal_vmin
        numx0 = x0_list.size
        
#        print("x = ", x_lin)
        print("x0 = ", x0_list)
        UnitStep = lambda x: (np.sign(x) + 1) / 2
        g1 = lambda x, x0, s0, xmin = xmin: np.log10(UnitStep(x - x0) + \
            UnitStep(x0 - x) * (x0 - xmin)/ (x + 10**s0 * (-x + x0) - xmin))
        g2 = lambda x, x0, s0, xmax = xmax: np.log10(UnitStep(x0 - x) + \
            UnitStep(x - x0) * (x + 10**s0 * (-x + x0) - xmax)/ (x0 - xmax))
#        print("g1 = ", np.array([[x, g1(x, x0_list[0], s)] for x in x_lin]))
#        print("g2 = ", np.array([[x, g2(x, x0_list[0], s)] for x in x_lin]))
        
        g = lambda x, x0, s1, s2: g1(x, x0, s1) + g2(x, x0, s2) * UnitStep(x - x0)
        
        s_list = np.array([[s, sc]] + [[sc,sc]] * (numx0 - 2) + [[sc, s]])
        g_total = lambda x, sign = 1, x0 = x0_list, s_list = s_list: \
            np.array([sign * g(x, x0_list[i], s_list[i,0], s_list[i,1]) \
            for i in range(x0_list.size)]).prod(axis = 0)
        
        g_lin = g_total(x_lin)
#        print("g = ", g_lin)
        
        xT_guess = (x0_list[:-1] + x0_list[1:]) / 2
        print("xT_guess = ", xT_guess)
        bounds = np.array([(x0_list[i], x0_list[i + 1]) for i in range(x0_list.size - 1)])
        print("bounds = ", bounds)
        print("gtot = ", g_total(xT_guess))
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
        y_diff_sum = y_diff.sum()
        num_steps = np.array([max(1, np.floor(x_num_steps * np.abs(yd)/y_diff_sum)) for yd in y_diff])
        y_list = np.array([np.linspace(y_turns[i], y_turns[i+1], num_steps[i]) for i in range(num_steps.size)])
        x_list = np.array([g_inverse_list(y_list[i], x_turns[i], x_turns[i+1]) for i in range(y_list.size)])
        x_list = np.concatenate(x_list)
        y_list = np.concatenate(y_list)
        x_list = x_list[np.array([x_list[i] != x_list[i+1] for i in range(x_list.size-1)] + [True])]
        y_list = y_list[np.array([y_list[i] != y_list[i+1] for i in range(y_list.size-1)] + [True])]
        print("x_list = ", x_list)        
        print("y_list = ", y_list)
        self.vmin_star_list = x_list

        if plot:
            plt.close()
            plt.plot(x_lin, g_lin)
            plt.plot(x_turns, y_turns, 'o')
            plt.plot(x_list, y_list, '*')
            plt.xlim([xmin, xmax])
            plt.ylim([-s * sc**(numx0 - 1), s * sc**(numx0 - 1)])
            plt.show()
        
        return

        











