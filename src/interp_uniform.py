# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 16:20:52 2015

@author: Andreea
"""
from __future__ import print_function
from __future__ import division
import numpy as np
from math import trunc


class interp1d:
    '''
    This is specific to the case of uniform spacing in x!
    This is based off of scipy.interpolate.interp1d, but lacks most of the checks and
    robustness of the scipy version.
    It is a bit faster, but must be run correctly (undefined behavior otherwise).
    Input for initialization:
        two numpy.arrays x and y sorted in increasing order by x.
    '''
    def __init__(self, x, y, kind="linear", fill_value=np.nan):
        self.x = x
        self.y = y
        self.len = len(self.x)
        self.spacing = x[1] - x[0]
        self.fill_value = fill_value
        if kind == "linear":
            self._call = self._call_linear
        else:
            raise NotImplementedError("This kind is unsupported!")

    def __call__(self, x_new):
        return self._evaluate(x_new)

    def _check_bounds(self, x_new):
        below_bounds = x_new < self.x[0]
        above_bounds = x_new > self.x[-1]
        if below_bounds:
            raise ValueError("A value in x_new is below the interpolation range.")
        if above_bounds:
            raise ValueError("A value in x_new is above the interpolation range.")
        out_of_bounds = np.logical_or(below_bounds, above_bounds)
        return out_of_bounds

    def _evaluate(self, x_new):
        out_of_bounds = self._check_bounds(x_new)
        if out_of_bounds:
            y_new = self.fill_value
        else:
            y_new = self._call(x_new)
        return y_new

    def _call_linear(self, x_new):
        if x_new == self.x[0]:
            return self.y[0]
        x_new_index = trunc((x_new - self.x[0])/self.spacing)
        if x_new != self.x[x_new_index]:
            x_new_index += 1
        lo = x_new_index - 1
        hi = x_new_index
        slope = (self.y[hi] - self.y[lo]) / (self.x[hi] - self.x[lo])
        y_new = slope*(x_new - self.x[lo]) + self.y[lo]
        return y_new
