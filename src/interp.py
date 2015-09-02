"""
Copyright (c) 2015 Andreea Georgescu

Created on Thu Nov 20 22:52:11 2014

This file is a shortened version of the scipy.interpolate.interp1d algorithm,
restricted to linear interpolation and lacks most of the checks and robustness
of the scipy version. It is faster, but must be run correctly
(undefined behavior otherwise).


This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.


SciPy license:

Copyright (c) 2001, 2002 Enthought, Inc.
All rights reserved.

Copyright (c) 2003-2013 SciPy Developers.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

    Redistributions of source code must retain the above copyright notice,
    this list of conditions and the following disclaimer.
    Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the documentation
    and/or other materials provided with the distribution.
    Neither the name of Enthought nor the names of the SciPy Developers may
    be used to endorse or promote products derived from this software without
    specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS”
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

from __future__ import print_function
from __future__ import division
import numpy as np


class interp1d:
    """
    This is based off of scipy.interpolate.interp1d, but is restricted to linear
    interpolation and lacks most of the checks and robustness of the scipy version.
    It is faster, but must be run correctly (undefined behavior otherwise).
    Input for initialization:
        x, y: ndarray
            Two numpy.arrays x and y sorted in increasing order by x.
    """
    def __init__(self, x, y, kind="linear"):
        self.x = x
        self.y = y
        self.len = len(self.x)
        if kind == "linear":
            self._call = self._call_linear
        else:
            raise NotImplementedError("This kind is unsupported!")

    def __call__(self, x_new):
        return self._evaluate(x_new)

    def _check_bounds(self, x_new):
        below_bounds = x_new < self.x[0]
        above_bounds = x_new > self.x[-1]

        if below_bounds.any():
            raise ValueError("A value in x_new is below the interpolation range.")
        if above_bounds.any():
            raise ValueError("A value in x_new is above the interpolation range.")
        return

    def _evaluate(self, x_new):
        self._check_bounds(x_new)
        y_new = self._call(x_new)
        return y_new

    def _call_linear(self, x_new):
        x_new_index = np.searchsorted(self.x, x_new)
        x_new_index = x_new_index.clip(1, self.len-1).astype(int)
        lo = x_new_index - 1
        hi = x_new_index
        x_lo = self.x[lo]
        x_hi = self.x[hi]
        y_lo = self.y[lo]
        y_hi = self.y[hi]
        slope = (y_hi - y_lo) / (x_hi - x_lo)
        y_new = slope*(x_new - x_lo) + y_lo
        return y_new
