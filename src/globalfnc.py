# -*- coding: utf-8 -*-
"""
Created on Fri Nov 21 02:49:53 2014

@author: Andreea
"""

from __future__ import division
import numpy as np
from math import erf
import math


fermiGeV = 1./0.1973269602 #Natural[GeV femto Meter]
kilogram = 1e-9/1.782661758e-36
SpeedOfLight = 299792.458 #km/s
AtomicMassUnit = 0.931494028
ProtonMass = 1.00727646677 * AtomicMassUnit
mPhiRef = 1000.

vobs = 232.
v0bar = 220.
vesc = 533.

def HelmFF(ER, A, mT):
    #print 'HelmFF'
    q = np.sqrt(2e-6 * mT * ER)
    s = 0.9
    ha = 0.52
    hc = 1.23 * A**(1./3) - 0.6 #half-density radius
    cpa = 7./3. * (np.pi * ha)**2
    rsq = hc**2 + cpa
    r1tmp = rsq - 5. * s**2
    r1 = map(lambda r,s: np.sqrt(r) if r > 0 else np.sqrt(s), r1tmp, rsq)
    x = np.abs(q * r1 * fermiGeV)
    y = q * s * fermiGeV
    f = np.array(map(lambda i: 3.0 * (np.sin(i) - i * np.cos(i))/i**3 if i > 5.e-8 \
        else 1. - i**2 * (-0.1 + i**2 * (1./200. + i**2 * (-1./15120. + i**2/1330560.))), x))
    return f**2 * np.exp(-y**2)

def GaussianFFSD(ER, A, mT):
    #print 'GaussianFFSD'
    q = np.sqrt(2e-6 * mT * ER)
    R = 0.92 * A**(1./3) + 2.68 - 0.78 * np.sqrt((A**(1./3) - 3.8)**2 + 0.2)
    x = np.abs(q * R * fermiGeV)
    return np.exp(-x**2 / 4.)
    

FFSI_options = {'HelmFF' : HelmFF,
}
FFSD_options = {'GaussianFFSD' : GaussianFFSD,
}
FF_options = {'SI' : FFSI_options,
              'SD66' : FFSD_options,
}

def VMin(ER, mT, mx, delta):
    muT = mx * mT / (mx + mT)
    return SpeedOfLight / np.sqrt(2.e6 * ER * mT) * np.abs(delta + ER * mT / muT)
    
def ERecoilBranch(vmin, mT, mx, delta, sign):
    muT = mx * mT /(mx + mT)
    return 1.e6 / SpeedOfLight**2 * muT**2 * vmin**2 / (2.*mT) * \
        (1. + sign * np.sqrt(1. - 2.*delta / (muT * vmin**2) * SpeedOfLight**2 * 1.e-6))**2

def eta0Maxwellian(vmin, vobs, v0bar, vesc):
    x = vmin/v0bar
    y = vobs/v0bar
    z = vesc/v0bar
    eta = map(lambda i: -2. * np.exp(-z**2) / np.sqrt(np.pi) - erf(i-y) / (2.*y) + erf(i+y) / (2.*y) \
        if i + y <= z \
        else np.exp(-z**2) * (i - y - z) / (np.sqrt(np.pi) * y) - erf(i-y) / (2.*y) + erf(z) / (2.*y) \
        if i - y <= z < i + y \
        else 0, x)
    return eta / (-2. * np.exp(-z**2 * z) / np.sqrt(np.pi) + erf(z)) / v0bar
        
def MaximumGapC0scaled(x, mu_over_x):
    if mu_over_x < 1.:
        return 1.
    elif 1. <= mu_over_x < 2.:
        return 1. - np.exp(-x) * (1. + mu_over_x * x - x)
    else:
        l = np.array([ ((k - mu_over_x) * x)**(k-1) * np.exp(-k * x) / math.factorial(k) * \
            (x * (mu_over_x - k) + k) for k in range(math.trunc(np.floor(mu_over_x)))])
        return 1. - l.sum()
        