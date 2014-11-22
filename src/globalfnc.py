# -*- coding: utf-8 -*-
"""
Created on Fri Nov 21 02:49:53 2014

@author: Andreea
"""


fermiGeV = 1./0.1973269602 #Natural[GeV femto Meter]
kilogram = 1e-9/1.782661758e-36
SpeedOfLight = 299792.458 #km/s
AtomicMassUnit = 0.931494028
ProtonMass = 1.00727646677 * AtomicMassUnit
mPhiRef = 1.

vobs = 232.
v0bar = 220.
vesc = 533.

import numpy as np

def HelmFF(ER, A, mT):
    print 'HelmFF'
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
    print 'GaussianFFSD'
    q = np.sqrt(2e-6 * mT * ER)
    R = 0.92 * A**(1./3) + 2.68 - 0.78 * np.sqrt((A**(1./3) - 3.8)**2 + 0.2)
    x = np.abs(q * R * fermiGeV)
    return np.exp(-x**2 / 4)
    

FFSI_options = {'HelmFF' : HelmFF,
}
FFSD_options = {'GaussianFFSD' : GaussianFFSD,
}
FF_options = {'SI' : FFSI_options,
              'SD66' : FFSD_options,
}