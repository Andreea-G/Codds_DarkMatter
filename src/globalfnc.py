# -*- coding: utf-8 -*-
"""
Created on Fri Nov 21 02:49:53 2014

@author: Andreea
"""

from __future__ import division
import numpy as np
from math import erf
import math
pi = np.pi

T = True
F = False


fermiGeV = 1./0.1973269602 #Natural[GeV femto Meter]
kilogram = 1e-9/1.782661758e-36
SpeedOfLight = 299792.458 #km/s
AtomicMassUnit = 0.931494028
ProtonMass = 1.00727646677 * AtomicMassUnit
mPhiRef = 1000.
rho = 0.3
conversion_factor = rho * SpeedOfLight**2 * 1e5 * 3600 * 24
        
v0bar = default_v0bar = 220.
vobs = default_vobs = default_v0bar + 12.
vesc = default_vesc = 533.

MaximumGapLimit_exper = ["superCDMS", \
        "LUX2013zero", "LUX2013one", "LUX2013three", "LUX2013five", "LUX2013many", \
        "XENON10", "CDMSlite2013CoGeNTQ", "CDMSSi2012"]
GaussianLimit_exper = ["KIMS2012"]
DAMARegion_exper = ["DAMA2010NaSmRebinned", "DAMA2010ISmRebinned"]

def FileNameTail(fp, fn, mPhi):
    if mPhi == 1000.:
        mPhi_string = ""
    else:
        mPhi_string = "_mPhi" + str(math.trunc(mPhi))
    fnfp = fn/fp
    fnfp_string = "_fnfp"
    if fnfp == 1.:
        fnfp_string = "_fnfp1"
    elif abs(fnfp) < 1:
        if fnfp < 0:
            fnfp_string += "_neg0" + str(math.trunc(round(10 * abs(fnfp)))) # TODO: int(.)
        else:
            fnfp_string += "0" + str(math.trunc(round(10 * abs(fnfp))))
    else:
        if fnfp < 0:
            fnfp_string += "_neg" + str(math.trunc(round(abs(fnfp))))
        else:
            fnfp_string += str(math.trunc(round(abs(fnfp))))
    return mPhi_string + fnfp_string
      
def OutputDirectory(output_main_dir, scattering_type, mPhi, delta):
    out_dir = output_main_dir
    if mPhi == 1000.:
        out_dir += "Contact"
    else:
        out_dir += "LongRange"
    out_dir += scattering_type + "_delta_"
    if delta < 0:
        out_dir += "neg"
    out_dir += str(math.trunc(abs(delta))) + "/"
    return out_dir

def Gaussian(x, mu, sigma):
    return np.exp(-(x-mu)**2 / (2 * sigma**2)) / (np.sqrt(2 * pi) * sigma)

def GPoisson(x, nu, sigma):
    eps = 1.e-4
    n = 1
    add = nu * np.exp(-(x-1.)**2 / (2 * sigma**2))
    summation = 0
    nfact = 1 #factorial
    while add > eps * (summation+add):
        summation += add
        n += 1
        nfact *= n
        add = 1. * nu**n / nfact * np.exp(-(x-n)**2 / (2. * n * sigma**2)) / np.sqrt(n)
    result = summation * np.exp(-nu) / np.sqrt(2 * np.pi) / sigma
#    print("GPoisson: ", result)
    return result

def HelmFF(ER, A, mT):
    q = np.sqrt(2e-6 * mT * ER)
    s = 0.9
    ha = 0.52
    hc = 1.23 * A**(1./3) - 0.6 #half-density radius
    cpa = 7./3. * (pi * ha)**2
    rsq = hc**2 + cpa
    r1tmp = rsq - 5. * s**2
    r1 = list(map(lambda r,s: np.sqrt(r) if r > 0 else np.sqrt(s), r1tmp, rsq))
    x = np.abs(q * r1 * fermiGeV)
    y = q * s * fermiGeV
    f = np.array(list(map(lambda i: 3.0 * (np.sin(i) - i * np.cos(i))/i**3 if i > 5.e-8 \
        else 1. - i**2 * (-0.1 + i**2 * (1./200. + i**2 * (-1./15120. + i**2/1330560.))), x)))
    return f**2 * np.exp(-y**2)

def GaussianFFSD(ER, A, mT):
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
              'SD44' : FFSD_options,
}

def VMin(ER, mT, mx, delta):
    muT = mx * mT / (mx + mT)
    return SpeedOfLight * 1.e-3 / np.sqrt(2. * ER * mT) * abs(delta + ER * mT / muT)
    
def ERecoilBranch(vmin, mT, mx, delta, sign):
    muT = mx * mT /(mx + mT)
    return 1.e6 / SpeedOfLight**2 * muT**2 * vmin**2 / (2.*mT) * \
        (1. + sign * np.sqrt(1. - 2.*delta / (muT * vmin**2) * SpeedOfLight**2 * 1.e-6))**2

def eta0Maxwellian(vmin, vobs, v0bar, vesc):
    x = vmin/v0bar
    y = vobs/v0bar
    z = vesc/v0bar
    erfz = erf(z)
    sqrt_pi = np.sqrt(pi)
    exp_z_sq = np.exp(-z**2)
    exp_z_sq_z = np.exp(-z**2 * z)
    eta = list(map(lambda i: -2. * exp_z_sq / sqrt_pi - erf(i-y) / (2.*y) + erf(i+y) / (2.*y) \
        if i + y <= z \
        else exp_z_sq * (i - y - z) / (sqrt_pi * y) - erf(i-y) / (2.*y) + erfz / (2.*y) \
        if i - y <= z < i + y \
        else 0, x))
    return eta / (-2. * exp_z_sq_z / sqrt_pi + erfz) / v0bar

def eta1Maxwellian(vmin, vobs, v0bar, vesc):
    x = vmin/v0bar
    y = vobs/v0bar
    z = vesc/v0bar
    delta_v = 15.   # 30 * cos(60 * pi / 180)
    erfz = erf(z)
    sqrt_pi = np.sqrt(pi)
    exp_z_sq = np.exp(-z**2)
    exp_z_sq_z = np.exp(-z**2 * z)
    erf_z = erf(z)
    eta = list(map(lambda i: (np.exp(-(i+y)**2) + np.exp(-(i-y)**2)) / (sqrt_pi * y) + \
        (erf(i-y) - erf(i+y)) / (2 * y**2) \
        if i + y <= z \
        else exp_z_sq * (-i + z) / (sqrt_pi * y**2) + np.exp(-(i-y)**2) / (sqrt_pi * y) + \
        (erf(i-y) - erf_z) / (2 * y**2) \
        if i - y <= z < i + y \
        else 0, x))
    return eta / (-2. * exp_z_sq_z / sqrt_pi + erfz) * delta_v / v0bar**2

def MaximumGapC0scaled(x, mu_over_x):
    if mu_over_x < 1.:
        return 1.
    elif 1. <= mu_over_x < 2.:
        return 1. - np.exp(-x) * (1. + mu_over_x * x - x)
    else:
        l = np.array([ ((k - mu_over_x) * x)**(k-1) * np.exp(-k * x) / math.factorial(k) * \
            (x * (mu_over_x - k) + k) for k in range(math.trunc(np.floor(mu_over_x)))])
        return 1. - l.sum()
        
