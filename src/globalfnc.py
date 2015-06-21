# -*- coding: utf-8 -*-
"""
Created on Fri Nov 21 02:49:53 2014

@author: Andreea
"""

from __future__ import division
from __future__ import absolute_import
from __future__ import print_function
import numpy as np
import math
from scipy.special import erf, erfinv
from scipy.stats import chi2
pi = np.pi

T = True    # short-hand notation
F = False

PRECISSION = 1.e-3
# Unit conversions
fermiGeV = 1./0.1973269602  # Natural[GeV femto Meter]
kilogram = 1e-9/1.782661758e-36  # kg in GeV  (units: [GeV/kg])
SpeedOfLight = 299792.458  # km/s
AtomicMassUnit = 0.931494028
ProtonMass = 1.00727646677 * AtomicMassUnit
mPhiRef = 1000.
rho = 0.3  # GeV/cm**3
conversion_factor = rho * SpeedOfLight**2 * 1e5 * 3600 * 24
ConfidenceLevel = 0.9

v0bar = default_v0bar = 220.
vobs = default_vobs = default_v0bar + 12.
vesc = default_vesc = 533.

if False:   # alternative velocities; kept for reference
    v0bar = 230 - 3 * 24.4
    vobs = v0bar + 12
    vesc = 544 - 3 * 39

''' List of experiment names corresponding to each type of statistical analysis
'''
MaximumGapLimit_exper = ["superCDMS",
                         "LUX2013zero", "LUX2013one", "LUX2013three", "LUX2013five", "LUX2013many",
                         "XENON10", "XENON100", "CDMSlite2013CoGeNTQ", "CDMSSi2012"]
GaussianLimit_exper = ["KIMS2012", "PICASSO"]
DAMARegion_exper = ["DAMA2010NaSmRebinned", "DAMA2010ISmRebinned"]
DAMALimit_exper = ["DAMA2010NaSmRebinned_TotRateLimit"]
Poisson_exper = ["SIMPLEModeStage2"]
FoxMethod_exper = ["CDMSSi2012", "CDMSSiGeArtif", "CDMSSiArtif"]


''' Colors for plotting
'''
Color = {"superCDMS": 'peru',
         "LUX2013zero": 'magenta', "LUX2013one": 'magenta', "LUX2013three": 'magenta',
         "LUX2013five": 'magenta', "LUX2013many": 'magenta',
         "XENON10": 'orange', "XENON100": 'royalblue',
         "CDMSlite2013CoGeNTQ": 'cyan', "CDMSSi2012": 'red',
         "KIMS2012": 'purple', "PICASSO": 'darkturquoise',
         "DAMA2010NaSmRebinned_TotRateLimit": 'green',
         "DAMA2010NaSmRebinned": 'green',
         "DAMA2010ISmRebinned": 'green',
         "SIMPLEModeStage2": 'saddlebrown'
         }


def confidence_level(sigma):
    return erf(sigma/np.sqrt(2))


def sigma_dev(CL):
    return np.sqrt(2) * erfinv(CL)


def chi_squared(dof, CL=ConfidenceLevel):
    return chi2.ppf(CL, dof)


def chi_squared1(CL=ConfidenceLevel):
    # special case for dof = 1
    return sigma_dev(CL)**2


def import_file(full_path_to_module):
    ''' Imports Python module from file.
    '''
    import os
    import sys
    directory, module_name = os.path.split(full_path_to_module)
    module_name = os.path.splitext(module_name)[0]

    path = list(sys.path)
    sys.path.insert(0, directory)
    try:
        module = __import__(module_name)
        return module
    finally:
        sys.path[:] = path  # restore


def FileNameTail(fp, fn, mPhi):
    ''' Gives a file name-tail that is added to each output file, to distinguish
    between different input parameters.
    '''
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
            fnfp_string += "_neg0" + str(math.trunc(round(10 * abs(fnfp))))  # TODO: int(.)
        else:
            fnfp_string += "0" + str(math.trunc(round(10 * abs(fnfp))))
    else:
        if fnfp < 0:
            fnfp_string += "_neg" + str(math.trunc(round(abs(fnfp))))
        else:
            fnfp_string += str(math.trunc(round(abs(fnfp))))
    return mPhi_string + fnfp_string


def OutputDirectory(output_main_dir, scattering_type, mPhi, delta):
    ''' Gives the name of the output directory for the given input parameters.
    '''
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


def Output_file_name(exper_name, scattering_type, mPhi, mx, fp, fn, delta, HALO_DEP,
                     filename_tail, OUTPUT_MAIN_DIR, quenching):
    ''' Gives the name of the output file name for the given input parameters.
    '''
    output_dir = OutputDirectory(OUTPUT_MAIN_DIR, scattering_type, mPhi, delta)
    output_file_no_extension = "./" + output_dir + "UpperLimit_" + exper_name

    if HALO_DEP:
        output_file_no_extension += "_mxsigma"
        if vesc != default_vesc:
            output_file_no_extension += "_vesc" \
                + str(math.trunc(round(vesc)))
        if vobs != default_vobs:
            output_file_no_extension += "_vobs" \
                + str(math.trunc(round(vobs)))
    else:
        output_file_no_extension += "_mx_" + str(mx) + "GeV"

    output_file_no_extension += FileNameTail(fp, fn, mPhi) + filename_tail

    if quenching is not None:
        output_file_no_extension += "_q" + str(quenching)
    print(output_file_no_extension)
    return output_file_no_extension


def Gaussian(x, mu, sigma):
    ''' Gaussian resolution function.
    '''
    return np.exp(-(x-mu)**2 / (2 * sigma**2)) / (np.sqrt(2 * pi) * sigma)


def GPoisson(x, nu, sigma):
    ''' Resolution function for the Xe experiment; it is a combination of
    Poisson and Gaussian resolution.
        nu is the expected # of events.
    '''
    eps = 1.e-4
    n = 1
    add = nu * np.exp(-(x-1.)**2 / (2 * sigma**2))
    summation = 0
    nfact = 1  # factorial
    while add > eps * (summation+add):
        summation += add
        n += 1
        nfact *= n
        add = 1. * nu**n / nfact * np.exp(-(x-n)**2 / (2. * n * sigma**2)) / np.sqrt(n)
    result = summation * np.exp(-nu) / np.sqrt(2 * np.pi) / sigma
#    print("GPoisson: ", result)
    return result


def HelmFF(ER, A, mT):
    ''' Helm Form Factor. See http://arxiv.org/abs/hep-ph/0608035.
    Input:
        ER: recoil energy
        A: target mass number
        mT: target nuclide mass
    '''
    q = np.sqrt(2e-6 * mT * ER)
    s = 0.9
    ha = 0.52
    hc = 1.23 * A**(1./3) - 0.6  # half-density radius
    cpa = 7./3. * (pi * ha)**2
    rsq = hc**2 + cpa
    r1tmp = rsq - 5. * s**2
    r1 = list(map(lambda r, s: np.sqrt(r) if r > 0 else np.sqrt(s), r1tmp, rsq))
    x = np.abs(q * r1 * fermiGeV)
    y = q * s * fermiGeV
    f = np.array(list(map(lambda i: 3.0 * (np.sin(i) - i * np.cos(i))/i**3 if i > 5.e-8
                 else 1. - i**2 * (-0.1 + i**2 * (1./200. + i**2 * (-1./15120. + i**2/1330560.))), x)))
    return f**2 * np.exp(-y**2)


def GaussianFFSD(ER, A, mT):
    ''' Gaussian Form Factor for spin-dependent interactions.
    Inpput:
        ER: recoil energy
        A: target atomic mass
        mT: target nuclide mass
    '''
    q = np.sqrt(2e-6 * mT * ER)
    R = 0.92 * A**(1./3) + 2.68 - 0.78 * np.sqrt((A**(1./3) - 3.8)**2 + 0.2)
    x = np.abs(q * R * fermiGeV)
    return np.exp(-x**2 / 4.)

''' Form factor options are used to select the correct FF depending on the type of
interaction (spin-independent, spin-dependent with axial-vector or pseudo-scalar coupling).
'''
FFSI_options = {'HelmFF': HelmFF,
                }
FFSD_options = {'GaussianFFSD': GaussianFFSD,
                }
FF_options = {'SI': FFSI_options,
              'SDPS': FFSD_options,
              'SDAV': FFSD_options,
              }


def ERecoil_ratio(mT_1, mT_2, mx, quenching_1, quenching_2):
    mu_1 = mT_1 * mx / (mT_1 + mx)
    mu_2 = mT_2 * mx / (mT_2 + mx)
    return mu_2**2 / mu_1**2 * mT_1 / mT_2 * quenching_2 / quenching_1


def VMin(ER, mT, mx, delta):
    ''' Minimum velocity for a given recoil energy ER, target mass mT, DM mass mx,
    and mass split delta.
    '''
    muT = mx * mT / (mx + mT)
    return SpeedOfLight * 1.e-3 / np.sqrt(2. * ER * mT) * abs(delta + ER * mT / muT)


def VminDelta(mT, mx, delta):
    muT = mx * mT / (mx + mT)
    return SpeedOfLight / 500. * np.sqrt(delta / 2. / muT) if delta > 0 \
        else np.array([0] * len(mT))


def ERecoilBranch(vmin, mT, mx, delta, sign):
    ''' Recoil energy for given vmin.
    Input:
        vmin: minimum DM velocity vmin
        mT: target mass
        mx: DM mass
        delta: mass split
        sign: +1 or -1, corresponding to the upper and lower branch, respectively.
    '''
    muT = mx * mT / (mx + mT)
    return 1.e6 * muT**2 / (2.*mT) * \
        (vmin / SpeedOfLight +
         sign * np.sqrt((vmin / SpeedOfLight)**2 - 2.*delta / muT * 1.e-6))**2


def dERecoildVmin(vmin, mT, mx, delta, sign):
    ''' Derivative of recoil energy ER with respect to velocity vmin
    Input:
        vmin: minimum DM velocity vmin
        mT: target mass
        mx: DM mass
        delta: mass split
        sign: +1 or -1, corresponding to the upper and lower branch, respectively.
    '''
    muT = mx * mT / (mx + mT)
    sqrt_factor = np.sqrt(1. - 2.*delta / (muT * vmin**2) * SpeedOfLight**2 * 1.e-6)
    return sign * muT**2 * vmin / mT * (1. + sign * sqrt_factor)**2 / sqrt_factor


def eta0Maxwellian(vmin, vobs, v0bar, vesc):
    ''' Velocity integral eta0 in a Standard Halo Model with Maxwellian velocity
    distribution.
    Input:
        vmin: minumum DM velocity
        vobs: observed velocity
        v0bar: velocity dispersion
        vesc: excape velocity
    '''
    x = vmin/v0bar
    y = vobs/v0bar
    z = vesc/v0bar
    erfz = erf(z)
    sqrt_pi = np.sqrt(pi)
    exp_z_sq = np.exp(-z**2)
    exp_z_sq_z = np.exp(-z**2) * z
    eta = list(map(lambda i: -2. * exp_z_sq / sqrt_pi - erf(i-y) / (2.*y) +
                   erf(i+y) / (2.*y)
                   if i + y <= z
                   else exp_z_sq * (i - y - z) / (sqrt_pi * y) - erf(i-y) / (2.*y) +
                   erfz / (2.*y)
                   if i - y <= z < i + y
                   else 0, x))
    return eta / (-2. * exp_z_sq_z / sqrt_pi + erfz) / v0bar


def eta1Maxwellian(vmin, vobs, v0bar, vesc):
    ''' Same as eta0Maxwellian, but this is the modulation velocity integral
    eta1 = d eta0 / d vobs * delta_v.
        delta_v = v_Earth * cos(gamma) where the velocity of the Earth is
    v_Earth = 30 km/s and is inclined at an angle of gamma = 60 deg wrt the galactic plane.
    '''
    x = vmin/v0bar
    y = vobs/v0bar
    z = vesc/v0bar
    delta_v = 15.   # 30 * cos(60 * pi / 180)
    erfz = erf(z)
    sqrt_pi = np.sqrt(pi)
    exp_z_sq = np.exp(-z**2)
    exp_z_sq_z = np.exp(-z**2) * z
    erf_z = erf(z)
    eta = list(map(lambda i: (np.exp(-(i+y)**2) + np.exp(-(i-y)**2)) / (sqrt_pi * y) +
                   (erf(i-y) - erf(i+y)) / (2 * y**2)
                   if i + y <= z
                   else exp_z_sq * (-i + z) / (sqrt_pi * y**2) +
                   np.exp(-(i-y)**2) / (sqrt_pi * y) + (erf(i-y) - erf_z) / (2 * y**2)
                   if i - y <= z < i + y
                   else 0, x))
    return eta / (-2. * exp_z_sq_z / sqrt_pi + erfz) * delta_v / v0bar**2


def MaximumGapC0scaled(x, mu_over_x):
    ''' Scaled probability C0 of the maximum gap size being smaller than a particular
    value of x, where x is the size of the maximum gap in the random experiment
    and mu is the total expected number of events.
    Based on Maximum Gap Method (Eq 2 of PHYSICAL REVIEW D 66, 032005, 2002)
    '''
    if mu_over_x < 1.:
        return 1.
    elif 1. <= mu_over_x < 2.:
        return 1. - np.exp(-x) * (1. + mu_over_x * x - x)
    else:
        l = np.array([((k - mu_over_x) * x)**(k-1) * np.exp(-k * x) / math.factorial(k) *
                     (x * (mu_over_x - k) + k)
                     for k in range(1, math.trunc(np.floor(mu_over_x)))])
        return 1. - l.sum()


if __name__ == "__main__":
    vmin = np.linspace(1, 1000, 1000)
    mx = 1.3
    delta = -200
    mT = np.array([65.134, 66.995, 67.9278, 68.8571, 70.7203])
    ER_plus = np.array([ERecoilBranch(vmin, mT_i, mx, delta, 1) for mT_i in mT])
    ER_minus = np.array([ERecoilBranch(vmin, mT_i, mx, delta, -1) for mT_i in mT])
#    print("vmin =", vmin)
#    print("ER_plus =", ER_plus)
#    print("ER_minus =", ER_minus)
    import matplotlib.pyplot as plt
    plt.close()
    for ER in ER_plus:
        plt.plot(vmin, ER)
    for ER in ER_minus:
        plt.plot(vmin, ER)
    plt.plot(vmin, 2.7 * np.ones_like(vmin))
    plt.plot(vmin, 3 * np.ones_like(vmin))
    ER = np.linspace(2, 6, 100)
    plt.plot(490 * np.ones_like(ER), ER)
    plt.plot(495 * np.ones_like(ER), ER)
#    plt.xlim([450, 550])
#    plt.ylim([2, 5])
    plt.show()
