"""
Calculate rates with various fitting functions
"""

import numpy as np
from scipy.special import eval_chebyt


R = 8.314


def single_arrhenius(a_par, n_par, ea_par,
                     t_ref, temp):
    """ calc value with single arrhenius function
    """
    rate_constants = a_par * ((temp / t_ref)**n_par) * np.exp(-ea_par/(R*temp))
    return rate_constants


def double_arrhenius(a_par1, n_par1, ea_par1,
                     a_par2, n_par2, ea_par2,
                     t_ref, temp):
    """ calc value with single arrhenius function
    """
    rate_constants = (
        a_par1 * ((temp / t_ref)**n_par1) * np.exp(-ea_par1/(R*temp)) +
        a_par2 * ((temp / t_ref)**n_par2) * np.exp(-ea_par2/(R*temp))
    )
    return rate_constants


def arrhenius(params, t_ref, temp):
    """ simplified function whcih will call single or double arrhenis
        based on the number of params that are passed in
        only does single and double fits.
        params must be [a1, n1, ea1] or [a1, n2, ea1, a2, n2, ea2]
    """
    assert len(params) in (3, 6)

    if len(params) == 3:
        rate_constants = single_arrhenius(
            params[0], params[1], params[2],
            t_ref, temp)
    else:
        rate_constants = double_arrhenius(
            params[0], params[1], params[2],
            params[3], params[4], params[5],
            t_ref, temp)

    return rate_constants


def lindemann(highp_ks, lowp_ks, pressures, temps):
    """ calculate pressure-dependence constants according to Lindemann
        model; no value for high
    """
    rate_constant_dct = {}
    for pressure in pressures:
        rate_constant_dct[pressure] = lindemann_rate_constants(
            highp_ks, lowp_ks, pressure, temps)

    return rate_constant_dct


def lindemann_rate_constants(highp_ks, lowp_ks, pressure, temps):
    """ calculate pressure-dependence constants according to Lindemann
        model
    """
    # Calculate the pr term
    pr_term = _pr_term(highp_ks, lowp_ks, pressure, temps)

    # Calculate Lindemann rate constants
    rate_constants = highp_ks * (pr_term / (1.0 + pr_term))

    return rate_constants


def troe(highp_ks, lowp_ks, pressures, temps,
         alpha, ts1, ts3, ts2=None):
    """ calculate pressure-dependence constants according to Troe
        model; no value for high
    """
    rate_constant_dct = {}
    for pressure in pressures:
        rate_constant_dct[pressure] = troe_rate_constants(
            highp_ks, lowp_ks, pressure, temps,
            alpha, ts1, ts3, ts2)

    return rate_constant_dct


def troe_rate_constants(highp_ks, lowp_ks, pressure, temp,
                        alpha, ts1, ts3, ts2=None):
    """ calculate pressure-dependence constants according to Troe
        model
    """
    # Calculate the pr term and broadening factor
    pr_term = _pr_term(highp_ks, lowp_ks, pressure, temp)
    f_term = _f_broadening_term(alpha, pr_term, temp, ts1, ts3, ts2)

    # Calculate Troe rate constants
    rate_constants = highp_ks * (pr_term / (1.0 + pr_term)) * f_term

    return rate_constants


def plog(plog_dct, t_ref, pressures, temps):
    """ calculate the rate constant using a dictionary of plog params
    """
    rate_constant_dct = {}
    for pressure in pressures:
        rate_constant_dct[pressure] = plog_rate_constants(
            plog_dct, t_ref, pressure, temps)

    return rate_constant_dct


def plog_rate_constants(plog_dct, t_ref, pressure, temps):
    """ calculate the rate constant using a dictionary of plog params
    """
    plog_pressures = plog_dct.keys()
    plog_pressures.remove('high')

    # Check if pressure is in plog dct; use plog pressure for numerical stab
    pressure_defined = False
    for plog_pressure in plog_pressures:
        if np.isclose(pressure, plog_pressure, atol=1.0e-3):
            pressure_defined = True
            plog_params = plog_dct[plog_pressure]

    # If pressure equals value use, arrhenius expression
    if pressure_defined:
        rate_constants = arrhenius(plog_params, t_ref, temps)
    # Find which two PLOG pressures our pressure of interest sits between
    else:
        for i, _ in enumerate(plog_pressures):
            if i != len(plog_pressures):
                if plog_pressures[i] < pressure < plog_pressures[i+1]:
                    plow = plog_pressures[i]
                    phigh = plog_pressures[i+1]
                    plow_params = plog_dct[plow]
                    phigh_params = plog_dct[phigh]
                    break

        klow = arrhenius(plow_params, t_ref, temps)
        khigh = arrhenius(phigh_params, t_ref, temps)
        pres_term = (
            (np.log10(pressure) - np.log10(plow)) /
            (np.log10(phigh) - np.log10(plow))
        )
        logk = (
            np.log10(klow) +
            ((np.log10(khigh) - np.log10(klow)) * pres_term)
        )

        rate_constants = 10**(logk)

    return rate_constants


def chebyshev(alpha, tmin, tmax, pmin, pmax, pressures, temps):
    """ computes the rate constants using the chebyshev polynomials
    """
    rate_constant_dct = {}
    for pressure in pressures:
        rate_constant_dct[pressure] = chebyshev_rate_constants(
            temps, pressure, alpha, tmin, tmax, pmin, pmax)

    return rate_constant_dct


def chebyshev_rate_constants(temps, pressure, alpha, tmin, tmax, pmin, pmax):
    """ computes the rate constants using the chebyshev polynomials
    """
    alpha_nrows, alpha_ncols = alpha.shape

    rate_constants = np.zeros(len(temps))
    for i, temp in enumerate(temps):
        ctemp = (
            (2.0 * temp**(-1) - tmin**(-1) - tmax**(-1)) /
            (tmax**(-1) - tmin**(-1))
        )
        cpress = (
            (2.0 * np.log10(pressure) - np.log10(pmin) - np.log10(pmax)) /
            (np.log10(pmax) - np.log10(pmin))
        )

        logk = 0.0
        for j in range(alpha_nrows):
            for k in range(alpha_ncols):
                logk += (
                    alpha[j][k] *
                    eval_chebyt(j, ctemp) *
                    eval_chebyt(k, cpress)
                )

        rate_constants[i] = 10**(logk)

    return rate_constants


def _pr_term(highp_rateks, lowp_rateks, pressure, temp):
    """ calculate the corrective pr term used for Lindemann and Troe
        pressure-dependent forms
    """
    pr_term = (lowp_rateks / highp_rateks) * (pressure / (R * temp))
    return pr_term


def _f_broadening_term(alpha, pr_term, temp, ts1, ts3, ts2):
    """ calculate the F broadening factor used for Troe
        pressure-dependent forms
    """

    # Calculate Fcent term
    f_cent = ((1.0 - alpha) * np.exp(-temp / ts3) +
              alpha * np.exp(-temp / ts1))
    if ts2 is not None:
        f_cent += np.exp(-ts2 / temp)

    # Calculate the Log F term
    c_val = -0.4 - 0.67 * np.log10(f_cent)
    n_val = 0.75 - 1.27 * np.log10(f_cent)
    d_val = 0.14
    val = ((np.log10(pr_term) + c_val) /
           (n_val - d_val * (np.log10(pr_term) + c_val)))**2
    logf = (1.0 + val)**(-1) * np.log10(f_cent)

    # Calculate F broadening term
    f_term = 10**(logf)

    return f_term
