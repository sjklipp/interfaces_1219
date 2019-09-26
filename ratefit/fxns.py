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
    rate_ks = a_par * ((temp / t_ref)**n_par) * np.exp(-ea_par/(R*temp))
    return rate_ks


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


def lindemann(highp_rateks, lowp_rateks, pressure, temp):
    """ calculate pressure-dependence constants according to Lindemann
        model
    """
    # Calculate the pr term
    pr_term = _pr_term(highp_rateks, lowp_rateks, pressure, temp)

    # Calculate Lindemann rate constants
    rate_constants = highp_rateks * (pr_term / (1.0 + pr_term))

    return rate_constants


def troe(highp_rateks, lowp_rateks, pressure, temp,
         alpha, t1, t3, t2=0.0):
    """ calculate pressure-dependence constants according to Troe
        model
    """
    # Calculate the pr term and broadening factor
    pr_term = _pr_term(highp_rateks, lowp_rateks, pressure, temp)
    f_term = _f_broadening_term(alpha, pr_term, temp, t1, t3, t2)

    # Calculate Troe rate constants
    rate_constants = highp_rateks * (pr_term / (1.0 + pr_term)) * f_term

    return rate_constants


def plog(ra, temps, pressure, p1, p2):
    """ calculate the rate constant at intermediate pressure
    """
    pres = (
        (np.log10(pressure) - np.log10(p1)) / 
        (np.log10(p2) - np.log10(p1))
    )
    logk = np.log10(k1) + ((np.log10(k2) - np.log10(k1)) * pres_term)
    rate_constants = 10**(logk)

    return rate_constants


# def chebyshev(temps, pressures, alpha, tmin, tmax):
#     """ computes the rate constants using the chebyshev polynomials
#     """
# 
#     logk = 0.0
#     for i, temp in enumerate(temps):
#         for j, pressure in enumerate(pressures):
#             ctemp = (
#                 ((2.0 * temp**(-1)) - tmin**(-1) - tmax**(-1)) /
#                 (tmax**(-1) - tin**(-1))
#             )
#             cpress = (
#                 ((2.0 * np.log10(pressure)) - np.log10(pmin) - np.log10(pmax)) /
#                 (np.log10(pmax) - np.log10(pmin))
#             )
#             logk += (
#                 alpha[i][j] * 
#                 eval_chebyt(i, ctemp) * 
#                 eval_chebyt(j, cpress) * 
#             )
# 
#     rate_constants = 10**(logk)

    return rate_constants


def _pr_term(highp_rateks, lowp_rateks, pressure, temp):
    """ calculate the corrective pr term used for Lindemann and Troe
        pressure-dependent forms
    """
    pr_term = (lowp_rateks / highp_rateks) * (pressure / (R * temp))
    return pr_term


def _f_broadening_term(alpha, pr, temp, t1, t3, t2):
    """ calculate the F broadening factor used for Troe
        pressure-dependent forms
    """

    # Calculate Fcent term
    f_cent = ((1.0 - alpha) * np.exp(-temp / t3) +
              alpha * np.exp(-temp / t1) +
              np.exp(-t2 / temp))

    # Calculate the Log F term
    c_val = -0.4 - 0.67 * np.log10(f_cent)
    n_val = 0.75 - 1.27 * np.log10(f_cent)
    d_val = 0.14
    val = ((np.log10(pr) + c_val) /
           (n_val - d_val * (np.log10(pr) + c_val)))**2
    logf = (1.0 + val)**(-1) * np.log10(f_cent)

    # Calculate F broadening term
    f_term = 10**(logf)

    return f_term
