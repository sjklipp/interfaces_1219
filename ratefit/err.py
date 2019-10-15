"""
Calculate rates with various fitting functions
"""

import numpy as np


def calc_sse_and_mae(calc_ks, fit_ks):
    """ (1) get the sum of square error (SSE) useful when determining
            which double plog routine will be used to initialize
            the nonlinear solver
        (2) also get the mean absolute error (MAE), which is written
            to the plog file
    """

    # Only run if there are more than 2 rate constants
    sse = 0.0
    abs_err = []
    if len(calc_ks) > 2:
        for calc_k, fit_k in zip(calc_ks, fit_ks):
            sse += (np.log(calc_k) - np.log(fit_k))**2.0
            abs_err.append(np.abs((calc_k - fit_k) / calc_k))
        abs_err = np.array(abs_err, dtype=np.float64)
        mean_abs_err = np.mean(abs_err) * 100.0
        max_abs_err = np.max(abs_err) * 100.0
    else:
        sse = None
        mean_abs_err = None
        max_abs_err = None

    return sse, mean_abs_err, max_abs_err


def assess_pressure_dependence(pdep_dct, temps, temp_compare,
                               tolerance=20.0, plow=None, phigh=None):
    """ Assess how much the rate constants change from
        a low-pressure to high-pressure regime
    """
    # Get a list of the sorted pressures
    pressures = [pressure for pressure in pdep_dct
                 if pressure != 'high']
    pressures.sort()

    # Set the lowest pressure if not specified by user
    if plow is None:
        plow = min(pressures)
    if phigh is None:
        phigh = max(pressures)

    # Get idxs of rate constants corresponding to the temps for comparison
    # Get idxs corresponding to k(T) vals in each k(T, P) pdep_dct entry
    temps.sort()
    temp_idxs = [np.where(temps == temp)[0][0] for temp in temp_compare]

    # See changes
    is_pressure_dependent = False
    for idx in temp_idxs:
        ktp_low = pdep_dct[plow][idx]
        ktp_high = pdep_dct[phigh][idx]
        ktp_dif = (abs(ktp_low - ktp_high) / ktp_low) * 100.0
        if ktp_dif > tolerance:
            is_pressure_dependent = True

    return is_pressure_dependent
