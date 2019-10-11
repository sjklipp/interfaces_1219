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


def pdep_dependence(pdep_dct, plow=None, tol1=20, tol2=20):
    """ Assess how much the rate constants change from
        a low-pressure regime
    """
    # Get a list of the sorted pressures
    pressures = [pressure for pressure in pdep_dct
                 if pressure != 'high']
    pressures.sort()

    # Set the lowest pressure if not specified by user
    if plow is None:
        plow = min(pressures)

    # See changes
    within_tolerance = True
    for i, pressure in enumerate(pressures):
        if i != len(pressures) + 1:
            ktplow = pdep_dct[pressures[i]]
            ktphigh = pdep_dct[pressures[i+1]]
            abs_err = np.abs(ktplow - ktphigh) / ktplow
            mean_abs_err = np.mean(abs_err) * 100.0
            max_abs_err = np.max(abs_err) * 100.0
            if mean_abs_err > tol1 and max_abs_err > tol2:
                within_tolerance = False

    return within_tolerance
