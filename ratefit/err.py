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
        mean_avg_err = np.mean(abs_err)*100.0
        max_avg_err = np.max(abs_err)*100.0
    else:
        sse = None
        mean_avg_err = None
        max_avg_err = None

    return sse, mean_avg_err, max_avg_err
