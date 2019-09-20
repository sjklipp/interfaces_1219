""" fit rate constants to Arrhenius expressions
"""

# import sys
import numpy as np
from scipy.optimize import leastsq

# put in QCEngine gas constant
R = 8.314


def to_single_arrhenius(temps, rate_constants, t_ref):
    """ this subroutine takes in a vector of rate constants and
        returns the Arrhenius parameters, as well as
        the T-range over which they were fit"""

    # consider several cases depending on the number of valid rate constants
    # no k is positive, so return all zeros
    if len(rate_constants) == 0:
        a_fit, n_fit, ea_fit = 0.0, 0.0, 0.0

    # if num(k) > 0 is 1: set A = k
    elif len(rate_constants) == 1:
        a_fit, n_fit, ea_fit = rate_constants[0], 0.0, 0.0

    # if num(k) > 0 is 2,3: fit A and Ea
    elif (len(rate_constants) == 2) or (len(rate_constants) == 3):
        # Build vectors and matrices used for the fitting
        a_vec = np.ones(len(temps))
        ea_vec = (-1.0 / R) * (1.0 / temps)
        coeff_mat = np.array([a_vec, ea_vec], dtype=np.float64)
        coeff_mat = coeff_mat.transpose()
        k_vec = np.log(rate_constants)
        # Perform the least-squares fit
        theta = np.linalg.lstsq(coeff_mat, k_vec, rcond=None)[0]
        # Set the fitting parameters
        a_fit, n_fit, ea_fit = np.exp(theta[0]), 0.0, theta[1]

    # if num(k) > 0 is more than 3: fit A, n, and Ea
    elif len(rate_constants) > 3:
        # Build vectors and matrices used for the fitting
        a_vec = np.ones(len(temps))
        n_vec = np.log(temps / t_ref)
        ea_vec = (-1.0 / R) * (1.0 / temps)
        coeff_mat = np.array([a_vec, n_vec, ea_vec], dtype=np.float64)
        coeff_mat = coeff_mat.transpose()
        k_vec = np.log(rate_constants)
        # Perform the least-squares fit
        theta = np.linalg.lstsq(coeff_mat, k_vec, rcond=None)[0]
        # Set the fitting parameters
        a_fit, n_fit, ea_fit = np.exp(theta[0]), theta[1], theta[2]

    # Pack the parameters into a list
    fit_params = [a_fit, n_fit, ea_fit]

    return fit_params


def to_double_arrhenius_dsarrfit(temps, rate_constants, t_ref):
    """ the main subroutine for obtaining the sum of two Arrhenius expressions.
        (1) Basically, the subroutine first determines whether the curve is
            bending up or down at high temperature.
        (2) Based upon the curvature, one of two possible cycles are initiated
            to provide a good starting guess for the nonlinear solver."""


    # Check if the fitting was successful
    if len(arrfit_guess) == 6:
        best_guess = arrfit_guess
    else:
        best_guess = None

    return best_guess


def to_double_arrhenius_scipy(sgl_a, sgl_n, sgl_ea,
                              temps, rate_constants, t_ref):
    """ perform a double Arrhenius fit with python
    """

    # Build a guess vector
    guess_params = [(sgl_a / 2.0), (sgl_n + 0.1), sgl_ea,
                    (sgl_a / 2.0), (sgl_n - 0.1), sgl_ea]

    # Perform a new least-squares fit
    plsq = leastsq(_mod_arr_residuals, guess_params,
                   args=(rate_constants, temps, t_ref),
                   ftol=1.0E-9, xtol=1.0E-9, maxfev=100000)

    # Assign parameters

    return plsq[0]


def _mod_arr_residuals(guess_params, rate_constant, temp, t_ref):
    """ this subroutine computes the residual used by the nonlinear solver
        in fit_double_arrhenius_python
    """

    # compute the fitted rate constant
    k_fit1 = np.exp(
        np.log(guess_params[0]) +
        guess_params[1] * np.log(temp/t_ref) -
        guess_params[2]/(R * temp)
    )
    k_fit2 = np.exp(
        np.log(guess_params[3]) +
        guess_params[4] * np.log(temp/t_ref) -
        guess_params[5]/(R * temp)
    )
    k_fit = k_fit1 + k_fit2

    # calculate error
    err = np.log10(rate_constant) - np.log10(k_fit)

    return err


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
