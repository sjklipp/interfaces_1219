""" fit rate constants to Arrhenius expressions
"""

import numpy as np
import scipy.optimize.leastsq

# put in QCEngine gas constant
R = 8.314


def get_valid_temps_rate_constants(temps, rate_constants):
    """ this subroutine takes in a array of rate constants and
        returns the subset of this array that is positive,
        along with the corresponding Temperature array """

    # start using only the temps at which the rate constant is well defined
    valid_t, valid_k = [], []
    for temp, rate_constant in zip(temps, rate_constants):
        if rate_constant > 0.0 and min(temps) <= temp <= max(temps):
            valid_t.append(temp)
            valid_k.append(rate_constant)

    # Convert the lists to numpy arrays
    valid_t = np.array(valid_t, dtype=np.float64)
    valid_k = np.array(valid_k, dtype=np.float64)

    return valid_t, valid_k


def fit_single_arrhenius(temps, rate_constants, t_ref):
    """ this subroutine takes in a vector of rate constants and
        returns the Arrhenius parameters, as well as
        the T-range over which they were fit"""

    # obtain temperatures at which the rate constant is well defined
    valid_t, valid_k = get_valid_temps_rate_constants(temps, rate_constants)

    # consider several cases depending on the number of valid rate constants
    # no k is positive, so return all zeros
    if not valid_k:
        a_fit, n_fit, ea_fit = 0.0, 0.0, 0.0
        fit_range = [0, 0]

    # if num(k) > 0 is 1: set A = k
    elif len(valid_k) == 1:
        a_fit, n_fit, ea_fit = valid_k, 0.0, 0.0
        fit_range = [temps.index(min(valid_t)), temps.index(max(valid_t))]

    # if num(k) > 0 is 2,3: fit A and Ea
    elif (len(valid_k) == 2) or (len(valid_k) == 3):
        # Build vectors and matrices used for the fitting
        a_vec = np.ones(len(valid_t))
        ea_vec = (-1.0 / R) * valid_t
        coeff_mat = np.array([a_vec, ea_vec], dtype=np.float64)
        coeff_mat = coeff_mat.transpose()
        k_vec = np.log(valid_k)
        # Perform the least-squares fit
        theta = np.linalg.lstsq(coeff_mat, k_vec)[0]
        # Set the fitting parameters
        a_fit, n_fit, ea_fit = np.exp(theta[0]), 0.0, theta[1]
        fit_range = [temps.index(min(valid_t)), temps.index(max(valid_t))]

    # if num(k) > 0 is more than 3: fit A, n, and Ea
    elif len(valid_k) > 3:
        # Build vectors and matrices used for the fitting
        a_vec = np.ones(len(valid_t))
        n_vec = np.log(valid_t / t_ref)
        ea_vec = (-1.0 / R) * valid_t
        coeff_mat = np.array([a_vec, n_vec, ea_vec], dtype=np.float64)
        coeff_mat = coeff_mat.transpose()
        k_vec = np.log(valid_k)
        # Perform the least-squares fit
        theta = np.linalg.lstsq(coeff_mat, k_vec)[0]
        # Set the fitting parameters
        a_fit, n_fit, ea_fit = np.exp(theta[0]), theta[1], theta[2]
        fit_range = [temps.index(min(valid_t)), temps.index(max(valid_t))]

    return a_fit, n_fit, ea_fit, fit_range


def fit_double_arrhenius_dsarrfit(temps, rate_constants, t_ref):
    """ the main subroutine for obtaining the sum of two Arrhenius expressions.
        (1) Basically, the subroutine first determines whether the curve is
            bending up or down at high temperature.
        (2) Based upon the curvature, one of two possible cycles are initiated
            to provide a good starting guess for the nonlinear solver."""

    # obtain temperatures at which the rate constant is well defined
    valid_t, valid_k = get_valid_temps_rate_constants(temps, rate_constants)

    # obtain initial fitting coefficients from single fitting
    s_a, s_n, s_ea, s_fit_range = fit_single_arrhenius(
        temps, rate_constants, t_ref)

    # build initial guess vector for the double fitting
    first_guess = [(s_a / 2.0), s_n, s_ea,
                   (s_a / 2.0), s_n, s_ea]

    # obtain fitted rate constants
    fit_ks = s_a * (valid_t / t_ref)**s_n * np.exp(-s_ea / (R * valid_t))

    # calculate the sum-of-squares error
    first_sse = 0.0

    # write the input for ds arrfit and run the program
    arrfit_inp_str = write_arrfit_inp(temps, ks, locat_T)
    with open(arr_fit_inp, 'w') as input_file:
        input_file.write(arrfit_inp_str)

    # run the program
    run_arrfit(path)

    # parse the arrfit output
    with open(path, 'r') as output_file:
        out_str = output_file.read()
    arrfit_guess = parse_arrfit(out_str)

    # Check if the fitting was successful
    if len(arrfit_guess) == 6:
        best_guess = arrfit_guess
    else:
        best_guess = None

    return best_guess


def fit_double_arrhenius_python(temps, rate_constants, t_ref):
    """ perform a double Arrhenius fit with python
    """

    # obtain temperatures at which the rate constant is well defined
    valid_t, valid_k = get_valid_temps_rate_constants(temps, rate_constants)

    # obtain initial fitting coefficients from single fitting
    s_a, s_n, s_ea, s_fit_range = fit_single_arrhenius(
        temps, rate_constants, t_ref)

    # Build a guess vector
    init_guess = [(s_a / 2.0), (s_n + 0.1), s_ea,
                  (s_a / 2.0), (s_n - 0.1), s_ea]

    # Perform a new least-squares fit
    plsq = scipy.optimize.leastsq(mod_arr_residuals, init_guess,
                                  args=(valid_k, valid_t, data),
                                  ftol=1.0E-9, xtol=1.0E-9, maxfev=100000)
    [a1, n1, ea1, a2, n2, ea2] = plsq[0]

    # compute new rate constants using fitted parameters
    fit_ks = (
        a1 * (valid_t / t_ref)**n1 * np.exp(-ea1 / (R * valid_t)) +
        a2 * (valid_t / t_ref)**n2 * np.exp(-ea2 / (R * valid_t))
    )

    # Compute new SSE and assess if this is the best guess
    new_sse = 0.0
    if new_sse < first_sse:
        best_guess = [a1, n1, ea1, a2, n2, ea2]
    else:
        print("nonlinear solver was not better than single exponential!")
        best_guess = init_guess
            
    # calculate sse for the final guess 

    return best_guess, fit_range


def get_valid_temps_rate_constants(temps, rate_constants):
    """ this subroutine takes in a array of rate constants and
        returns the subset of this array that is positive,
        along with the corresponding Temperature array """

    # start using only the temps at which the rate constant is well defined
    valid_t, valid_k = [], []
    for temp, rate_constant in zip(temps, rate_constants):
        if rate_constant > 0.0 and min(temps) <= temp <= max(temps):
            valid_t.append(temp)
            valid_k.append(rate_constant)

    # Convert the lists to numpy arrays
    valid_t = np.array(valid_t, dtype=np.float64)
    valid_k = np.array(valid_k, dtype=np.float64)

    return valid_t, valid_k


def calc_sse_and_mae(local_ks, fit_ks):
    """ (1) get the sum of square error (SSE) useful when determining
            which double plog routine will be used to initialize
            the nonlinear solver
        (2) also get the mean absolute error (MAE), which is written
            to the plog file
    """

    # Only run if there are more than 2 rate constants
    sse = 0.0
    mae = []
    if len(local_ks) > 2:
        for local_k, fit_k in zip(local_ks, fit_ks):
            sse += (np.log(local_k) - np.log(fit_k))**2.0
            mae.append(np.abs((local_k - fit_k) / local_k))
        mae = np.array(mae, dtype=np.float64)
        mae = [np.mean(mae)*100.0, np.max(mae)*100.0]
    else:
        sse = 0.0
        mae = [0.0, 0.0]

    return sse, mae


def mod_arr_residuals(p, target, temp, t_ref):
    """ this subroutine computes the residual used by the nonlinear solver
        in fit_double_arrhenius_python
    """

    # compute the fitted rate constant
    k_fit1 = np.exp(
        np.log(p[0]) + p[1] * np.log((temp/t_ref)) - p[2]/(R * temp))
    k_fit2 = np.exp(
        np.log(p[3]) + p[4] * np.log((temp/t_ref)) - p[5]/(R * temp))
    k_fit = k_fit1 + k_fit2

    # calculate error
    err = np.log10(target) - np.log10(k_fit)

    return err
