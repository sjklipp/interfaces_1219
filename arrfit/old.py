""" fit rate constants to Arrhenius expressions
"""

# import sys
import numpy as np
# import scipy.optimize.leastsq

# put in QCEngine gas constant
R = 8.314


def single_arrhenius_fit(temps, rate_constants, t_ref):
    """ this subroutine takes in a vector of rate constants and
        returns the Arrhenius parameters, as well as
        the T-range over which they were fit"""

    # consider several cases depending on the number of valid rate constants
    # no k is positive, so return all zeros
    if len(rate_constants) == 0:
        a_fit, n_fit, ea_fit = 0.0, 0.0, 0.0
        fit_range = [0, 0]

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


def double_arrhenius_fit_dsarrfit(temps, rate_constants, t_ref):
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


def double_arrhenius_fit_scipy(temps, rate_constants, t_ref):
    """ perform a double Arrhenius fit with python
    """

    # obtain temperatures at which the rate constant is well defined
    valid_t, valid_k = get_valid_temps_rate_constants(temps, rate_constants)

    # obtain initial fitting coefficients from single fitting
    s_a, s_n, s_ea, s_fit_range = single_arrhenius_fit(
        temps, rate_constants, t_ref)

    # rate constants using fitted parameters
    sgl_fit_ks = single_arrhenius(s_a, s_n, s_ea,
                                  t_ref, temps)

    # Compute new SSE and assess if this is the best guess
    sgl_fit_sse = calc_sse_and_mae(valid_k, sgl_fit_ks)

    # Build a guess vector
    init_guess = [(s_a / 2.0), (s_n + 0.1), s_ea,
                  (s_a / 2.0), (s_n - 0.1), s_ea]

    # Perform a new least-squares fit
    plsq = scipy.optimize.leastsq(mod_arr_residuals, init_guess,
                                  args=(valid_k, valid_t, data),
                                  ftol=1.0E-9, xtol=1.0E-9, maxfev=100000)
    [a_par1, n_par1, ea_par1, a_par2, n_par2, ea_par2] = plsq[0]

    # compute new rate constants using fitted parameters
    dbl_fit_ks = double_arrhenius(a_par1, n_par1, ea_par1,
                                  a_par2, n_par2, ea_par2,
                                  t_ref, temps)

    # Compute new SSE and assess if this is the best guess
    dbl_fit_sse = calc_sse_and_mae(valid_k, dbl_fit_ks)

    # Assess if double fit is better
    if dbl_fit_sse < sgl_fit_sse:
        best_guess = [a_par1, n_par1, ea_par1, a_par2, n_par2, ea_par2]
    else:
        print("nonlinear solver was not better than single exponential!")
        best_guess = init_guess
            
    return best_guess, fit_range


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


def get_valid_temps_rate_constants(temps, rate_constants,
                                   tmin=None, tmax=None):
    """ this subroutine takes in a array of rate constants and
        returns the subset of this array that is positive,
        along with the corresponding Temperature array """

    # Convert temps and rate constants to floats
    temps = [float(temp) for temp in temps]
    rate_constants = [float(rate_constant)
                      if rate_constant != '***' else rate_constant
                      for rate_constant in rate_constants]

    # Set tmin and tmax
    if tmin is None:
        tmin = min(temps)
    if tmax is None:
        tmax = max(temps)
    assert tmin in temps and tmax in temps

    # Grab the temperature, rate constant pairs which correspond to
    # temp > 0, temp within tmin and tmax, rate constant defined (not ***)
    valid_t, valid_k = [], []
    for temp, rate_constant in zip(temps, rate_constants):
        if rate_constant == '***':
            continue
        else:
            if float(rate_constant) > 0.0 and tmin <= temp <= tmax:
                valid_t.append(temp)
                valid_k.append(rate_constant)

    # Convert the lists to numpy arrays
    valid_t = np.array(valid_t, dtype=np.float64)
    valid_k = np.array(valid_k, dtype=np.float64)

    return valid_t, valid_k


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
