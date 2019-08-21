""" fit rate constants to Arrhenius expressions
"""

import numpy as np

# put in QCEngine gas constant
R = 8.314


def get_valid_Tk(data, k):
    """ this subroutine takes in a array of rate constants and
        returns the subset of this array that is positive,
        along with the corresponding Temperature array """

    # start by using only the temperatures at which the rate constant is well defined
    local_t = []
    local_k = []
    for (t, T) in enumerate(data.temperature):
        if (k[t] > 0.0) and (T >= data.Tmin) and (T <= data.Tmax):
            local_t.append(T)
            local_k.append(k[t])

    local_t = np.array(local_t, dtype=np.float64)
    local_k = np.array(local_k, dtype=np.float64)

    return local_t, local_k


def fit_single_arrhenius(t, k):
    """ this subroutine takes in a vector of rate constants and
        returns the Arrhenius parameters, as well as
        the T-range over which they were fit"""

    # obtain temperatures at which the rate constant is well defined
    local_t, local_k = get_valid_Tk(t, k)

    # there are many cases to consider, depending upon the number of valid k's
    # no k is positive, so return all zeros
    if len(local_k) == 0:
        A_fit, n_fit, Ea_fit = 0.0, 0.0, 0.0
        fit_range = [0, 0]

    # if num(k) > 0 is 1: set A = k
    elif len(local_k) == 1:
        A_fit, n_fit, Ea_fit = local_k, 0.0, 0.0
        fit_range = [t.index(min(local_t)), t.index(max(local_t))]

    # if num(k) > 0 is 2,3: fit A and Ea
    elif (len(local_k) == 2) or (len(local_k) == 3):
        X = np.array([np.ones(len(local_t)), -1.0 / R / local_t],
                      dtype=np.float64)
        X = X.transpose()
        theta = np.linalg.lstsq(X, np.log(local_k))[0]
        A_fit = np.exp(theta[0]), 0.0, Ea_fit = theta[1]
        fit_range = [t.index(min(local_t)), t.index(max(local_t))]

    # if num(k) > 0 is more than 3: fit A, n, and Ea
    elif len(local_k) > 3:
        X = np.array([np.ones(len(local_t)), np.log(local_t/data.T0), -1.0 / R / local_t],
                     dtype=np.float64)
        X = X.transpose()
        theta = np.linalg.lstsq(X, np.log(local_k))[0]
        A_fit, n_fit, Ea_fit = np.exp(theta[0]), theta[1], theta[2]
        fit_range = [t.index(min(local_t)), t.index(max(local_t))]

    return A_fit, n_fit, Ea_fit, fit_range
