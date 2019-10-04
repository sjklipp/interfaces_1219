""" functions operating on the thermo block string
"""


import numpy as np


RC = 1.98720425864083e-3  # in kcal/mol.K


# functions which calculate quantiies using data from the thermo section #
def mechanism(nasa_dct, temps):
    """ Loop over the the Mech1 thermo entries
    """

    mech_thermo_dct = {}
    for name, thermo_dstr in nasa_dct.items():
        h_t, cp_t, s_t, g_t, = [], [], [], []
        for temp in temps:
            h_t.append(enthalpy(thermo_dstr, temp))
            cp_t.append(heat_capacity(thermo_dstr, temp))
            s_t.append(entropy(thermo_dstr, temp))
            g_t.append(gibbs(thermo_dstr, temp))

        mech_thermo_dct[name] = [h_t, cp_t, s_t, g_t]

    return mech_thermo_dct


def enthalpy(thm_dstr, temp):
    """ Calculate the Enthalpy [H(T)] of a species using the
        coefficients of its NASA polynomial
    """

    cfts = _coefficients_for_specific_temperature(thm_dstr, temp)

    if cfts is not None:
        h_t = (
            cfts[0] +
            ((cfts[1] * temp) / 2.0) +
            ((cfts[2] * temp**2) / 3.0) +
            ((cfts[3] * temp**3) / 4.0) +
            ((cfts[4] * temp**4) / 5.0) +
            (cfts[5] / temp)
        )
        h_t *= (RC * temp)
    else:
        h_t = 0.0

    return h_t


def heat_capacity(thm_dstr, temp):
    """ Calculate the Heat Capacity [Cp(T)] of a species using the
        coefficients of its NASA polynomial
    """
    cfts = _coefficients_for_specific_temperature(thm_dstr, temp)

    if cfts is not None:
        cp_t = (
            cfts[0] +
            (cfts[1] * temp) +
            (cfts[2] * temp**2) +
            (cfts[3] * temp**3) +
            (cfts[4] * temp**4)
        )
        cp_t *= RC
    else:
        cp_t = 0.0

    return cp_t


def entropy(thm_dstr, temp):
    """ Calculate the Entropy [S(T)] of a species using the
        coefficients of its NASA polynomial
    """
    cfts = _coefficients_for_specific_temperature(thm_dstr, temp)

    if cfts is not None:
        s_t = (
            (cfts[0] * np.log(temp)) +
            (cfts[1] * temp) +
            ((cfts[2] * temp**2) / 2.0) +
            ((cfts[3] * temp**3) / 3.0) +
            ((cfts[4] * temp**4) / 4.0) +
            (cfts[6])
        )
        s_t *= RC
    else:
        s_t = 0.0

    return s_t


def gibbs(thm_dstr, temp):
    """ Calculate the Gibbs Free Energy [H(T)] of a species using the
        coefficients of its NASA polynomial
    """

    h_t = enthalpy(thm_dstr, temp)
    s_t = entropy(thm_dstr, temp)
    if enthalpy is not None and entropy is not None:
        g_t = h_t - (s_t * temp)
    else:
        g_t = None

    return g_t


def _coefficients_for_specific_temperature(thm_dstr, temp):
    """ return the set of coefficients of the polynomial (low or high)
        that should be used for a given temperature
    """

    temps = temperatures(thm_dstr)
    if temps[0] < temp < temps[1]:
        cfts = low_coefficients(thm_dstr)
    elif temps[1] < temp < temps[2]:
        cfts = high_coefficients(thm_dstr)
    else:
        cfts = None
        # raise ValueError('Temperature outside range of NASA polynomial')

    return cfts
