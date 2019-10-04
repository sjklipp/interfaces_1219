""" functions operating on the reactions block string
"""


import numpy as np
from qcelemental import constants as qcc
import ratefit


# Constants and Conversion factors
# NAVO = qcc.constants.avogadro_constant
NAVO = 6.0221409e+23
CAL2KCAL = qcc.conversion_factor('cal/mol', 'kcal/mol')
J2KCAL = qcc.conversion_factor('J/mol', 'kcal/mol')
KJ2KCAL = qcc.conversion_factor('kJ/mol', 'kcal/mol')
KEL2KCAL = qcc.conversion_factor('kelvin', 'kcal/mol')


def mechanism(rxn_dct, units, t_ref, temps, pressures):
    """ calculate the reactions rates for a whole block via a dict
    """

    ktp_dct = {}
    for name, rxn_dstr in rxn_dct.items():
        ktp_dct[name] = reaction(
            rxn_dstr, t_ref, units, temps, pressures=pressures)

    return ktp_dct


def reaction(rxn_str, t_ref, rxn_units, temps, pressures=None):
    """ calculate the rate constant using the rxn_string
    """
    rate_constants = {}

    # Accepts a params dictionary
    # Read the parameters from the reactions string
    highp_params = high_p_parameters(rxn_str)
    lowp_params = low_p_parameters(rxn_str)
    troe_params = troe_parameters(rxn_str)
    chebyshev_params = chebyshev_parameters(rxn_str)
    plog_params = plog_parameters(rxn_str)
    # print('\nlocated params')
    # print(highp_params)
    # print(lowp_params)
    # print(troe_params)
    # print(chebyshev_params)
    # print(plog_params)

    # Calculate high_pressure rates
    highp_params = _update_params_units(highp_params, rxn_units)
    highp_ks = ratefit.fxns.arrhenius(highp_params, t_ref, temps)
    rate_constants['high'] = highp_ks

    # Calculate pressure-dependent rate constants based on discovered params
    # Either (1) Plog, (2) Chebyshev, (3) Lindemann, or (4) Troe
    # Update units if necessary
    if any(params is not None
           for params in (plog_params, chebyshev_params, lowp_params)):
        assert pressures is not None

    pdep_dct = {}
    if plog_params is not None:
        updated_plog_params = []
        for params in plog_params:
            updated_plog_params.append(
                [params[0]] + _update_params_units(params[1:], rxn_units))
        pdep_dct = _plog(updated_plog_params, pressures, temps, t_ref)

    elif chebyshev_params is not None:
        pdep_dct = _chebyshev(chebyshev_params, pressures, temps)

    elif lowp_params is not None:
        lowp_params = _update_params_units(lowp_params, rxn_units)
        lowp_ks = ratefit.fxns.arrhenius(lowp_params, t_ref, temps)
        if troe_params is not None:
            pdep_dct = _troe(troe_params, highp_ks, lowp_ks, pressures, temps)
        else:
            pdep_dct = ratefit.fxns.lindemann(
                highp_ks, lowp_ks, pressures, temps)

    if pdep_dct:
        for key, val in pdep_dct.items():
            rate_constants[key] = val

    return rate_constants


def _update_params_units(params, rxn_units):
    """ change the units if necessary
        only needed for highp, lowp, and plog
    """
    # Determine converstion factor for Ea Units
    ea_units = rxn_units[0]
    if ea_units == 'cal/mole':
        ea_conv_factor = CAL2KCAL
    elif ea_units == 'joules/mole':
        ea_conv_factor = J2KCAL
    elif ea_units == 'kjoules/mole':
        ea_conv_factor = KJ2KCAL
    elif ea_units == 'kelvin':
        ea_conv_factor = KEL2KCAL
    else:
        ea_conv_factor = 1.0

    # Determine converstion factor for A Units
    if rxn_units[1] == 'molecules':
        a_conv_factor = NAVO
    else:
        a_conv_factor = 1.0

    # update units of params
    if params is not None:
        params[2] *= ea_conv_factor
        # if len(params) == 6:
        #     params[5] *= ea_conv_factor

        params[0] *= a_conv_factor
        # if len(params) == 6:
        #     params[3] *= a_conv_factor

    return params


def _plog(plog_params, pressures, temps, t_ref):
    """ calc plog
    """
    plog_dct = {}
    for params in plog_params:
        plog_dct[params[0]] = params[1:]
    pdep_dct = ratefit.fxns.plog(plog_dct, t_ref, pressures, temps)
    return pdep_dct


def _chebyshev(chebyshev_params, pressures, temps):
    """ calc chebyshev
    """
    tmin = chebyshev_params[0][0]
    tmax = chebyshev_params[0][1]
    pmin = chebyshev_params[1][0]
    pmax = chebyshev_params[1][1]
    alpha = np.array(chebyshev_params[3])
    pdep_dct = ratefit.fxns.chebyshev(
        alpha, tmin, tmax, pmin, pmax, pressures, temps)
    return pdep_dct


def _troe(troe_params, highp_ks, lowp_ks, pressures, temps):
    """ calc troe
    """
    if len(troe_params) == 3:
        ts2 = None
    elif len(troe_params) == 4:
        ts2 = troe_params[3]
    pdep_dct = ratefit.fxns.troe(
        highp_ks, lowp_ks, pressures, temps,
        troe_params[0], troe_params[1], troe_params[2], ts2=ts2)
    return pdep_dct
