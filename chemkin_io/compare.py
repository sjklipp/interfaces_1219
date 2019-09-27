"""
compare two mechanism files
"""

import ratefit


# move this function to mechparser/reactions.py
def calculate_rate_constants(rxn_str, t_ref, temps, pressures=None):
    """ calculate the rate constant using the rxn_string
    """

    rate_constants = {}

    # Read the coefficients from the reactions string
    highp_cfts = highp_coefficients(rxn_str) 
    lowp_cfts = lowp_coefficients(rxn_str) 
    troe_cfts = troe_coefficients(rxn_str) 
    chebyshev_cfts = chebyshev_coefficients(rxn_str)
    plog_cfts = plog_coefficients(rxn_str) 

    # Calculate the high_pressure rates
    if len(highp_cfts) == 6:
        highp_ks = ratefit.fxns.single_arrhenius(
            highp_cfts[0], highp_cfts[1], highp_cfts[2]
            t_ref, temp)
    else:
        highp_ks = ratefit.fxns.double_arrhenius(
            highp_cfts[0], highp_cfts[1], highp_cfts[2]
            highp_cfts[3], highp_cfts[4], highp_cfts[5]
            t_ref, temp)

    rate_constants['high'] = highp_ks

    # Calculate the pressure dependent rates using either Lindemann or Troe fits
    if lowp_coefficients:
        assert pressures is not None
        if len(lowp_cfts) == 6:
            lowp_ks = ratefit.fxns.single_arrhenius(
                lowp_cfts[0], lowp_cfts[1], lowp_cfts[2]
                t_ref, temp)
        else:
            lowp_ks = ratefit.fxns.double_arrhenius(
                lowp_cfts[0], lowp_cfts[1], lowp_cfts[2]
                lowp_cfts[3], lowp_cfts[4], lowp_cfts[5]
                t_ref, temp)
        for pressure in pressures:
            if not troe_coefficients 
                pdep_ks = lindemann(
                    highp_ks, lowp_ks, pressure, temp)
            else:
                pdep_ks = troe(
                    highp_ks, lowp_ks, pressure, temp, 
                    troe_cfts[0], troe_cfts[1], troe_cfts[2], t2=..)

            rate_constants[pressure] = pdep_ks
    
    # Calculate the pressure dependent rates using each plog expression
    for key, val in plog_cfts.items():

    return rate_constants


def calculate_equilibrium_constant(rcts):
    """ use the thermo parameters to obtain the equilibrium
        constant
    """
    if _assess_reverse_reaction_rates():

    return equilK


def _assess_reverse_reactio_ratesn():
    """ assess whether the reaction should be flipped
    """
    return

