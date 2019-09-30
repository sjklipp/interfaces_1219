"""
compare thermo
"""

import itertools
import numpy as np
from chemkin_io import mechparser
# from chemkin_io.mechparser.compare.therm as comptherm


RC = 1.98720425864083e-3  # Gas Constant in kcal/mol.K


def calculate_reaction_rates(mech1_rxn_dct, mech2_rxn_dct,
                             mech1_units, mech2_units,
                             mech2_thermo_dct,
                             t_ref, temps, pressures):
    """ calculate the reactions rates for two mech files
    """

    reaction_dct = {}

    for m1_name, m1_rxn_dstr in mech1_rxn_dct.items():
        m1_ktp_dct = mechparser.reaction.calculate_rate_constants(
            m1_rxn_dstr, t_ref, mech1_units, temps, pressures=pressures)

        # Check what combination of m2 matches with m1
        m2_name_match, flip_rxn = _assess_reaction_match(
            m1_name, mech2_rxn_dct)

        # Calculate reaction rates
        if m2_name_match:
            m2_rxn_dstr = mech2_rxn_dct[m2_name_match]
            m2_ktp_dct = mechparser.reaction.calculate_rate_constants(
                m2_rxn_dstr, t_ref, mech2_units, temps, pressures=pressures)
            if flip_rxn:
                m2_ktp_dct = _reverse_reaction_rates(
                    m2_ktp_dct, mech2_thermo_dct, m2_name_match, temps)
        else:
            m2_ktp_dct = None

        # Add entry to overal thermo dictionary
        reaction_dct[m1_name] = {
            'm1': m1_ktp_dct,
            'm2': m2_ktp_dct
        }

    return reaction_dct


def _assess_reaction_match(m1_key, m2_dct):
    """ assess whether the reaction should be flipped
    """

    [m1_rcts, m1_prds] = m1_key
    m1_rct_comb = list(itertools.combinations(m1_rcts, len(m1_rcts)))
    m1_prd_comb = list(itertools.combinations(m1_prds, len(m1_prds)))

    for rxn in m2_dct:
        [m2_rcts, m2_prds] = rxn
        if m2_rcts in m1_rct_comb and m2_prds in m1_prd_comb:
            flip_rxn = False
            m2_key = rxn
            break
        elif m2_rcts in m1_prd_comb and m2_prds in m1_rct_comb:
            m2_key = rxn
            flip_rxn = True
            break
        else:
            m2_key = ()
            flip_rxn = None

    ret = m2_key, flip_rxn

    return ret


def _reverse_reaction_rates(ktp_dct, thermo_dct, rxn, temps):
    """ use the equilibrium constant to reverse the reaction rates
    """

    [rct_idxs, prd_idxs] = rxn
    k_equils = _calculate_equilibrium_constant(
        thermo_dct, rct_idxs, prd_idxs, temps)

    rev_ktp_dct = {}
    for pressure, rates in ktp_dct.items():
        rev_rates = []
        for rate, k_equil in zip(rates, k_equils):
            rev_rates.append(rate / k_equil)
        rev_ktp_dct[pressure] = rev_rates

    return rev_ktp_dct


def _calculate_equilibrium_constant(thermo_dct, rct_idxs, prd_idxs, temps):
    """ use the thermo parameters to obtain the equilibrium
        constant
    """

    k_equils = []
    for temp in temps:
        rct_gibbs = 0.0
        for rct in rct_idxs:
            rct_gibbs += mechparser.thermo.calculate_gibbs(
                thermo_dct[rct], temp)
        prd_gibbs = 0.0
        for prd in prd_idxs:
            prd_gibbs += mechparser.thermo.calculate_gibbs(
                thermo_dct[prd], temp)

        rxn_gibbs = prd_gibbs - rct_gibbs

        k_equils.append(
            np.exp(-rxn_gibbs / (RC * temp)))

    return k_equils


# Functions to build dictionaries
def build_name_dcts(mech1_str, mech2_str):
    """ builds the reaction dictionaries indexed by names
    """

    mech1_reaction_block = mechparser.util.clean_up_whitespace(
        mechparser.mechanism.reaction_block(mech1_str))
    mech1_reaction_dct = mechparser.reaction.dct_name_idx(
        mech1_reaction_block)

    mech2_reaction_block = mechparser.util.clean_up_whitespace(
        mechparser.mechanism.reaction_block(mech2_str))
    mech2_reaction_dct = mechparser.reaction.dct_name_idx(
        mech2_reaction_block)

    return mech1_reaction_dct, mech2_reaction_dct


def build_inchi_dcts(mech1_str, mech2_str,
                     mech1_csv_str, mech2_csv_str):
    """ builds new reaction dictionaries indexed by inchis
    """

    mech1_reaction_block = mechparser.util.clean_up_whitespace(
        mechparser.mechanism.reaction_block(mech1_str))
    mech1_name_dct = mechparser.mechanism.species_name_inchi_dct(
        mech1_csv_str)
    mech1_reaction_dct = mechparser.reaction.dct_inchi_idx(
        mech1_reaction_block, mech1_name_dct)

    mech2_reaction_block = mechparser.util.clean_up_whitespace(
        mechparser.mechanism.reaction_block(mech2_str))
    mech2_name_dct = mechparser.mechanism.species_name_inchi_dct(
        mech2_csv_str)
    mech2_reaction_dct = mechparser.reaction.dct_inchi_idx(
        mech2_reaction_block, mech2_name_dct)

    return mech1_reaction_dct, mech2_reaction_dct
# def get_mech_names_for_species(mech1_csv_str, mech2_csv_str, ich):
#     """ build dictionaries to get the name for a given InCHI string
#     """
#     mech1_inchi_dct = mechparser.mechanism.species_inchi_name_dct(
#         mech1_csv_str)
#     mech2_inchi_dct = mechparser.mechanism.species_inchi_name_dct(
#         mech2_csv_str)
#
#
