"""
Take data dictionaries from mechanisms and combine them under a common index
"""

import itertools
import numpy as np


def thermo(m1_thermo_dct, m2_thermo_dct, temps):
    """ Loop over the the Mech1 thermo entries
    """

    total_thermo_dct = {}

    # Build full thermo dictionary with common index
    # First loop through m1: add common and m1-unique species
    for m1_name, m1_vals in m1_thermo_dct.items():
        m2_vals = m2_thermo_dct.get(m1_name, None)
        total_thermo_dct[m1_name] = {
            'm1': m1_vals,
            'm2': m2_vals
        }

    # Now add the entries where m2 exists, but m1 does not
    uni_m2_names = [name
                    for name in m2_thermo_dct.keys()
                    if name not in m1_thermo_dct]

    for m2_name in uni_m2_names:
        total_thermo_dct[m2_name] = {
            'm1': None,
            'm2': m2_thermo_dct[m2_name]
        }

    return total_thermo_dct


def combine_reaction_rates(m1_rxn_dct, m2_rxn_dct,
                           m1_units, m2_units,
                           m2_thermo_dct,
                           t_ref, temps, pressures):
    """ calculate the reactions rates for two mech files
    """

    total_ktp_dct = {}

    # Get the thermo dictionaries for each mechanism
    m1_ktp_dct = mechparser.reaction.calculate_rate_constants(
        m1_rxn_dct, m1_units, t_ref, temps, pressures=pressures)
    m2_ktp_dct = mechparser.reaction.calculate_rate_constants(
        m2_rxn_dct, m2_units, t_ref, temps, pressures=pressures)

    # Build full rates dictionary with common index
    # First loop through m1: add common and m1-unique species
    for m1_name, m1_ktp in m1_ktp_dct.items():

        # Check what (if/any) combination of m2 matches with m1
        m2_name_match, reverse_rates = _assess_reaction_match(
            m1_name, m2_ktp_dct)

        # Calculate reaction rates, reverse if needed
        if m2_name_match:
            if not reverse_rates:
                m2_ktp = m2_ktp_dct[m1_name]
            else:
                m2_ktp = _reverse_reaction_rates(
                    m2_ktp_dct, m2_thermo_dct, m2_name_match, temps)
        else:
            m2_ktp_dct = None

        # Add entry to overal thermo dictionary
        total_ktp_dct[m1_name] = {
            'm1': m1_ktp,
            'm2': m2_ktp
        }

        # Now add the entries where m2 exists, but m1 does not
        # add the code to do this

    return total_ktp_dct


# Thermo functions
# Functions to build dictionaries
def build_name_dcts(mech1_str, mech2_str):
    """ builds the thermo dictionaries indexed by names
    """

    mech1_thermo_block = mechparser.util.clean_up_whitespace(
        mechparser.mechanism.thermo_block(mech1_str))
    mech1_thermo_dct = mechparser.thermo.dct_name_idx(
        mech1_thermo_block)

    mech2_thermo_block = mechparser.util.clean_up_whitespace(
        mechparser.mechanism.thermo_block(mech2_str))
    mech2_thermo_dct = mechparser.thermo.dct_name_idx(
        mech2_thermo_block)

    return mech1_thermo_dct, mech2_thermo_dct


def build_inchi_dcts(mech1_str, mech2_str,
                     mech1_csv_str, mech2_csv_str):
    """ builds new thermo dictionaries indexed by inchis
    """

    mech1_thermo_block = mechparser.util.clean_up_whitespace(
        mechparser.mechanism.thermo_block(mech1_str))
    mech1_name_dct = mechparser.mechanism.species_name_inchi_dct(
        mech1_csv_str)
    mech1_thermo_dct = mechparser.thermo.dct_inchi_idx(
        mech1_thermo_block, mech1_name_dct)

    mech2_thermo_block = mechparser.util.clean_up_whitespace(
        mechparser.mechanism.thermo_block(mech2_str))
    mech2_name_dct = mechparser.mechanism.species_name_inchi_dct(
        mech2_csv_str)
    mech2_thermo_dct = mechparser.thermo.dct_inchi_idx(
        mech2_thermo_block, mech2_name_dct)

    return mech1_thermo_dct, mech2_thermo_dct


def mech_name_for_species(mech1_csv_str, mech2_csv_str, ich):
    """ build dictionaries to get the name for a given InCHI string
    """

    mech1_inchi_dct = mechparser.mechanism.species_inchi_name_dct(
        mech1_csv_str)
    mech2_inchi_dct = mechparser.mechanism.species_inchi_name_dct(
        mech2_csv_str)

    if ich in mech1_inchi_dct:
        mech1_name = mech1_inchi_dct[ich]
    else:
        mech1_name = 'Not in Mechanism'
    if ich in mech2_inchi_dct:
        mech2_name = mech2_inchi_dct[ich]
    else:
        mech2_name = 'Not in Mechanism'

    return mech1_name, mech2_name

# Rate functions
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




