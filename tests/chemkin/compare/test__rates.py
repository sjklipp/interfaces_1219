"""
compare thermo
"""

import os
import itertools
import chemkin_io


# Get mechanism information
def _read_file(file_name):
    with open(file_name, encoding='utf8', errors='ignore') as file_obj:
        file_str = file_obj.read()
    return file_str


PATH = os.path.dirname(os.path.realpath(__file__))

MECH1_PATH = os.path.join(PATH, '../mechparser/data/syngas')
MECH1_STR = _read_file(os.path.join(
    MECH1_PATH, 'mechanism_green.txt'))
MECH1_CSV_STR = _read_file(os.path.join(
    MECH1_PATH, 'smiles_green.csv'))

MECH2_PATH = os.path.join(PATH, '../mechparser/data/syngas')
MECH2_STR = _read_file(os.path.join(
    MECH2_PATH, 'mechanism_elte_2016_ijck_alt.txt'))
MECH2_CSV_STR = _read_file(os.path.join(
    MECH2_PATH, 'smiles_elte_alt.csv'))

INDEX = 'inchi'

# Temperatures to run
TEMPS = [500.0, 1000.0, 1500.0]
# TEMPS = [300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0,
#          1300.0, 1400.0, 1500.0, 1600.0, 1700.0, 1800.0, 1900.0, 2000.0,
#          2300.0, 2400.0, 2500.0, 2600.0, 2700.0, 2800.0, 2900.0, 3000.0]


# Functions to build dictionaries
def build_name_dcts(mech1_str, mech2_str):
    """ builds the reaction dictionaries indexed by names
    """

    mech1_reaction_block = chemkin_io.mechparser.util.clean_up_whitespace(
        chemkin_io.mechparser.mechanism.reaction_block(mech1_str))
    mech1_reaction_dct = chemkin_io.mechparser.reaction.dct_name_idx(
        mech1_reaction_block)

    mech2_reaction_block = chemkin_io.mechparser.util.clean_up_whitespace(
        chemkin_io.mechparser.mechanism.reaction_block(mech2_str))
    mech2_reaction_dct = chemkin_io.mechparser.reaction.dct_name_idx(
        mech2_reaction_block)

    return mech1_reaction_dct, mech2_reaction_dct


def build_inchi_dcts(mech1_str, mech2_str,
                     mech1_csv_str, mech2_csv_str):
    """ builds new reaction dictionaries indexed by inchis
    """

    mech1_reaction_block = chemkin_io.mechparser.util.clean_up_whitespace(
        chemkin_io.mechparser.mechanism.reaction_block(mech1_str))
    mech1_name_dct = chemkin_io.mechparser.mechanism.species_name_inchi_dct(
        mech1_csv_str)
    mech1_reaction_dct = chemkin_io.mechparser.reaction.dct_inchi_idx(
        mech1_reaction_block, mech1_name_dct)

    mech2_reaction_block = chemkin_io.mechparser.util.clean_up_whitespace(
        chemkin_io.mechparser.mechanism.reaction_block(mech2_str))
    mech2_name_dct = chemkin_io.mechparser.mechanism.species_name_inchi_dct(
        mech2_csv_str)
    mech2_reaction_dct = chemkin_io.mechparser.reaction.dct_inchi_idx(
        mech2_reaction_block, mech2_name_dct)

    return mech1_reaction_dct, mech2_reaction_dct


def get_reaction_units(mech1_str, mech2_str):
    """ get reaction units
    """

    mech1_reaction_block = chemkin_io.mechparser.util.clean_up_whitespace(
        chemkin_io.mechparser.mechanism.reaction_block(mech1_str))
    mech1_units = chemkin_io.mechparser.reaction.units(
        mech1_reaction_block)

    mech2_reaction_block = chemkin_io.mechparser.util.clean_up_whitespace(
        chemkin_io.mechparser.mechanism.reaction_block(mech2_str))
    mech2_units = chemkin_io.mechparser.reaction.units(
        mech2_reaction_block)

    return mech1_units, mech2_units


def calculate_reaction_rates(mech1_rxn_dct, mech2_rxn_dct,
                             mech1_units, mech2_units,
                             t_ref, temps, pressures):
    """ calculate the reactions rates for two mech files
    """

    reaction_dct = {}

    for m1_name, m1_rxn_dstr in mech1_rxn_dct.items():
        m1_ktp_dct = chemkin_io.mechparser.reaction.calculate_rate_constants(
            m1_rxn_dstr, t_ref, mech1_units, temps, pressures=pressures)

        # Check what combination of m2 matches with m1
        m2_name_match, flip_rxn = _assess_reaction_match(
             m1_name, mech2_rxn_dct)

        # Calculate reaction rates
        if m2_name_match:
            m2_rxn_dstr = mech2_rxn_dct[m2_name_match]
            m2_ktp_dct = chemkin_io.mechparser.reaction.calculate_rate_constants(
                m2_rxn_dstr, t_ref, mech2_units, temps, pressures=pressures)
            if flip_rxn:
                m2_ktp_dct = _reverse_reaction_rates(
                    m2_ktp_dct, m2_thermo_dct, m2_name_match,)
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
    m1_rct_comb = list(itertools.combinations(m1_rcts))
    m1_prd_comb = list(itertools.combinations(m1_prds))

    for rxn in m2_dct:
        [m2_rcts, m2_prds] = rxn
        if m2_rcts in m1_rct_comb and m2_prds in m1_prd_comb:
            flip_rxn = False
            m2_key = rxn
        elif m2_rcts in m1_prd_comb and m2_prds in m1_rct_comb:
            m2_key = rxn
            flip_rxn = True
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
        for rate, k_equil in rates, k_equils:
            rev_rates.append(rate / k_equil)
        rev_ktp_dct[pressure] = rev_rates

    return rev_ktp_dct


def _calculate_equilibrium_constant(thermo_dct, rct_idxs, prd_idxs, temp):
    """ use the thermo parameters to obtain the equilibrium
        constant
    """

    rct_gibbs = 0.0
    for rct in rct_idxs:
        rct_gibbs += chemkin_io.mechparser.thermo.calculate_gibbs(
            thermo_dct[rct], temp)
    prd_gibbs = 0.0
    for prd in prd_idxs:
        prd_gibbs += chemkin_io.mechparser.thermo.calculate_gibbs(
            thermo_dct[prd], temp)

    rxn_gibbs = prd_gibbs - rct_gibbs

    k_equil = np.exp(-rxn_gibbs / (RC * temp))

    return k_equil



if __name__ == '__main__':

    if INDEX == 'name':
        M1_RXN_DCT, M2_RXN_DCT = build_name_dcts(
            MECH1_STR, MECH2_STR)
    elif INDEX == 'inchi':
        M1_RXN_DCT, M2_RXN_DCT = build_inchi_dcts(
            MECH1_STR, MECH2_STR, MECH1_CSV_STR, MECH2_CSV_STR)

    print('\nMech1 thermo dct')
    for key, val in M1_RXN_DCT.items():
        print(key)
        print(val)
    print('\n\nMech2 thermo dct')
    for key, val in M2_RXN_DCT.items():
        print(key)
        print(val)

    M1_UNITS = chemkin_io.mechparser.reaction.units(
        MECH1_STR)
    M2_UNITS = chemkin_io.mechparser.reaction.units(
        MECH2_STR)

    # THM_VALS_DCT = calculate_mech_thermo(
    #     M1_THM_DCT, M2_THM_DCT, TEMPS)

    # print('\n\n\nCombined thermo vals')
    # for idx in THM_VALS_DCT:
    #     if INDEX == 'name':
    #         m1_name = idx
    #         m2_name = idx
    #     elif INDEX == 'inchi':
    #         print('\n\nInCHI: ', idx)
    #         m1_name, m2_name = get_mech_name_for_species(
    #             MECH1_CSV_STR, MECH2_CSV_STR, idx)
    #     print('M1 Name: ', m1_name)
    #     print('M2 Name: ', m2_name)
    #     m1_vals = THM_VALS_DCT[idx]['m1']
    #     m2_vals = THM_VALS_DCT[idx]['m2']
    #     print('\nM1 Enthalpy', m1_vals[0])
    #     print('M2 Enthalpy', m2_vals[0])
    #     print('M1 Entropy', m1_vals[1])
    #     print('M2 Entropy', m2_vals[1])
    #     print('M1 Gibbs', m1_vals[2])
    #     print('M2 Gibbs', m2_vals[2])
    #     print('M1 Heat Capacity', m1_vals[3])
    #     print('M2 Heat Capacity', m2_vals[3])
