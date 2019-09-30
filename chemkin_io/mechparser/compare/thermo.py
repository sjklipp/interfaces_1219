"""
compare thermo
"""

from chemkin_io import mechparser


RC = 1.98720425864083e-3  # Gas Constant in kcal/mol.K


def calculate_mech_thermo(m1_thermo_dct, m2_thermo_dct, temps):
    """ Loop over the the Mech1 thermo entries
    """

    thermo_dct = {}

    # Calculate all thermo quanties with mech1 keys existing
    for m1_name, m1_thermo_dstr in m1_thermo_dct.items():

        # Calculate the thermo values for mech1
        m1_enthalpy, m1_entropy, m1_gibbs, m1_cp = [], [], [], []
        for temp in temps:
            m1_enthalpy.append(
                mechparser.thermo.calculate_enthalpy(
                    m1_thermo_dstr, temp))
            m1_entropy.append(
                mechparser.thermo.calculate_entropy(
                    m1_thermo_dstr, temp))
            m1_gibbs.append(
                mechparser.thermo.calculate_gibbs(
                    m1_thermo_dstr, temp))
            m1_cp.append(
                mechparser.thermo.calculate_heat_capacity(
                    m1_thermo_dstr, temp))

        # Calculate the thermo values for mech2
        m2_enthalpy, m2_entropy, m2_gibbs, m2_cp = [], [], [], []
        if m1_name in m2_thermo_dct:
            m2_thermo_dstr = m2_thermo_dct[m1_name]
            for temp in temps:
                m2_enthalpy.append(
                    mechparser.thermo.calculate_enthalpy(
                        m2_thermo_dstr, temp))
                m2_entropy.append(
                    mechparser.thermo.calculate_entropy(
                        m2_thermo_dstr, temp))
                m2_gibbs.append(
                    mechparser.thermo.calculate_gibbs(
                        m2_thermo_dstr, temp))
                m2_cp.append(
                    mechparser.thermo.calculate_heat_capacity(
                        m2_thermo_dstr, temp))
        else:
            m2_enthalpy, m2_entropy, m2_gibbs, m2_cp = None, None, None, None

        # Add entry to overal thermo dictionary
        thermo_dct[m1_name] = {
            'm1': [m1_enthalpy, m1_entropy, m1_gibbs, m1_cp],
            'm2': [m2_enthalpy, m2_entropy, m2_gibbs, m2_cp],
        }

    # Now add the entries where m2 exists, but m1 does not
    uni_m2_names = [name
                    for name in m2_thermo_dct.keys()
                    if name not in m1_thermo_dct]
    for name in uni_m2_names:

        # Set values for mech1
        m1_enthalpy, m1_entropy, m1_gibbs, m1_cp = None, None, None, None

        # Calculate values for mech2
        m2_enthalpy, m2_entropy, m2_gibbs, m2_cp = [], [], [], []
        for temp in temps:
            m2_enthalpy.append(
                mechparser.thermo.calculate_enthalpy(
                    m2_thermo_dct[name], temp))
            m2_entropy.append(
                mechparser.thermo.calculate_entropy(
                    m2_thermo_dct[name], temp))
            m2_gibbs.append(
                mechparser.thermo.calculate_gibbs(
                    m2_thermo_dct[name], temp))
            m2_cp.append(
                mechparser.thermo.calculate_heat_capacity(
                    m2_thermo_dct[name], temp))

        # Add entry to overal thermo dictionary
        thermo_dct[name] = {
            'm1': [m1_enthalpy, m1_entropy, m1_gibbs, m1_cp],
            'm2': [m2_enthalpy, m2_entropy, m2_gibbs, m2_cp],
        }

    return thermo_dct


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
