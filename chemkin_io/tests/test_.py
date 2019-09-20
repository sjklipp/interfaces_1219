""" test chemkion_io.mech
"""
from __future__ import unicode_literals
from builtins import open
import os
import pandas
import automol
import chemkin_io


def _read_file(file_name):
    with open(file_name, encoding='utf8', errors='ignore') as file_obj:
        file_str = file_obj.read()
    return file_str


PATH = os.path.dirname(os.path.realpath(__file__))
NATGAS_PATH = os.path.join(PATH, 'data/natgas')
HEPTANE_PATH = os.path.join(PATH, 'data/heptane')
NATGAS_MECH_STR = _read_file(os.path.join(NATGAS_PATH, 'mechanism.txt'))

HEPTANE_MECH_STR = _read_file(os.path.join(HEPTANE_PATH, 'mechanism.txt'))
HEPTANE_TAB = pandas.read_csv(os.path.join(HEPTANE_PATH, 'smiles.csv'))
HEPTANE_TAB['inchi'] = list(map(automol.smiles.inchi, HEPTANE_TAB['smiles']))


def test__species_block():
    """ test chemkin_io.species_block
    """
    mech_str = NATGAS_MECH_STR
    block_str = chemkin_io.species_block(mech_str)
    assert len(block_str.splitlines()) == 131


def test__reaction_block():
    """ test chemkin_io.reaction_block
    """
    mech_str = NATGAS_MECH_STR
    block_str = chemkin_io.reaction_block(mech_str)
    assert len(block_str.splitlines()) == 1834


def test__thermo_block():
    """ test chemkin_io.thermo_block
    """
    mech_str = NATGAS_MECH_STR
    block_str = chemkin_io.thermo_block(mech_str)
    for x in block_str:
        print(x)
    assert len(block_str.splitlines()) == 522


def test__species__names():
    """ test mech.species.names
    """
    mech_str = NATGAS_MECH_STR
    block_str = chemkin_io.species_block(mech_str)
    spc_names = chemkin_io.species.names(block_str)
    assert len(spc_names) == 130


def test__reaction__data_strings():
    """ test mech.reaction.data_strings
    """
    mech_str = NATGAS_MECH_STR
    block_str = chemkin_io.reaction_block(mech_str)
    rxn_strs = chemkin_io.reaction.data_strings(block_str)
    assert len(rxn_strs) == 1678

    rct_names_lst = list(
        map(chemkin_io.reaction.DataString.reactant_names, rxn_strs))
    prd_names_lst = list(
        map(chemkin_io.reaction.DataString.product_names, rxn_strs))
    coeffs_lst = list(
        map(chemkin_io.reaction.DataString.high_p_coefficients, rxn_strs))

    # make sure we don't have any None's
    assert all(rct_names_lst)
    assert all(prd_names_lst)
    assert all(coeffs_lst)

    # generating inchis for testing elsewhere
    mech_str = HEPTANE_MECH_STR
    block_str = chemkin_io.reaction_block(mech_str)
    rxn_strs = chemkin_io.reaction.data_strings(block_str)
    rct_names_lst = list(
        map(chemkin_io.reaction.DataString.reactant_names, rxn_strs))
    prd_names_lst = list(
        map(chemkin_io.reaction.DataString.product_names, rxn_strs))

    ich_dct = dict(zip(HEPTANE_TAB['name'], HEPTANE_TAB['inchi']))
    mult_dct = dict(zip(HEPTANE_TAB['name'], HEPTANE_TAB['mult']))
    ich_dct['OHV'] = None
    mult_dct['OHV'] = None
    ich_dct['CHV'] = None
    mult_dct['CHV'] = None

    rct_ichs_lst = list(
        list(map(ich_dct.__getitem__, names)) for names in rct_names_lst)
    rct_mults_lst = list(
        list(map(mult_dct.__getitem__, names)) for names in rct_names_lst)
    prd_ichs_lst = list(
        list(map(ich_dct.__getitem__, names)) for names in prd_names_lst)
    prd_mults_lst = list(
        list(map(mult_dct.__getitem__, names)) for names in prd_names_lst)

    def _accept(ichs):
        return (
            all(ichs) and len(ichs) <= 2 and not
            any(ich in automol.convert.inchi.HARDCODED_INCHI_DCT
                for ich in ichs)
        )

    r_ich_lst = []
    p_ich_lst = []
    r_mlt_lst = []
    p_mlt_lst = []
    for rct_ichs, prd_ichs, rct_mults, prd_mults in zip(
            rct_ichs_lst, prd_ichs_lst, rct_mults_lst, prd_mults_lst):
        if _accept(rct_ichs) and _accept(prd_ichs):
            r_idxs = automol.inchi.argsort(rct_ichs)
            p_idxs = automol.inchi.argsort(prd_ichs)
            r_ichs = [rct_ichs[idx] for idx in r_idxs]
            r_mults = [rct_mults[idx] for idx in r_idxs]
            p_ichs = [prd_ichs[idx] for idx in p_idxs]
            p_mults = [prd_mults[idx] for idx in p_idxs]

            r_ich = automol.inchi.standard_form(automol.inchi.join(r_ichs))
            p_ich = automol.inchi.standard_form(automol.inchi.join(p_ichs))

            assert (list(map(automol.inchi.standard_form, r_ichs)) ==
                    list(map(automol.inchi.standard_form,
                             automol.inchi.split(r_ich))))
            assert (list(map(automol.inchi.standard_form, p_ichs)) ==
                    list(map(automol.inchi.standard_form,
                             automol.inchi.split(p_ich))))

            r_mlt = '_'.join(map(str, r_mults))
            p_mlt = '_'.join(map(str, p_mults))

            print(r_ich, p_ich)
            print(r_mlt, p_mlt)

            r_ich_lst.append(r_ich)
            p_ich_lst.append(p_ich)
            r_mlt_lst.append(r_mlt)
            p_mlt_lst.append(p_mlt)

    df = pandas.DataFrame.from_dict({
        'prod_inchi': p_ich_lst,
        'prod_mults': p_mlt_lst,
        'reac_inchi': r_ich_lst,
        'reac_mults': r_mlt_lst,
    })

    df.to_csv('x.csv', index=False)


def test__thermo__data_strings():
    """ test mech.thermo.data_strings
    """
    mech_str = NATGAS_MECH_STR
    block_str = chemkin_io.thermo_block(mech_str)
    thm_strs = chemkin_io.thermo.data_strings(block_str)
    for x in thm_strs:
        print(thm_strs)
    assert len(thm_strs) == 130


if __name__ == '__main__':
    # test__species_block()
    # test__reaction_block()
    # test__thermo_block()
    # test__species__names()
    # test__reaction__data_strings()
    test__thermo__data_strings()
