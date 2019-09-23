"""
  Tests the varecof_io.writer functions
"""

import varecof_io.writer

# For test__corr_potentials writer
NPOT = 5
BND_IDX = [1, 3]
BND_SYM = ['C', 'O']
RVALS = [1.5958, 1.6958, 1.7958, 1.8958, 1.9958,
         2.0958, 2.1958, 2.2958, 2.3958, 2.4958]

POTENTIALS = [
    [0.052, 0.175, 0.430, 0.724, 0.996,
     1.199, 1.308, 1.317, 1.243, 1.113],
    [-0.722, -0.517, -0.372, -0.277, -0.218,
     -0.181, -0.153, -0.126, -0.096, -0.064],
    [0.224, 0.329, 0.556, 0.823, 1.071,
     1.255, 1.348, 1.346, 1.263, 1.127]
]
SPECIES_NAME = 'mol'
POT_LABELS = ['basis+relaxed', 'basis', 'relaxed']
FORTRAN_COMPILER = 'gfortran'
SPECIES_CORR_POTENTIALS = ['mol']


# For test__compile_correction_potential
MAKE_PATH = '.'


def test__corr_potentials_writer():
    """ tests varecof_io.writer.corr_potentials.species
              varecof_io.writer.corr_potentials.dummy
              varecof_io.writer.corr_potentials.auxiliary
              varecof_io.writer.corr_potentials.makefile
    """

    # Write the species_corr.f string
    species_corr_str = varecof_io.writer.corr_potentials.species(
        RVALS, POTENTIALS, BND_IDX, BND_SYM,
        species_name=SPECIES_NAME, pot_labels=POT_LABELS)
    # print('\n\nspecies_corr.f:')
    # print(species_corr_str)
    with open('mol_corr.f', 'w') as mol_corr_file:
        mol_corr_file.write(species_corr_str)

    # Write the dummy_corr.f string
    dummy_corr_str = varecof_io.writer.corr_potentials.dummy()
    # print('\n\ndummy_corr.f:')
    # print(dummy_corr_str)
    with open('dummy_corr.f', 'w') as dummy_corr_file:
        dummy_corr_file.write(dummy_corr_str)

    # Write the pot_aux.f string
    pot_aux_str = varecof_io.writer.corr_potentials.auxiliary()
    # print('\n\npot_aux.f:')
    # print(pot_aux_str)
    with open('pot_aux.f', 'w') as pot_aux_file:
        pot_aux_file.write(pot_aux_str)

    # Write the makefile string
    makefile_str = varecof_io.writer.corr_potentials.makefile(
        FORTRAN_COMPILER,
        species_corr_potentials=SPECIES_CORR_POTENTIALS)
    # print('\n\nmakefile:')
    # print(makefile_str)
    with open('makefile', 'w') as makefile_file:
        makefile_file.write(makefile_str)


def test__compile_correction_potential():
    """ test varecof_io.writer.corr_potentials.compile_correction_pot
    """
    varecof_io.writer.corr_potentials.compile_corr_pot(
        MAKE_PATH)


if __name__ == '__main__':
    test__corr_potentials_writer()
    test__compile_correction_potential()
