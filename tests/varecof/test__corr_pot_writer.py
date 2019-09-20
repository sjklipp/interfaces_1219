"""
  Tests the varecof_io.writer functions
"""

import varecof_io.writer

NPOT = 5
BND_LBL = ['C', 'O']
BND_IDX = [1, 3]
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


def test__potentials_writer():
    """ tests varecof_io.writer.potentials.species
              varecof_io.writer.potentials.dummy
              varecof_io.writer.potentials.auxiliary
              varecof_io.writer.potentials.makefile
    """
    
    # Write the species_corr.f string
    species_corr_str = varecof_io.writer.potentials.species()
    print('\n\nspecies_corr.f:')
    print(species_corr_str)
    
    # Write the dummy_corr.f string
    dummy_corr_str = varecof_io.writer.potentials.dummy()
    print('\n\ndummy_corr.f:')
    print(dummy_corr_str)
    
    # Write the pot_aux.f string
    auxiliary_str = varecof_io.writer.potentials.auxiliary()
    print('\n\npot_aux.f:')
    print(auxiliary_str)
    
    # Write the makefile string
    makefile_str = varecof_io.writer.potentials.makefile()
    print('\n\nmakefile:')
    print(makefile_str)


def test__compile_correction_potential():
    """ test varecof_io.writer.potentials.compile_correction_pot
    """
    varecof_io.writer.potentials.compile_correction_pot(PATH)


if __name__ == '__main__':
    test__potentials_writer()
    test__compile_correction_potential()
