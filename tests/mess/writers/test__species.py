""" Tests the writing of the energy transfer section
"""

import mess_io.writers


def test__atom_writer():
    """ Writes a string containing all info for an atom in MESS style
    """

    # Set the name and electronic levels for the atom
    atom_name = 'O'
    atom_elec_levels = ((1, 0.00), (3, 150.0), (9, 450.0))

    # Use the writer to create a string for the atom section
    atom_section_str = mess_io.writers.write_atom(atom_name, atom_elec_levels)

    # Print the atom section string
    print(atom_section_str)


def test__molecule_writer():
    """ Writes a string containing all info for a molecule in MESS style
    """

    # Set the information for a molecule
    mol_geom = (('O', (1.911401284, 0.16134481659, -0.05448080419)),
                ('N', (4.435924209, 0.16134481659, -0.05448080419)),
                ('N', (6.537299661, 0.16134481659, -0.05448080419)))
    mol_symfactor = 1.000
    mol_freqs = (100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0)
    mol_elec_levels = ((1, 0.0), (3, 50.0))
    mol_zero_e = -35.0

    # Get the string for the core using the geometry
    mol_core = mess_io.writers.write_core_rigidrotor(mol_geom, mol_symfactor, interp_emax='')

    # Use the writer to create a string for the molecule section
    molecule_section_str = mess_io.writers.write_molecule(mol_core, mol_freqs, mol_zero_e, mol_elec_levels,
                                                          hind_rot='', anharm='', tunnel='')

    # Print the molecule section string
    print(molecule_section_str)


if __name__ == '__main__':
    test__atom_writer()
    test__molecule_writer()
