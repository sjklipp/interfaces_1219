""" Tests the writing of the species section
"""

import mess_io.writer


MOLECULE_MESS_STRING = """RRHO
  Core RigidRotor
    SymmetryFactor          1.0
  End
  Geometry[angstrom]        3
    O   1.911401284  0.16134481659  -0.05448080419
    N   4.435924209  0.16134481659  -0.05448080419
    N   6.537299661  0.16134481659  -0.05448080419
  Frequencies[1/cm]         9
    100.00  200.00  300.00  400.00  500.00
    600.00  700.00  800.00  900.00
  ZeroEnergy[kcal/mol]      -35.0
  ElectronicLevels[1/cm]    2
    1  0.0
    3  50.0
End
"""

ATOM_MESS_STRING = """Atom
  Name O
  ElectronicLevels[1/cm]    3
    1  0.0
    3  150.0
    9  450.0
"""

def test__species_writer():
    """ Writes the MESS input for a Well
    """

    # Set a label for the well
    species_label = 'TEST'

    # Set the data string to the global molecule section
    species_data = MOLECULE_MESS_STRING

    # Use the writer to create a string for well section
    species_section_str = mess_io.writer.write_species(species_label, species_data)

    # Print the well section string
    print(species_section_str)


def test__well_writer():
    """ Writes the MESS input for a Well
    """

    # Set a label for the well
    well_label = 'W1'

    # Set the data string to the global molecule section
    well_data = MOLECULE_MESS_STRING

    # Use the writer to create a string for well section
    well_section_str = mess_io.writer.write_well(well_label, well_data)

    # Print the well section string
    print(well_section_str)


def test__bimolecular_writer():
    """ Writes the MESS input for a bimolecular set
    """

    # Set a label for the bimolecular set
    bimol_label = 'R1'

    # Set labels for the two species in the bimolecular set
    species1_label = 'Mol1'
    species2_label = 'Mol2'

    # Set the data strings to the global atom and molecule strings
    species1_data = ATOM_MESS_STRING
    species2_data = MOLECULE_MESS_STRING

    # Set the ground energy variable
    ground_energy = 50.0

    # Use the writer to create a string for the molecule section
    bimolecular_str = mess_io.writer.write_bimolecular(bimol_label,
                                                       species1_label, species1_data,
                                                       species2_label, species2_data,
                                                       ground_energy)

    # Print the bimol section string
    print(bimolecular_str)


def test__ts_sadpt_writer():
    """ ts sadpt writer
    """

    # Set the data string to the global molecule section
    ts_data = MOLECULE_MESS_STRING

    # Set labels for TS
    ts_label = 'B1'
    reac_label = 'R1'
    prod_label = 'P1'

    # Use the writer to create a string for the ts sadpt section
    ts_sadpt_str = mess_io.writer.write_ts_sadpt(ts_label, reac_label, prod_label, ts_data)

    # Print the ts sadpoint section
    print(ts_sadpt_str)


def test__ts_irc_writer():
    """ ts irc writer
    """

    # Set the number of points along the irc
    nirc = 21

    # Loop over all the points of the irc and build MESS strings
    irc_pt_strings = []
    for i in range(21):
        irc_pt_string = '! IRC Point {0}\n'.format(str(i+1))
        irc_pt_string += MOLECULE_MESS_STRING
        irc_pt_strings.append(irc_pt_string)

    # Set labels for TS
    ts_label = 'B1'
    reac_label = 'R1'
    prod_label = 'P1'

    # Use the writer to create a string for the ts irc section
    ts_irc_str = mess_io.writer.write_ts_irc(ts_label, reac_label, prod_label, irc_pt_strings)

    # Print the ts sadpoint section
    print(ts_irc_str)


if __name__ == '__main__':
    test__species_writer()
    test__well_writer()
    test__bimolecular_writer()
    test__ts_sadpt_writer()
    test__ts_irc_writer()
