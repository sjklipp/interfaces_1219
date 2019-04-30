"""
Builds a MESS input for a messpf run
"""

import os
from lib import globalkeys
from lib import species_head
from lib import species
from lib import molecule_from_data


# Set the name of the MESS input file to be created
mess_file_name = 'pf.inp'

# Write the global keys section
globalkeys(
    filename=mess_file_name,
    messtype='pf',
    temperatures=[200,300,400]
)

# Writes a string for the head of a species section
species_head(
    filename=mess_file_name
)

# Writes the data component of a species section
species(
    filename=mess_file_name,
    label='S1',
    data=molecule_from_data(
        core='rigidrotor',
        geom=(('O', (1.911401284, 0.16134481659, -0.05448080419)),
              ('N', (4.435924209, 0.16134481659, -0.05448080419)),
              ('N', (6.537299661, 0.16134481659, -0.05448080419))),
        symfactor=1.000,
        freqs=(100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0),
        elec_levels=((1, 0.0), (3, 50.0)),
        zero_e=-35.0,
    )
)
