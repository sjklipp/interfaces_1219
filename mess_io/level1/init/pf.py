"""
Builds a MESS input for a messpf run
"""

from lib import global_keys
from lib import species_head
from lib import species
from lib import molecule
from lib import core


# Set the name of the MESS input file to be created
mess_file_name = 'messpf.inp'

# Write the global keys section
global_keys(
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
    data=molecule(
        core=core(
            'rigidrotor',
            geom1=(('O', (1.911401284, 0.16134481659, -0.05448080419)),
                   ('N', (2.435924209, 0.16134481659, -0.05448080419)),
                   ('N', (3.537299661, 0.16134481659, -0.05448080419))),
            sym_factor=1.000
        ),
        freqs=(100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0),
        zero_energy=-35.0,
        elec_levels=((1, 0.0), (3, 50.0))
    )
)

# Writes a final 'End' for the Model keyword
with open(mess_file_name, 'a') as f:
    f.write('End')
