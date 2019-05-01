"""
Builds a MESS input for a reaction
"""

import os
from lib import rxn_chan_head
from lib import species_head
from lib import species_sep
from lib import global_keys
from lib import energy_transfer
from lib import well
from lib import bimolecular
from lib import ts_sadpt
from lib import atom
from lib import molecule
from lib import core
from lib import hr
from lib import energy_from_path
from lib import geom_from_path
from lib import freqs_from_path
from lib import hr_from_path


# Set the name of the MESS input file to be created
mess_file_name = 'messrxn.inp'

# Write the global keys section
global_keys(
    filename=mess_file_name,
    messtype='reaction',
    temperatures=[200, 300, 400],
    pressures=[0.1, 1.0, 10.0]
)

# Write the energy transfer section
energy_transfer(
    filename=mess_file_name,
    exp_factor=150.0,
    exp_power=50.0,
    exp_cutoff=80.0,
    eps1=100.0,
    eps2=200.0,
    sig1=10.0,
    sig2=20.0,
    mass1=15.0,
    mass2=25.0
)

# Writes a string for the head of a reaction channel section
rxn_chan_head(
    filename=mess_file_name
)

# Set paths to the directory to file containing a well
ref_mol_path='./data/well'
mol_path='./data/well'
well(
    filename=mess_file_name,
    label='R1',
    data=molecule(
        core=core(
            'rigidrotor',
            geom1=geom_from_path(os.path.join(mol_path, 'mol.xyz')),
            sym_factor=2.000
        ),
        freqs=freqs_from_path(os.path.join(mol_path, 'mol.freqs')),
        zero_energy=energy_from_path(
            ref_elec_path=os.path.join(ref_mol_path, 'ref.ene'),
            ref_zpve_path=os.path.join(ref_mol_path, 'ref.zpve'),
            spec1_elec_path=os.path.join(mol_path, 'mol.ene'),
            spec1_zpve_path=os.path.join(mol_path, 'mol.zpve')
        ),
        elec_levels=((1, 0.0),)
    )
)

# Writes a string for a string to seperate different species sections
species_sep(
    filename=mess_file_name
)

# Set paths to the directory to file containing a bimolecular set
spec1_path='./data/bimol/s1'
spec2_path='./data/bimol/s2'
ref_mol_path='./data/bimol'
bimolecular(
    filename=mess_file_name,
    bimol_label='P1',
    species1_label='Mol1',
    species1_data=molecule(
        core=core(
            'rigidrotor',
            geom1=geom_from_path(os.path.join(spec1_path, 'mol.xyz')),
            sym_factor=3.000
        ),
        freqs=freqs_from_path(os.path.join(spec1_path, 'mol.freqs')),
        zero_energy=0.0,
        elec_levels=((1, 0.0),)
    ),
    species2_label='Atom2',
    species2_data=atom(
        name='O',
        elec_levels=((1, 0.0), (3, 150.0), (5, 450.0))
    ),
    ground_energy=energy_from_path(
        ref_elec_path=os.path.join(ref_mol_path, 'ref.ene'),
        ref_zpve_path=os.path.join(ref_mol_path, 'ref.zpve'),
        spec1_elec_path=os.path.join(spec1_path, 'mol.ene'),
        spec1_zpve_path=os.path.join(spec1_path, 'mol.zpve'),
        spec2_elec_path=os.path.join(spec2_path, 'atom.ene'),
        spec2_zpve_path=''
    )
)

# Writes a string for a string to seperate different species sections
species_sep(
    filename=mess_file_name
)

# Set paths to the directory to file containing a transition state
ts_path='./data/ts'
ref_mol_path='./data/ts'
ts_sadpt(
    filename=mess_file_name,
    ts_label='B1',
    reac_label='R1',
    prod_label='P1',
    ts_data=molecule(
        core=core(
            'rigidrotor',
            geom1=geom_from_path(ts_path, 'mol.xyz'),
            sym_factor=3.000
        ),
        freqs=freqs_from_path(ts_path, 'mol.freqs'),
        zero_energy=energy_from_path(
            ref_elec=(ref_mol_path, 'ref.ene'),
            ref_zpve=(ref_mol_path, 'ref.zpve'),
            spec1_elec=(ts_path, 'mol.ene'),
            spec1_zpve=(ts_path, 'mol.zpve')
        ),
        elec_levels=((1, 0.0),),
        hind_rot=hr(hr_from_path(ts_path, 'mol.hr'))
    )
)

# Writes a final 'End' for the Model keyword
with open(mess_file_name, 'a') as f:
    f.write('End')
