"""
Builds a MESS input for a reaction
"""

import os
from lib import globalkeys
from lib import species_head
from lib import well
from lib import molecule_from_path


# Set the name of the MESS input file to be created
mess_file_name = 'reaction.inp'

# Write the global keys section
globalkeys(
    filename=mess_file_name,
    messtype='reaction',
    pressures=[200,300,400],
    temperatures=[200,300,400]
)

# Writes a string for the head of a species section
species_head(
    filename=mess_file_name
)

# Set paths to the directory to file containing a well
mol_path='./data/well'
ref_mol_path='./data/well'
well(
    filename=mess_file_name,
    label='S1',
    data=molecule_from_path(
        core='rigidrotor',
        geom_path=os.path.join(mol_path, 'mol.xyz'),
        energy_path=os.path.join(mol_path, 'mol.ene'),
        ref_energy_path=os.path.join(ref_mol_path, 'ref.ene'),
        freqs_path=os.path.join(mol_path, 'mol.freq'),
        sym_factor_path=os.path.join(mol_path, 'mol.symnum'),
        elec_levels_path=os.path.join(mol_path, 'mol.levels'),
    )
    ref_energy_path=os.path.join(ref_mol_path, 'ref.ene')
)


# Set paths to the directory to file containing a bimolecular set
spec1_path='./data/bimol/s1'
spec2_path='./data/bimol/s2'
ref_mol_path='./data/bimol'
bimolecular(
    filename=mess_file_name,
    label='P1',
    species1_label='Mol1',
    species1_data=molecule_from_path(
        core='rigidrotor',
        geom_path=os.path.join(spec1_path, 'mol.xyz'),
        energy_path=os.path.join(spec1_path, 'mol.ene'),
        freqs_path=os.path.join(spec1_path, 'mol.freq'),
        sym_factor_path=os.path.join(spec1_path, 'mol.symnum'),
        elec_levels_path=os.path.join(spec1_path, 'mol.levels'),
    ),
    species2_label='Atom2',
    species2_data=molecule_from_path(
        name='O',
        energy_path=os.path.join(spec2_path, 'atom.ene'),
        elec_levels_path=os.path.join(spec2_path, 'atom.levels')
    ),
    ref_energy_path=os.path.join(ref_mol_path, 'ref.ene')
)

