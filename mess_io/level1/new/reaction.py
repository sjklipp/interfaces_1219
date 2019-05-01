"""
Builds a MESS input for a reaction
"""

# Writes main sections of the MESS input
from mess_io.writer import write_global_pf as globalpf
from mess_io.writer import write_global_reaction as globalrxn
from mess_io.writer import write_energy_transfer as etransfer
from mess_io.writer import write_species as species
from mess_io.writer import write_well as well
from mess_io.writer import write_bimolecular as bimolecular
from mess_io.writer import write_ts_sadpt as tssadpt
from mess_io.writer import write_molecule as molecule
from mess_io.writer import write_atom as atom
# Various strings to write
import mess_io.writer.stringslib.RXN_CHAN_HEAD_STR 
import mess_io.writer.stringslib.SPECIES_HEAD_STR 
import mess_io.writer.tringslib.SPECIES_SEC_SEP_STR
# Helper build functions
from util import build_core
from util import build_hr
# Helps obtain data from paths to files
from util import energy_from_path
from util import geom_from_path
from util import freqs_from_path
from util import hr_from_path


# Add dictionary to call to clean up the labeling and the file paths?

# Initialize input file name and empty file input string
inp_str = ''
mess_file_name = 'messrxn.inp'

# Write the global keys section
input_str += globalrxn(temperatures=[200, 300, 400],
                       pressures=[0.1, 1.0, 10.0])

# Write the energy transfer section
input_str += etransfer(exp_factor=150.0, 
                       exp_power=50.0, 
                       exp_cutoff=80.0,
                       eps1=100.0, 
                       eps2=200.0,
                       sig1=10.0,
                       sig2=20.0,
                       mass1=15.0, 
                       mass2=25.0)

# Writes a string for the head of a reaction channel section
input_str += RXN_CHAN_HEAD_STR 

# Set paths to the directory to file containing a well
ref_mol_path='./data/well'
mol_path='./data/well'
input_str += well(
    label='R1',
    data=molecule(
        core=build_core(
            'rigidrotor',
            geom1=geom_from_path(mol_path, 'mol.xyz'),
            sym_factor=2.000),
        freqs=freqs_from_path(mol_path, 'mol.freqs'),
        zero_energy=energy_from_path(
            ref_elec=(ref_mol_path, 'ref.ene'),
            ref_zpve=(ref_mol_path, 'ref.zpve'),
            spec1_elec=(mol_path, 'mol.ene'),
            spec1_zpve=(mol_path, 'mol.zpve')),
        elec_levels=((1, 0.0),))
)

# Writes a string for a string to seperate different species sections
input_str += SPECIES_SEC_SEP_STR 

# Set paths to the directory to file containing a bimolecular set
spec1_path='./data/bimol/s1'
spec2_path='./data/bimol/s2'
ref_mol_path='./data/bimol'
input_str += bimolecular(
    bimol_label='P1',
    species1_label='Mol1',
    species1_data=molecule(
        core=build_core(
            'rigidrotor',
            geom1=geom_from_path(spec1_path, 'mol.xyz'),
            sym_factor=3.000),
        freqs=freqs_from_path(spec1_path, 'mol.freqs'),
        zero_energy=0.0,
        elec_levels=((1, 0.0),)),
    species2_label='Atom2',
    species2_data=atom(
        name='O',
        elec_levels=((1, 0.0), (3, 150.0), (5, 450.0))),
    ground_energy=energy_from_path(
        ref_elec=(ref_mol_path, 'ref.ene'),
        ref_zpve=(ref_mol_path, 'ref.zpve'),
        spec1_elec=(spec1_path, 'mol.ene'),
        spec1_zpve=(spec1_path, 'mol.zpve'),
        spec2_elec=(spec2_path, 'atom.ene'),
        spec2_zpve=('', ''))
)

# Writes a string for a string to seperate different species sections
input_str += SPECIES_SEC_SEP_STR 

# Set paths to the directory to file containing a transition state
ts_path='./data/ts'
ref_mol_path='./data/ts'
input_str += ts_sadpt(
    ts_label='B1',
    reac_label='R1',
    prod_label='P1',
    ts_data=molecule(
        core=build_core(
            'rigidrotor',
            geom1=geom_from_path(ts_path, 'mol.xyz'),
            sym_factor=3.000),
        freqs=freqs_from_path(ts_path, 'mol.freqs'),
        zero_energy=energy_from_path(
            ref_elec=(ref_mol_path, 'ref.ene'),
            ref_zpve=(ref_mol_path, 'ref.zpve'),
            spec1_elec=(ts_path, 'mol.ene'),
            spec1_zpve=(ts_path, 'mol.zpve')),
        elec_levels=((1, 0.0),),
        hind_rot=build_hr(hr_from_path(ts_path, 'mol.hr')))
)

# Writes a final 'End' for the Model keyword
input_str += 'End'

# Writes the string to a file
mess_file_name = 'messrxn.inp'
with open(mess_file_name, 'a') as f:
    f.write(input_str)
