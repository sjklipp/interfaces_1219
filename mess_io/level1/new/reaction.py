"""
Builds a MESS input for a reaction
"""

# import string lib because I am lazy
import mess_io.writer
# Writes main sections of the MESS input
from mess_io.writer import write_global_reaction as globalrxn
from mess_io.writer import write_energy_transfer as etransfer
from mess_io.writer import write_species as species
from mess_io.writer import write_well as well
from mess_io.writer import write_bimolecular as bimolecular
from mess_io.writer import write_ts_sadpt as tssadpt
from mess_io.writer import write_molecule as molecule
from mess_io.writer import write_atom as atom
from mess_io.writer import write_tunnel_eckart as eckart
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
input_str = ''
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
input_str += mess_io.writer.stringslib.RXN_CHAN_HEAD_STR 

# Set paths to the directory to file containing a well
ref_mol_path='./data/well'
mol_path='./data/well'
input_str += well(
    well_label='R1',
    well_data=molecule(
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
        elec_levels=((0.0, 1),))
)

# Writes a string for a string to seperate different species sections
input_str += mess_io.writer.stringslib.SPECIES_SEC_SEP_STR 

# Set paths to the directory to file containing a bimolecular set
spec1_path='./data/bimol/s1'
spec2_path='./data/bimol/s2'
ref_mol_path='./data/bimol'
input_str += bimolecular(
    bimol_label = 'P1',
    species1_label = 'Mol1',
    species1_data=molecule(
        core=build_core(
            'rigidrotor',
            geom1=geom_from_path(spec1_path, 'mol.xyz'),
            sym_factor=3.000),
        freqs=freqs_from_path(spec1_path, 'mol.freqs'),
        zero_energy=0.0,
        elec_levels=((0.0, 1),)),
    species2_label = 'Atom2',
    species2_data=atom(
        name='O',
        elec_levels=((0.0, 1), (150.0, 3), (450.0, 5))),
    ground_energy=energy_from_path(
        ref_elec=(ref_mol_path, 'ref.ene'),
        ref_zpve=(ref_mol_path, 'ref.zpve'),
        spec1_elec=(spec1_path, 'mol.ene'),
        spec1_zpve=(spec1_path, 'mol.zpve'),
        spec2_elec=(spec2_path, 'atom.ene'),
        spec2_zpve=('', ''))
)

# Writes a string for a string to seperate different species sections
input_str += mess_io.writer.stringslib.SPECIES_SEC_SEP_STR 

# Set paths to the directory to file containing a transition state
ts_path='./data/ts'
reac_path='./data/ts'
prod_path='./data/ts'
ref_mol_path='./data/ts'
input_str += tssadpt(
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
        elec_levels=((0.0, 1),),
        hind_rot=build_hr(hr_from_path(ts_path, 'mol.hr')),
        tunnel='',
#        tunnel=eckart(eckart_from_path(
#            ts_freqs_path=(ts_path, 'mol.freqs')
#            ts_energy_path=(ts_path, 'mol.ene')
#            reac_energy_path=(reac_path, 'mol.ene') 
#            prod_energy_path=(prod_path, 'mol.ene')):
        anharm='',
        rovib_coups='',
        rot_dists='')
)

# Writes a string to show end of last section and last 'End for Model keyword
input_str += mess_io.writer.stringslib.SPECIES_SEC_SEP_STR 
input_str += 'End'

# Writes the string to a file
mess_file_name = 'messrxn.inp'
with open(mess_file_name, 'w') as f:
    f.write(input_str)
