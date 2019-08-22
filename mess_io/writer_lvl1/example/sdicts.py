""" 
Base dictionary definitions for cleaning higher level MESS input calls. 
These can be used to input to level 1 writer calls.
"""

from mess_io.writer_lvl1 import MESS


# Well header dictionary 
well_dict = {MESS.HD_LBL: 'R1',
             MESS.REF_ENE_PATH: 'energy/ccf12_tzf12/mol.ene',
             MESS.REF_ZPVE_PATH: 'freq/mol.harmzpve'}

# Bimolecular header dictionary 
bimol_dict = {MESS.HD_LBL: 'P1',
              MESS.GENE: 'PATH',
              MESS.REF_ENE_PATH: 'energy/ccf12_tzf12/mol.ene',
              MESS.REF_ZPVE_PATH: 'freq/mol.harmzpve'}

# TS header dictionary 
ts_dict = {MESS.HD_LBL: 'B1',
           MESS.R_LBL: 'R1',
           MESS.P_LBL: 'P1',
           MESS.REF_ENE_PATH: 'energy/ccf12_tzf12/mol.ene',
           MESS.REF_ZPVE_PATH: 'freq/mol.harmzpve'}

# Molecule data dictionary
mol_data_dict = {MESS.CORE: 'rigidrotor',
                 MESS.MTYP: 'molecule',
                 MESS.SYMF: 1.000,
                 MESS.ELVL: [[0.0, 1]],
                 MESS.GEOM: 'PATH',
                 MESS.FREQ: 'PATH',
                 MESS.ZENE: 'PATH'}

# Molecule path dictionary
mol_path_dict = {MESS.GEOM_PATH: 'geom/Final/geom.xyz',
                 MESS.FREQ_PATH: 'freq/2_Proj/hrproj_freq.dat',
                 MESS.HROT_PATH: 'hindrot/mol.hrpot',
                 MESS.MOL_ENE_PATH: 'energy/ccf12_tzf12/mol.ene',
                 MESS.MOL_ZPVE_PATH: 'freq/mol.harmzpve'}

# Global keyword section
glob_keys_dict = {
    MESS.RTYP: 'reaction',
    MESS.TEMP: [300, 400, 500, 600, 700, 800, 900, 1000, 
                1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 
                2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000], 
    MESS.PRES: [0.1, 1.0, 10.0, 100.0, 1000.0]
}

# Energy-transfer keyword section
energy_trans_dict = {    
    MESS.EXP_FCT: 185.0, 
    MESS.EXP_POW: 0.75, 
    MESS.EXP_CUT: 10.0,
    MESS.EPS1: 376.6, 
    MESS.EPS2: 376.6,
    MESS.SIG1: 3.5,
    MESS.SIG2: 3.5,
    MESS.MASS1: 44.013, 
    MESS.MASS2: 44.013
}
