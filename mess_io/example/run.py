"""
Python dictionaries containing the information to be passed into the MESS writers
"""

import copy
from mess_io.writer_lvl1 import MESS
from mess_io.writer_lvl1 import build_mess_global_keys_str
from mess_io.writer_lvl1 import build_mess_energy_transfer_str
from mess_io.writer_lvl1 import build_mess_ts_str
from mess_io.writer_lvl1 import build_mess_wells_str
from mess_io.writer_lvl1 import build_mess_bimols_str
from sdicts import * 

# Set dirs to min and TS
MIN_DIR = '/lcrc/project/CMRP/mechanisms/info/1_EGDN/1_Min'
TS_DIR = '/lcrc/project/CMRP/mechanisms/info/1_EGDN/2_TS'
WELL_DIR = '/lcrc/project/CMRP/mechanisms/info/1_EGDN/3_Wells'

# Set the head dicts reference paths
ref_path_dict = {MESS.REF_PATH: MIN_DIR+'/01_EGDN'}
ts_dict.update(ref_path_dict)
well_dict.update(ref_path_dict)
bimol_dict.update(ref_path_dict)

# Global keyword section
glob_keys_dict = {
    MESS.RTYP: 'reaction',
    MESS.TEMP: [500, 600, 700, 800, 900, 1000, 
                1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 
                2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000], 
    MESS.PRES: [0.1, 1.0, 10.0, 100.0, 1000.0]
}

# Energy-transfer keyword section
energy_trans_dict = {    
    MESS.EXP_FCT: 300.0, 
    MESS.EXP_POW: 0.40, 
    MESS.EXP_CUT: 10.0,
    MESS.EPS1: 468.748, 
    MESS.SIG1: 4.350,
    MESS.MASS1: 152.007, 
    MESS.EPS2: 468.748, 
    MESS.SIG2: 4.350,
    MESS.MASS2: 152.007, 
}

# HONO Elimination Transition State
ts_c2_dict = copy.deepcopy(ts_dict)
ts_c2_data_dict = copy.deepcopy(mol_data_dict)
ts_c2_path_dict = copy.deepcopy(mol_path_dict)
ts_c2_dict.update({MESS.HD_LBL: 'B21'})
ts_c2_dict.update({MESS.P_LBL: 'W2'})
ts_c2_dict.update({MESS.TS_TYP: 'sadpt'})
ts_c2_data_dict.update({MESS.HROT: 'PATH'})
ts_c2_data_dict.update({MESS.TUNL: ['ECKART', [1467.80, 42.529, 66.302]]})
ts_c2_path_dict.update({MESS.MOL_PATH: TS_DIR+'/3_Elim/EGDN/HONO'})
ts_c2 = [
     (ts_c2_dict, ts_c2_data_dict, ts_c2_path_dict)
]
# HONO2 Elimination Transition State
ts_c3_dict = copy.deepcopy(ts_dict)
ts_c3_data_dict = copy.deepcopy(mol_data_dict)
ts_c3_path_dict = copy.deepcopy(mol_path_dict)
ts_c3_dict.update({MESS.HD_LBL: 'B31'})
ts_c3_dict.update({MESS.P_LBL: 'W3'})
ts_c3_dict.update({MESS.TS_TYP: 'sadpt'})
ts_c3_data_dict.update({MESS.HROT: 'PATH'})
ts_c3_data_dict.update({MESS.TUNL: ['ECKART', [1339.57, 48.084, 37.938]]})
ts_c3_path_dict.update({MESS.MOL_PATH: TS_DIR+'/3_Elim/EGDN/HONO2'})
ts_c3 = [
     (ts_c3_dict, ts_c3_data_dict, ts_c3_path_dict)
]

# EGDN Reactant Well
wellr_dict = copy.deepcopy(well_dict)
wellr_data = copy.deepcopy(mol_data_dict)
wellr_path = copy.deepcopy(mol_path_dict)
wellr_data.update({MESS.HROT: 'PATH'})
wellr_path.update({MESS.MOL_PATH: MIN_DIR+'/01_EGDN'})
wells_r = [
    (wellr_dict, wellr_data, wellr_path)
]
# EGDN-HONO...HONO Exit-Channel Well
well_c2_dict = copy.deepcopy(well_dict)
well_c2_data = copy.deepcopy(mol_data_dict)
well_c2_path = copy.deepcopy(mol_path_dict)
well_c2_dict.update({MESS.HD_LBL: 'W2'})
well_c2_path.update({MESS.FREQ_PATH: 'freq/1_NonProj/mol.harmfreq'})
well_c2_path.update({MESS.MOL_PATH: WELL_DIR+'/egdn_hono'})
wells_c2 = [
    (well_c2_dict, well_c2_data, well_c2_path),
]
# EGDN-HONO2...HONO2 Exit-Channel Well
well_c3_dict = copy.deepcopy(well_dict)
well_c3_data = copy.deepcopy(mol_data_dict)
well_c3_path = copy.deepcopy(mol_path_dict)
well_c3_dict.update({MESS.HD_LBL: 'W3'})
well_c3_path.update({MESS.FREQ_PATH: 'freq/1_NonProj/mol.harmfreq'})
well_c3_path.update({MESS.MOL_PATH: WELL_DIR+'/egdn_hono2'})
wells_c3 = [
    (well_c3_dict, well_c3_data, well_c3_path)
]

# EGDN-NO2 + NO2 Products
bimol_c1_dict = copy.deepcopy(bimol_dict)
bimol_c1_data_1 = copy.deepcopy(mol_data_dict)
bimol_c1_path_1 = copy.deepcopy(mol_path_dict)
bimol_c1_dict.update({MESS.HD_LBL: 'P1'})
bimol_c1_data_1.update({MESS.MOL_LBL: 'EGDN-NO2'})
bimol_c1_data_1.update({MESS.HROT: 'PATH'})
bimol_c1_data_1.update({MESS.ZENE: 0.00})
bimol_c1_data_1.update({MESS.ELVL: [[0.0, 2]]})
bimol_c1_path_1.update({MESS.MOL_PATH: MIN_DIR+'/02_EGDN-NO2'})
bimol_c1_data_2 = copy.deepcopy(mol_data_dict)
bimol_c1_path_2 = copy.deepcopy(mol_path_dict)
bimol_c1_data_2.update({MESS.MOL_LBL: 'NO2'})
bimol_c1_data_2.update({MESS.ZENE: 0.00})
bimol_c1_data_2.update({MESS.ELVL: [[0.0, 2]]})
bimol_c1_path_2.update({MESS.FREQ_PATH: 'freq/1_NonProj/mol.harmfreq'})
bimol_c1_path_2.update({MESS.MOL_PATH: MIN_DIR+'/13_NO2/wb'})
bimols_c1 = [
     (bimol_c1_dict, 
      (bimol_c1_data_1, bimol_c1_path_1),
      (bimol_c1_data_2, bimol_c1_path_2)),
]
# EGDN-HONO + HONO Products
bimol_c2_dict = copy.deepcopy(bimol_dict)
bimol_c2_data_1 = copy.deepcopy(mol_data_dict)
bimol_c2_path_1 = copy.deepcopy(mol_path_dict)
bimol_c2_dict.update({MESS.HD_LBL: 'P2'})
bimol_c2_data_1.update({MESS.MOL_LBL: 'EGDN-HONO'})
bimol_c2_data_1.update({MESS.HROT: 'PATH'})
bimol_c2_data_1.update({MESS.ZENE: 0.00})
bimol_c2_path_1.update({MESS.MOL_PATH: MIN_DIR+'/04_EGDN-HONO'})
bimol_c2_data_2 = copy.deepcopy(mol_data_dict)
bimol_c2_path_2 = copy.deepcopy(mol_path_dict)
bimol_c2_data_2.update({MESS.MOL_LBL: 'HONO'})
bimol_c2_data_2.update({MESS.HROT: 'PATH'})
bimol_c2_data_2.update({MESS.ZENE: 0.00})
bimol_c2_path_2.update({MESS.MOL_PATH: MIN_DIR+'/11_HONO/wb'})
bimols_c2 = [
     (bimol_c2_dict, 
      (bimol_c2_data_1, bimol_c2_path_1),
      (bimol_c2_data_2, bimol_c2_path_2)),
]
# EGDN-HONO2 + HONO2 Products
bimol_c3_dict = copy.deepcopy(bimol_dict)
bimol_c3_data_1 = copy.deepcopy(mol_data_dict)
bimol_c3_path_1 = copy.deepcopy(mol_path_dict)
bimol_c3_dict.update({MESS.HD_LBL: 'P3'})
bimol_c3_data_1.update({MESS.MOL_LBL: 'EGDN-HONO2'})
bimol_c3_data_1.update({MESS.HROT: 'PATH'})
bimol_c3_data_1.update({MESS.ZENE: 0.00})
bimol_c3_path_1.update({MESS.MOL_PATH: MIN_DIR+'/05_EGDN-HONO2'})
bimol_c3_data_2 = copy.deepcopy(mol_data_dict)
bimol_c3_path_2 = copy.deepcopy(mol_path_dict)
bimol_c3_data_2.update({MESS.MOL_LBL: 'HONO2'})
bimol_c3_data_2.update({MESS.HROT: 'PATH'})
bimol_c3_data_2.update({MESS.ZENE: 0.00})
bimol_c3_path_2.update({MESS.MOL_PATH: MIN_DIR+'/12_HONO2'})
bimols_c3 = [
     (bimol_c3_dict, 
      (bimol_c3_data_1, bimol_c3_path_1),
      (bimol_c3_data_2, bimol_c3_path_2)),
]

# Build the MESS input
mess_str = build_mess_global_keys_str(glob_keys_dict)
mess_str += build_mess_energy_transfer_str(energy_trans_dict)
mess_str += '!+++++++++++++++++++++++++++++++++++++++++++++++\n'
mess_str += '! EGDN REACTANT\n'
mess_str += '!+++++++++++++++++++++++++++++++++++++++++++++++\n'
mess_str += build_mess_wells_str(wells_r)
mess_str += '!+++++++++++++++++++++++++++++++++++++++++++++++\n'
mess_str += '! O-N FISSION CHANNEL\n'
mess_str += '!+++++++++++++++++++++++++++++++++++++++++++++++\n'
mess_str += build_mess_bimols_str(bimols_c1)
mess_str += '!+++++++++++++++++++++++++++++++++++++++++++++++\n'
mess_str += '! HONO ELIMINATION CHANNEL\n'
mess_str += '!+++++++++++++++++++++++++++++++++++++++++++++++\n'
mess_str += build_mess_ts_str(ts_c2)
mess_str += build_mess_wells_str(wells_c2)
mess_str += build_mess_bimols_str(bimols_c2)
mess_str += '!+++++++++++++++++++++++++++++++++++++++++++++++\n'
mess_str += '! HONO2 ELIMINATION CHANNEL ++++++++++\n'
mess_str += '!+++++++++++++++++++++++++++++++++++++++++++++++\n'
mess_str += build_mess_ts_str(ts_c3)
mess_str += build_mess_wells_str(wells_c3)
mess_str += build_mess_bimols_str(bimols_c3)
with open('init.inp', 'w') as messfile:
    messfile.write(mess_str)
