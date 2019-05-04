"""
Python dictionaries containing the information to be passed into the MESS writers
"""

from mess_io.writer_lvl1 import MESS
from mess_io.writer_lvl1 import build_mess_global_keys_str
from mess_io.writer_lvl1 import build_mess_energy_transfer_str
from mess_io.writer_lvl1 import build_mess_species_str
from mess_io.writer_lvl1 import build_mess_wells_str
from mess_io.writer_lvl1 import build_mess_bimols_str
from mess_io.writer_lvl1 import build_mess_ts_sadpt_str


# Dictionary to define the global keys
# globkeys = { }
glob_keys = {
    MESS.RTYP: 'reaction',
    MESS.TEMP: [100, 200, 300, 400, 500],
    MESS.PRES: [0.1, 1.0, 10.0]
}

# Dictionary to define the energy transfer keys
# energy_transfer = { }
energy_transfer = {    
    MESS.EXP_FCT: 150.0, 
    MESS.EXP_POW: 50.0, 
    MESS.EXP_CUT: 80.0,
    MESS.EPS1: 100.0, 
    MESS.EPS2: 200.0,
    MESS.SIG1: 10.0,
    MESS.SIG2: 20.0,
    MESS.MASS1: 15.0, 
    MESS.MASS2: 25.0
}

# List of Species Dictionaries
# species = [ ({}, {}), ]
species = [
    ({MESS.MOL_LBL: 'S1',
      MESS.CORE: 'rigidrotor',
      MESS.GEOM: (('O', (1.911401284, 0.16134481659, -0.05448080419)),
                  ('N', (2.435924209, 0.16134481659, -0.05448080419)),
                  ('N', (3.537299661, 0.16134481659, -0.05448080419))),
      MESS.SYMF: 1.000,
      MESS.FREQ: (100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0),
      MESS.ZENE: -35.0,
      MESS.ELVL: ((0.0, 1), (50.0, 3))}, {}),
    ({MESS.MOL_LBL: 'S1',
      MESS.CORE: 'rigidrotor',
      MESS.GEOM: (('O', (1.911401284, 0.16134481659, -0.05448080419)),
                  ('N', (2.435924209, 0.16134481659, -0.05448080419)),
                  ('N', (3.537299661, 0.16134481659, -0.05448080419))),
      MESS.SYMF: 1.000,
      MESS.FREQ: (100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0),
      MESS.ZENE: -35.0,
      MESS.ELVL: ((0.0, 1), (50.0, 3))}, {})
]

# List of Well Dictionaries
wells = [
    # wells = [ ({}, {}), ]
    # EGDN-NO2            
    ({MESS.MOL_LBL: 'W1',
      MESS.CORE: 'rigidrotor',
      MESS.SYMF: 1.000,
      MESS.ELVL: [[0.0, 1]],
      MESS.GEOM: 'PATH',
      MESS.FREQ: 'PATH',
      MESS.ZENE: 'PATH'},
     {MESS.MOL_PATH: './data/well',
      MESS.REF_PATH: './data/well',
      MESS.GEOM_PATH: 'mol.xyz',
      MESS.FREQ_PATH: 'mol.freqs',
      MESS.MOL_ENE_PATH: 'mol.ene',
      MESS.MOL_ZPVE_PATH: 'mol.zpve',
      MESS.REF_ENE_PATH: 'ref.ene',
      MESS.REF_ZPVE_PATH: 'ref.zpve'})
]

# List of Bimolecular Dictionaries
# bimolecular = [ ( {}, ({}, {}), ({}, {}) ), ]
bimols = [
     ({MESS.BI_LBL: 'R1',
       MESS.GENE: 'PATH',
       MESS.REF_PATH: './data/bimol',
       MESS.REF_ENE_PATH: 'ref.ene',
       MESS.REF_ZPVE_PATH: 'ref.zpve'},
      ({MESS.MOL_LBL: 'S1',
        MESS.MTYP: 'molecule',
        MESS.CORE: 'rigidrotor',
        MESS.GEOM: (('O', (1.911401284, 0.16134481659, -0.05448080419)),
                    ('N', (2.435924209, 0.16134481659, -0.05448080419)),
                    ('N', (3.537299661, 0.16134481659, -0.05448080419))),
        MESS.SYMF: 1.000,
        MESS.FREQ: (100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0),
        MESS.ZENE: '0.0',
        MESS.ELVL: ((1, 0.0), (3, 50.0))},
       {MESS.MOL_PATH: './data/bimol',
        MESS.MOL_ENE_PATH: 'mol1.ene',
        MESS.MOL_ZPVE_PATH: 'mol1.zpve'}),
      ({MESS.MOL_LBL: 'S2',
        MESS.MTYP: 'molecule',
        MESS.CORE: 'rigidrotor',
        MESS.GEOM: (('O', (1.911401284, 0.16134481659, -0.05448080419)),
                    ('N', (2.435924209, 0.16134481659, -0.05448080419)),
                    ('N', (3.537299661, 0.16134481659, -0.05448080419))),
        MESS.SYMF: 1.000,
        MESS.FREQ: (100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0),
        MESS.ZENE: 0.0,
        MESS.ELVL: ((0.0, 1), (50.0, 3))},
       {MESS.MOL_PATH: './data/bimol',
        MESS.MOL_ENE_PATH: 'mol2.ene',
        MESS.MOL_ZPVE_PATH: 'mol2.zpve'}),
    )
]

# List of Saddle-Point Transition State Dictionaries
# tras_states = [ ({}, ({}, {}), ]
trans_states = [
     ({MESS.TS_LBL: 'B1',
       MESS.R_LBL: 'R1',
       MESS.P_LBL: 'P1',
       MESS.MTYP: 'molecule',
       MESS.CORE: 'rigidrotor',
       MESS.GEOM: (('O', (1.911401284, 0.16134481659, -0.05448080419)),
                   ('N', (2.435924209, 0.16134481659, -0.05448080419)),
                   ('N', (3.537299661, 0.16134481659, -0.05448080419))),
       MESS.SYMF: 1.000,
       MESS.FREQ: (100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0),
       MESS.ZENE: '0.0',
       MESS.ELVL: ((1, 0.0), (3, 50.0)),
       MESS.HROT: 'PATH',
       MESS.ANHM: [ [1.00, 2.00, 3.00], [1.00, 2.00, 3.00], [4.00, 5.00, 6.00] ],
       MESS.RVCP: (-0.01035, -0.01035),
       MESS.RDIS: (('aaaa', -0.2623325156e-04),
                   ('bbaa', -0.2623325156e-04),
                   ('bbbb', -0.2623325156e-04))},
       {MESS.MOL_PATH: './data/ts',
       MESS.MOL_ENE_PATH: 'ts.ene',
       MESS.MOL_ZPVE_PATH: 'ts.zpve',
       MESS.FREQ_PATH: 'mol.freqs',
       MESS.REAC_PATH: './data/ts',
       MESS.REAC_ENE_PATH: 'reac.ene',
       MESS.REAC_ZPVE_PATH: 'reac.zpve',
       MESS.PROD_PATH: './data/ts',
       MESS.PROD_ENE_PATH: 'prod.ene',
       MESS.PROD_ZPVE_PATH: 'prod.zpve',
       MESS.HROT_PATH: 'mol.hr'}),
]


if __name__ == '__main__':
    mess_str = build_mess_global_keys_str(glob_keys)
    mess_str += build_mess_energy_transfer_str(energy_transfer)
    mess_str += build_mess_species_str(species)
    mess_str += build_mess_wells_str(wells)
    mess_str += build_mess_bimols_str(bimols)
    mess_str += build_mess_ts_sadpt_str(trans_states)
    mess_str += 'End'
    with open('mess.inp', 'w') as messfile:
        messfile.write(mess_str)

