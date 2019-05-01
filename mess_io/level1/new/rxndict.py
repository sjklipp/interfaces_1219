wells = [
    well1 = {
        'label': 'W1',
        'mol_paths': './data/well',
        # set types for required info
        'core': rigidrotor, multirotior,
        # various true/false for info
        'hindered_rotor': True,
        'anharmonicity': False,
        'rotdistortion': False 
    },               
]

wells = [
    well1 = {
        'label': 'W1',
        'mol_path': './data/well',
        'geom_path': (ts_path, 'mol.xyz'),
        'freqs_path': (ts_path, 'mol.freqs'),
        'energy_paths': ((ref_mol_path, 'ref.ene'),
                         (ref_mol_path, 'ref.zpve'),
                         (ts_path, 'mol.ene'),
                         (ts_path, 'mol.zpve')),
        'hind_rot_path': (ts_path, 'mol.hr'),
        'tunnel_path': ((ts_path, 'mol.freqs'),
                        (ts_path, 'mol.ene'),
                        (reac_path, 'mol.ene'), 
                        (prod_path, 'mol.ene'))
        'anharm_path': '',
        'rovib_coups_path': '',
        'rot_dists_path': ''
    }
)

