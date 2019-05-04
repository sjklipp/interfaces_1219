"""
class to store variables for mess input data
"""

class MESS():
    
    # Global Keywords
    RTYP = 'run_type'
    TEMP = 'tempertatures'
    PRES = 'pressures'
    
    # Energy Transfer
    EXP_FCT = 'exp_factor' 
    EXP_POW = 'exp_power' 
    EXP_CUT = 'exp_cutoff'
    EPS1 = 'eps1' 
    EPS2 = 'eps2'
    SIG1 = 'sig1'
    SIG2 = 'sig2'
    MASS1 = 'mass1' 
    MASS2 = 'mass2'
    
    # Data (Required)
    CORE = 'core'
    GEOM = 'geom'
    SYMF = 'sym_factor'
    FREQ = 'freqs'
    ELVL = 'elec_levels'
    ZENE = 'zero_energy'
    GENE = 'ground_energy'
    # Data (Optional)
    HROT = 'hindered_rotor'
    ANHM = 'anharm'
    RDIS = 'rot_dists'
    RVCP = 'rovib_coups'
    TUNL = 'tunnel'
    
    # Keywords
    INTRP_EMAX = 'interp_emax'
    QUANT_LVL_EMAX = 'quant_lvl_emax'
    
    # MESS Input Labels
    MTYP = 'mol_type'
    TS_LBL = 'ts_label'
    BI_LBL = 'bimolecular_label'
    MOL_LBL = 'mol_label'
    R_LBL = 'reactant_label'
    P_LBL = 'product_label'
    
    # Head Paths
    MOL_PATH = 'spec_path'
    REAC_PATH = 'reac_path'
    PROD_PATH = 'prod_path'
    REF_PATH = 'reference_path'
    # Data Paths
    GEOM_PATH = 'geom_path'
    SYMF_PATH = 'sym_factor_path'
    FREQ_PATH = 'freqs_path'
    ELVL_PATH = 'elec_levels_path'
    HROT_PATH = 'hindered_rotor_path'
    ANHM_PATH = 'anharm_path'
    RVCP_PATH = 'rovib_coups_path'
    RDIS_PATH = 'rot_dists_path'
    # Various Elec. Energy and ZPVE Paths
    MOL_ENE_PATH = 'spec_energy_path'
    MOL_ZPVE_PATH = 'spec_zpve_path'
    REAC_ENE_PATH = 'reac_energy_path'
    REAC_ZPVE_PATH = 'reac_zpve_path'
    PROD_ENE_PATH = 'prod_energy_path'
    PROD_ZPVE_PATH = 'prod_zpve_path'
    REF_ENE_PATH = 'reference_energy_path'
    REF_ZPVE_PATH = 'reference_zpve_path'
