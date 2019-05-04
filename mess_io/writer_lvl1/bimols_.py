"""
uses dictionaries defined in other file to build the MESS strings
"""

import mess_io
from mess_io.writer_lvl1.params import MESS
from mess_io.writer_lvl1.molecule_ import build_mess_molecule_str
from mess_io.writer_lvl1.util import geom_from_path
from mess_io.writer_lvl1.util import energy_from_path
from mess_io.writer_lvl1.util import freqs_from_path
from mess_io.writer_lvl1.util import hr_from_path
from mess_io.writer_lvl1.util import build_hr


def build_mess_bimols_str(bimols):
    """ Builds the MESS input strings for each of the species
    """
    
    # Initialize the species string with the species header
    bimols_str = mess_io.writer.stringslib.BIMOL_HEAD_STR

    # Loop over each bimolecular species
    for bimol in bimols:
    
        # Set dictionaries for clarity
        bimol_dict = bimol[0]
        spec1_data_dict = bimol[1][0]
        spec1_path_dict = bimol[1][1]
        spec2_data_dict = bimol[2][0]
        spec2_path_dict = bimol[2][1]
        
        # SET THE VALUES FOR REQUIRED INFORMATION #
            
        # Checks that the bimolecular label and grount-energy is set in the initial dictionary
        assert MESS.BI_LBL in bimol_dict
        assert MESS.GENE in bimol_dict
    
        # Check if the species1 and species2 molecule type is set to a molecule or atom
        assert spec1_data_dict[MESS.MTYP] == 'molecule' or spec1_data_dict[MESS.MTYP] == 'atom' 
        assert spec2_data_dict[MESS.MTYP] == 'molecule' or spec2_data_dict[MESS.MTYP] == 'atom' 

        # Build the species 1 string
        if spec1_data_dict[MESS.MTYP] == 'molecule':
            
            # Check for existance of required keys if species 1 is a molecule
            assert MESS.MOL_LBL in spec1_data_dict
            assert MESS.MTYP in spec1_data_dict
            assert MESS.GEOM in spec1_data_dict
            assert MESS.SYMF in spec1_data_dict
            assert MESS.FREQ in spec1_data_dict
            assert MESS.ZENE in spec1_data_dict
            assert MESS.ELVL in spec1_data_dict
            assert MESS.CORE in spec1_data_dict
            
            # Check that the zero energy is set to zero for the bimolecular set
            #assert spec1_data_dict[MESS.ZENE] == 0.0
            
            # Check that the core value is set for a molecule and is appropriate for species
            assert spec1_data_dict[MESS.CORE] == 'rigidrotor' or spec1_data_dict[MESS.CORE] == 'multirotor' 
            
            # Build the molecule string 
            spec1_str = build_mess_molecule_str(spec1_data_dict, spec1_path_dict)
        
        else:
            
            # Check for existance of required keys if species 1 is an atom
            assert MESS.NAME in spec1_data_dict
            assert MESS.ELVL in spec1_data_dict
            
            # Build the atom string 
            spec2_str = build_mess_atom_str(spec1_data_dict, spec1_path_dict)
        
        # Build the species 2 string
        if spec2_data_dict[MESS.MTYP] == 'molecule':
        
            # Check for existance of required keys if species 2 is a molecule
            assert MESS.MOL_LBL in spec2_data_dict
            assert MESS.MTYP in spec2_data_dict
            assert MESS.GEOM in spec2_data_dict
            assert MESS.SYMF in spec2_data_dict
            assert MESS.FREQ in spec2_data_dict
            assert MESS.ZENE in spec2_data_dict
            assert MESS.ELVL in spec2_data_dict
            assert MESS.CORE in spec2_data_dict
            
            # Check that the zero energy is set to zero for the bimolecular set
            #assert spec2_data_dict[MESS.ZENE] == 0.0
            
            # Check that the core value is set for a molecule and is appropriate for species
            assert spec2_data_dict[MESS.CORE] == 'rigidrotor' or spec1_data_dict[MESS.CORE] == 'multirotor' 
        
            # Build the molecule string 
            spec2_str = build_mess_molecule_str(spec2_data_dict, spec2_path_dict)

        else:
            
            # Check for existance of required keys if species 2 is an atom
            assert MESS.NAME in spec2_data_dict
            assert MESS.ELVL in spec2_data_dict
            
            # Build the atom string 
            spec2_str = build_mess_atom_str(spec2_data_dict, spec2_path_dict)

        # Ground-Energy set directly from dictionary or using paths
        if bimol_dict[MESS.GENE] == 'PATH':
            # Get the total energy for the reference, species 1, and species 2    
            ref_total_energy = energy_from_path(bimol_dict[MESS.REF_PATH], 
                                                bimol_dict[MESS.REF_ENE_PATH], 
                                                bimol_dict[MESS.REF_ZPVE_PATH])
            spec1_total_energy = energy_from_path(spec1_path_dict[MESS.MOL_PATH], 
                                                  spec1_path_dict[MESS.MOL_ENE_PATH], 
                                                  spec1_path_dict[MESS.MOL_ZPVE_PATH])
            spec2_total_energy = energy_from_path(spec2_path_dict[MESS.MOL_PATH], 
                                                  spec2_path_dict[MESS.MOL_ENE_PATH], 
                                                  spec2_path_dict[MESS.MOL_ZPVE_PATH])
            # Calculate the relative ground energy using the species energies
            ground_energy = (spec1_total_energy + spec2_total_energy) - ref_total_energy
            ground_energy *= 627.5095
        else:
            ground_energy = data_dict[MESS.GENE]
    
        # Build the bimolecular string with the label and molecule string
        bimol_str = mess_io.writer.write_bimolecular(bimol_dict[MESS.BI_LBL],
                                                     spec1_data_dict[MESS.MOL_LBL],
                                                     spec1_str, 
                                                     spec2_data_dict[MESS.MOL_LBL],
                                                     spec2_str, 
                                                     ground_energy)

        # Add to the species string the specie string and a seperator
        bimols_str += bimol_str
        bimols_str += mess_io.writer.stringslib.SPECIES_SEC_SEP_STR

        return bimols_str
