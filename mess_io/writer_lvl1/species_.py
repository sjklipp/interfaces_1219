"""
uses dictionaries defined in other file to build the MESS strings
"""

import mess_io
from mess_io.writer_lvl1.params import MESS
from mess_io.writer_lvl1.molecule_ import build_mess_molecule_str


def build_mess_species_str(species):
    """ Builds the MESS input strings for each of the species
    """
    
    # Initialize the species string with the species header
    species_str = mess_io.writer.stringslib.SPECIES_HEAD_STR

    # Loop over each species
    for specie in species:
    
        # Set dictionaries for clarity
        data_dict = specie[0]
        path_dict = specie[1]
        
        # SET THE VALUES FOR REQUIRED INFORMATION #
    
        # Check for existance of required keys
        assert MESS.MOL_LBL in data_dict
        assert MESS.GEOM in data_dict
        assert MESS.SYMF in data_dict
        assert MESS.FREQ in data_dict
        assert MESS.ZENE in data_dict
        assert MESS.ELVL in data_dict
        assert MESS.CORE in data_dict
       
        # Check that the core value is appropriate for species
        assert data_dict[MESS.CORE] == 'rigidrotor' or data_dict[MESS.CORE] == 'multirotor'

        # Build the molecule string with the core string and above data
        molecule_str = build_mess_molecule_str(data_dict, path_dict)
    
        # Build the species string with the label and molecule string
        specie_str = mess_io.writer.write_species(data_dict[MESS.MOL_LBL], molecule_str)

        # Add to the species string the specie string and a seperator
        species_str += specie_str
        species_str += mess_io.writer.stringslib.SPECIES_SEC_SEP_STR

    return species_str
