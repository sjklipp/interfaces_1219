"""
uses dictionaries defined in other file to build the MESS strings
"""

import mess_io
from mess_io.writer_lvl1.params import MESS
from mess_io.writer_lvl1.molecule_ import build_mess_molecule_str


def build_mess_wells_str(wells):
    """ Builds the MESS input strings for each of the wells
    """
    
    # Initialize the wells string with the wells header
    wells_str = mess_io.writer.stringslib.WELLS_HEAD_STR

    # Loop over each well
    for well in wells:
    
        # Set dictionaries for clarity
        data_dict = well[0]
        path_dict = well[1]
        
        # SET THE VALUES FOR REQUIRED INFORMATION #
    
        # Check for existance of required keys
        assert MESS.MOL_LBL in data_dict
        assert MESS.GEOM in data_dict
        assert MESS.SYMF in data_dict
        assert MESS.FREQ in data_dict
        assert MESS.ZENE in data_dict
        assert MESS.ELVL in data_dict
        assert MESS.CORE in data_dict
       
        # Check that the core value is appropriate for wells
        assert data_dict[MESS.CORE] == 'rigidrotor' or data_dict[MESS.CORE] == 'multirotor'

        # Build the molecule string with the core string and above data
        molecule_str = build_mess_molecule_str(data_dict, path_dict)
    
        # Build the wells string with the label and molecule string
        well_str = mess_io.writer.write_well(data_dict[MESS.MOL_LBL], molecule_str)

        # Add to the wells string the well string and a seperator
        wells_str += well_str
        wells_str += mess_io.writer.stringslib.SPECIES_SEC_SEP_STR

    return wells_str
