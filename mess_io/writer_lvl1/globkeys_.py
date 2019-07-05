"""
mess global keywords string
"""

import mess_io
from mess_io.writer_lvl1.params import MESS


def build_mess_global_keys_str(glob_keys_dict):
    """ builds the mess input string for the global keyword sections
    """

    # Check the MESS run type 
    assert glob_keys_dict[MESS.RTYP] == 'reaction' or glob_keys_dict[MESS.RTYP] == 'pf'

    # Checks if temperatures are declared as required
    assert MESS.TEMP in glob_keys_dict

    # Run reaction
    if glob_keys_dict[MESS.RTYP] == 'reaction':
        
        assert MESS.PRES in glob_keys_dict
        # NEED: Set the other global keywords if not set
        global_keys_str = mess_io.writer.write_global_reaction(glob_keys_dict[MESS.TEMP], glob_keys_dict[MESS.PRES])

    else:
   
        # NEED: Check if rel temp_inc and atom dist min in keys    
        # Use the writer to create a string for the global keyword section for messpf
        global_keys_str = mess_io.writer.write_global_pf(temps, rel_temp_inc=0.001, atom_dist_min=0.6)

    return global_keys_str
