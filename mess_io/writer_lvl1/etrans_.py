"""
mess energy transfer string
"""

import mess_io
from mess_io.writer_lvl1.params import MESS


def build_mess_energy_transfer_str(etrans_dict):
    """ builds the mess input string for the global keyword sections
    """

    # NEED: add checks to the 
    assert MESS.EXP_FCT in etrans_dict 
    assert MESS.EXP_POW in etrans_dict 
    assert MESS.EXP_CUT in etrans_dict
    assert MESS.EPS1 in etrans_dict    
    assert MESS.EPS2 in etrans_dict      
    assert MESS.SIG1 in etrans_dict    
    assert MESS.SIG2 in etrans_dict      
    assert MESS.MASS1 in etrans_dict    
    assert MESS.MASS2 in etrans_dict  

    # Use the writer to create a string for the energy transfer
    energy_trans_str = mess_io.writer.write_energy_transfer(etrans_dict[MESS.EXP_FCT],  
                                                            etrans_dict[MESS.EXP_POW], 
                                                            etrans_dict[MESS.EXP_CUT], 
                                                            etrans_dict[MESS.EPS1], 
                                                            etrans_dict[MESS.EPS2], 
                                                            etrans_dict[MESS.SIG1],  
                                                            etrans_dict[MESS.SIG2],    
                                                            etrans_dict[MESS.MASS1],   
                                                            etrans_dict[MESS.MASS2]) 

    return energy_trans_str
