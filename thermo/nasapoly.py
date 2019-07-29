"""
 Generates NASA Polynomial from MESS+THERMP+PAC99 outputs
"""

import os
from . import util


def get_pac99_polynomial(output_string):
    """ read in the polyn
    """
   
    output_lines = output_string.splitlines()
    for i in range(len(output_lines)):
        if 'THERMODYNAMIC DATA COEFFICIENTS, RECORD IMAGES' in output_lines[i]:
            block_start = i
            break    

    pac_lines = [line for line in output_lines[block_start+2: block_start+13]]

    pac_polynomial = '\n'.join(pac_lines)

    return pac_polynomial


def convert_pac_to_chemkin(pac_poly_str):
    """ convert the polynimal from pac format to chemkin polynomial
    """

    # Parse the lines of the pac string containing the desired coefficients
    lines = pac_poly_str.splitlines()
    las  = [0.0 for i in range(7)]
    has  = [0.0 for i in range(7)]
    las[0:5] = util.parse_line16(lines[6][1:81])
    las[5:7] = util.parse_line16(lines[7][49:81])
    has[0:5] = util.parse_line16(lines[9][1:81])
    has[5:7] = util.parse_line16(lines[10][49:81])

    # Build a string for the NASA polynomial in ChemKin format
    line2 = "% 15.8E% 15.8E% 15.8E% 15.8E% 15.8E    2\n"%(has[0], has[1], has[2], has[3], has[4])
    line3 = "% 15.8E% 15.8E% 15.8E% 15.8E% 15.8E    3\n"%(has[5], has[6], las[0], las[1], las[2])
    line4 = "% 15.8E% 15.8E% 15.8E% 15.8E                   4\n"%(las[3], las[4], las[5], las[6])
    poly_str = line2 + line3 + line4

    return poly_str
