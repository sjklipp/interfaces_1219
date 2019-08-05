"""
 Generates NASA Polynomial from MESS+THERMP+PAC99 outputs
"""

from . import util


def get_pac99_polynomial(output_string):
    """ read in the polyn
    """

    lines = output_string.splitlines()
    pac_polynomial = '\n'.join([lines[i] for i in range(11)])

    return pac_polynomial


def convert_pac_to_chemkin(pac_poly_str):
    """ convert the polynimal from pac format to chemkin polynomial
    """

    # Parse the lines of the pac string containing the desired coefficients
    lines = pac_poly_str.splitlines()
    las = [0.0 for i in range(7)]
    has = [0.0 for i in range(7)]
    las[0:5] = util.parse_line16(lines[6][0:80])
    las[5:7] = util.parse_line16(lines[7][48:80])
    has[0:5] = util.parse_line16(lines[9][0:80])
    has[5:7] = util.parse_line16(lines[10][48:80])

    # Build a string for the NASA polynomial in ChemKin format
    line2 = "% 15.8E% 15.8E% 15.8E% 15.8E% 15.8E    2\n"%(has[0], has[1], has[2], has[3], has[4])
    line3 = "% 15.8E% 15.8E% 15.8E% 15.8E% 15.8E    3\n"%(has[5], has[6], las[0], las[1], las[2])
    line4 = "% 15.8E% 15.8E% 15.8E% 15.8E                   4\n"%(las[3], las[4], las[5], las[6])
    poly_str = line2 + line3 + line4

    return poly_str
