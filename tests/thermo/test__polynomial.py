"""
Tests for interacting with the polynomial
"""

import os
from thermo import nasapoly

GIT_STR = """C2H3
  3 201704 C   2.00H   3.00    0.00    0.00    0.00 0   27.0452200     296391.000
  100.000   200.000 2  0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0        10698.000
  3.587567650D+00 3.894470300D-03 0.000000000D+00 0.000000000D+00 0.000000000D+00
  0.000000000D+00 0.000000000D+00 0.000000000D+00 3.437879620D+04 6.439970150D+00
  200.000  1000.000 5  0.0  1.0  2.0  3.0  4.0  0.0  0.0  0.0        10698.000
  2.881119522D+00 4.825191250D-03 1.818030931D-05-2.828286454D-08 1.209654946D-11
  0.000000000D+00 0.000000000D+00 0.000000000D+00 3.446352960D+04 9.703788500D+00
  1000.000  3000.000 5  0.0  1.0  2.0  3.0  4.0  0.0  0.0  0.0        10698.000
  2.984627599D+00 1.078391826D-02-5.158601830D-06 1.200731137D-09-1.103701620D-13
  0.000000000D+00 0.000000000D+00 0.000000000D+00 3.423078010D+04 7.923373280D+00"""


# Input data
FORMULA = 'CH4'
PAC99_OUTFILE_NAME = os.path.join(os.getcwd(), 'run', FORMULA+'.o97')


def test__polynomial():
    """ reads the polynomial from pac99 and converts it
    """

    # CH4
    print('\nCH4:\n')
    # Open the output from pac99
    with open(PAC99_OUTFILE_NAME, 'r') as pac99_file:
        pac99_str = pac99_file.read()

    # Get the pac99 polynomial
    pac99_poly_str = nasapoly.get_pac99_polynomial(pac99_str)
    print('\nPAC99 Polynomial:')
    print(pac99_poly_str)

    # Convert the pac99 polynomial to chemkin polynomial
    chemkin_poly_str = nasapoly.convert_pac_to_chemkin(pac99_poly_str)
    print('\nCHEMKIN Polynomial:')
    print(chemkin_poly_str)

    # from github
    print('\n\ngithub:\n')

    # Get the pac99 polynomial
    print('\nPAC99 Polynomial:')
    print(GIT_STR)

    # Convert the pac99 polynomial to chemkin polynomial
    chemkin_poly_str = nasapoly.convert_pac_to_chemkin(GIT_STR)
    print('\nCHEMKIN Polynomial:')
    print(chemkin_poly_str)

if __name__ == '__main__':
    test__polynomial()
