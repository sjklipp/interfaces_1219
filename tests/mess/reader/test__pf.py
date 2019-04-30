"""
 tests pf reader
"""

import mess_io.reader


def test__pf():
    """ Reads the output of a messpf file.
    """

    # Read the MESSPF output and store temps and Q fxns
    temps, logQ, dQ, dQ2 = mess_io.reader.readpf('pf.out')
   
    # Print the values from the MESSPF output
    for a, b, c, d in zip(temps, logQ, dQ, dQ2):
        print('{0} {1} {2} {3}'.format(a, b, c, d))


if __name__ == '__main__':
    test__pf()
