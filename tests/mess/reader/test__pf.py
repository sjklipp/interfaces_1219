"""
 tests pf reader
"""

import mess_io.reader


def test__pf():
    """ Reads the output of a messpf file.
    """

    # Read the MESSPF output and store temps and Q fxns
    temps, logq, dq_dt, dq2_dt2 = mess_io.reader.readpf('pf.out')

    # Print the values from the MESSPF output
    for val1, val2, val3, val4 in zip(temps, logq, dq_dt, dq2_dt2):
        print('{0} {1} {2} {3}'.format(val1, val2, val3, val4))


if __name__ == '__main__':
    test__pf()
