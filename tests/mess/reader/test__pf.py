"""
 tests pf reader
"""

import mess_io.reader


def test__pf():
    """ Reads the output of a messpf file.
    """

    # Read the pf output file as a string
    with open('pf.out', 'r') as pf_file:
        pf_out_str = pf_file.read()

    # Read the MESSPF output and store temps and Q fxns
    temps, logq, dq_dt, dq2_dt2 = mess_io.reader.pfs.partition_fxn(pf_out_str)

    # Print the values from the MESSPF output
    for val1, val2, val3, val4 in zip(temps, logq, dq_dt, dq2_dt2):
        print('{0} {1} {2} {3}'.format(val1, val2, val3, val4))


if __name__ == '__main__':
    test__pf()
