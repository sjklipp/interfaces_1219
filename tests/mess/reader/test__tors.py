"""
 tests torsional mode reader
"""

import mess_io.reader


MESS_LOG_FILE = 'mess.log'


def test__tors():
    """ Reads the output of a mess file.
    """

    # Get the output string
    with open(MESS_LOG_FILE, 'r') as mess_file:
        output_string = mess_file.read()

    # Read the freqs and zpes
    freqs = mess_io.reader.tors.read_freqs(output_string)
    zpes = mess_io.reader.tors.read_zpes(output_string)

    # Print the freqs and zpes
    for i, (freq, zpe) in enumerate(zip(freqs, zpes)):
        print('{0:4d}{1:10.4f}{2:10.4f}'.format(i+1, freq, zpe))


if __name__ == '__main__':
    test__tors()
