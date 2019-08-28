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

    # Read the freqs and zpves
    freqs = mess_io.reader.tors.freqs(output_string)
    zpves = mess_io.reader.tors.zpves(output_string)

    # Print the freqs and zpves
    for i, (freq, zpve) in enumerate(zip(freqs, zpves)):
        print('{0:4d}{1:10.4f}{2:10.4f}'.format(i+1, freq, zpve))


if __name__ == '__main__':
    test__tors()
