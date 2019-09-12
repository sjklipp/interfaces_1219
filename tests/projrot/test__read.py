"""
 tests reading of projrot output
"""

import projrot_io


def test__reader():
    """ test projrot_io.reader.rpht_output
    """

    # Set file names for various test files
    min_freqs_file = 'min.dat'
    one_imag_freqs_file = 'one_imag.dat'
    two_imag_freqs_file = 'two_imag.dat'

    # Obtain the real and imaginary frequencies from each file
    real_freqs1, imag_freqs1 = projrot_io.reader.rpht_output(
        min_freqs_file)
    real_freqs2, imag_freqs2 = projrot_io.reader.rpht_output(
        one_imag_freqs_file)
    real_freqs3, imag_freqs3 = projrot_io.reader.rpht_output(
        two_imag_freqs_file)

    # Print the frequencies from minimum file
    print('\nFrequencies from min.dat:')
    print('Real Frequencies:')
    for val in real_freqs1:
        print(val)
    print('Imaginary Frequencies:')
    for val in imag_freqs1:
        print(val)

    # Print the frequencies from one imag file
    print('\nFrequencies from one_imag.dat:')
    print('Real Frequencies:')
    for val in real_freqs2:
        print(val)
    print('Imaginary Frequencies:')
    for val in imag_freqs2:
        print(val)

    # Print the frequencies from two imag file
    print('\nFrequencies from two_imag.dat:')
    print('Real Frequencies:')
    for val in real_freqs3:
        print(val)
    print('Imaginary Frequencies:')
    for val in imag_freqs3:
        print(val)


if __name__ == '__main__':
    test__reader()
