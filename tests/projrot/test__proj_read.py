"""
 tests reading of projrot output
"""

import projrot_io


def test__reader():
    """ test projrot_io.reader.rpht_output
    """

    # Set file names for various test files
    min_freqs_file = 'freq_files/min.dat'
    one_imag_freqs_file = 'freq_files/one_imag.dat'
    two_imag_freqs_file = 'freq_files/two_imag.dat'

    # Obtain the real and imaginary frequencies from each file
    real_freqs1, imag_freqs1 = projrot_io.reader.rpht_output(
        min_freqs_file)
    real_freqs2, imag_freqs2 = projrot_io.reader.rpht_output(
        one_imag_freqs_file)
    real_freqs3, imag_freqs3 = projrot_io.reader.rpht_output(
        two_imag_freqs_file)

    assert real_freqs1 == [6000.0, 5000.0, 4000.0, 3000.0, 2000.0, 1000.0]
    assert not imag_freqs1
    assert real_freqs2 == [6000.0, 5000.0, 4000.0, 3000.0, 2000.0]
    assert imag_freqs2 == [1111.11]
    assert real_freqs3 == [6000.0, 5000.0, 4000.0, 3000.0]
    assert imag_freqs3 == [2222.22, 1111.11]


if __name__ == '__main__':
    test__reader()
