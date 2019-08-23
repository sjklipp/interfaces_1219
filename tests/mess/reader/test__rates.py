"""
 tests rates reader
"""

import mess_io.reader


def test__read_rates():
    """ Reads the rates from the output of a file.
    """

    # Set the reactant and product
    reactant = 'REACS'
    product = 'WR'

    # Read the MESS output file into a string
    with open('rate.out', 'r') as mess_file:
        output_string = mess_file.read()

    # Read the temperatures from the output
    temps, tunit = mess_io.reader.rates.get_temperatures(output_string)
    print('\nTemperatures:')
    print(temps)
    print(tunit)

    # Read the pressures from the output
    pressures, punit = mess_io.reader.rates.get_pressures(output_string)
    print('\nPressures:')
    print(pressures)
    print(punit)

    # Read the high-pressure rate constants
    highp_rates = mess_io.reader.read_highp_ks(
        output_string, reactant, product)
    print('\nHigh-Pressure Rate-Constants:')
    print(highp_rates)

    # Read pressure-dependent rate constants at two pressures
    p1_rates = mess_io.reader.read_pdep_ks(
        output_string, reactant, product, pressures[3], punit)
    print('\n{} Rate-Constants:'.format(pressures[3]))
    print(p1_rates)

    p2_rates = mess_io.reader.read_pdep_ks(
        output_string, reactant, product, pressures[4], punit)
    print('\n{} Rate-Constants:'.format(pressures[4]))
    print(p2_rates)


if __name__ == '__main__':
    test__read_rates()
