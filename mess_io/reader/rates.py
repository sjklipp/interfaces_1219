"""
 Obtain rate constants [k(T, P)] for a given reaction
"""


def read_highp_ks(output_string, reactant, product):
    """ Read the high-pressure rate constants
    """

    # Build the reaction string found in the MESS output
    reaction = reactant + '->' + product

    # Get the MESS output lines
    mess_lines = output_string.splitlines()

    # Find where the block of text where the high-pressure rates exist
    block_str = ('High Pressure Rate Coefficients ' +
                 '(Temperature-Species Rate Tables):')
    for i, line in enumerate(mess_lines):
        if block_str in line:
            block_start = i
            break

    # Get the high-pressure rate constants
    rate_constants = []
    for i in range(block_start, len(mess_lines)):
        if reaction in mess_lines[i]:
            rate_const_block_start = i
            rate_constants = grab_rate_constants(
                mess_lines, rate_const_block_start, reaction)
            break

    return rate_constants


def read_pdep_ks(output_string, reactant, product, pressure, pressure_unit):
    """ Read the pressure-dependent rate constants
    """

    # Build the reaction string found in the MESS output
    reaction = reactant + '->' + product

    # Get the MESS output lines
    mess_lines = output_string.splitlines()

    # Find where the block of text where the high-pressure rates exist
    block_str = ('Temperature-Species Rate Tables:')
    pressure_str = 'Pressure = ' + pressure + ' ' + pressure_unit

    for i, line in enumerate(mess_lines):
        if block_str in line:
            for j in range(i, len(mess_lines)):
                if 'Temperature-Pressure Rate Tables' in mess_lines[j]:
                    break
                if reaction in mess_lines[j]:
                    if pressure_str in mess_lines[j-2]:
                        rate_const_block_start = j
                        rate_constants = grab_rate_constants(
                            mess_lines, rate_const_block_start, reaction)

    return rate_constants


def grab_rate_constants(mess_lines, block_start, rxn):
    """ Obtain the rate constants for each pressure """

    # Find the column corresponding to the reaction
    rxn_col = 0
    rxn_headers = mess_lines[block_start].strip().split()
    for i, rxn_header in enumerate(rxn_headers):
        if rxn == rxn_header:
            rxn_col = i
            break

    # Parse the following lines and store the constants in a list
    rate_constants = []
    for i in range(block_start+1, len(mess_lines)):
        if mess_lines[i].strip() == '':
            print(mess_lines[i])
            break
        else:
            rate_constants.append(mess_lines[i].strip().split()[rxn_col])

    return rate_constants


def get_temperatures(output_string):
    """ Determine the temperatures run in the reaction
    """

    # Get the MESS output lines
    mess_lines = output_string.splitlines()

    # Find the block of lines where the temperatures can be read
    temp_str = 'Pressure-Species Rate Tables:'
    for i, line in enumerate(mess_lines):
        if temp_str in line:
            block_start = i
            break

    # Read the temperatures
    temps = []
    for i in range(block_start, len(mess_lines)):
        if 'Temperature =' in mess_lines[i]:
            tmp = mess_lines[i].strip().split()
            if tmp[2] not in temps:
                temps.append(tmp[2])
            else:
                temp_unit = tmp[3]
                break

    return temps, temp_unit


def get_pressures(output_string):
    """ Determine the pressures run in the reaction
    """

    # Get the MESS output lines
    mess_lines = output_string.splitlines()

    # Find the block of lines where the pressures can be read
    pressure_str = 'Pressure-Species Rate Tables:'
    for i, line in enumerate(mess_lines):
        if pressure_str in line:
            block_start = i
            break

    # Read the pressures
    pressures = []
    for i in range(block_start, len(mess_lines)):
        if 'P(' in mess_lines[i]:
            pressure_unit = mess_lines[i].strip().split('(')[1].split(')')[0]
            pressure_start = i+1
            for j in range(pressure_start, len(mess_lines)):
                if 'O-O' in mess_lines[j]:
                    break
                else:
                    tmp = mess_lines[j].strip().split()
                    pressures.append(tmp[0])
            break

    return pressures, pressure_unit
