"""
 Obtain rate constants [k(T, P)] for a given reaction
"""


def read_highp_ks(output_string, reaction, run_temps, run_pressures):
    """ Read the high-pressure rate constants
    """

    for i in range(len(RATE_LINES)):
        if 'High Pressure Rate Coefficients (Temperature-Species Rate Tables):' in RATE_LINES[i]:
            HIGH_PRESSURE_BLOCK_START = i
            break

    # Get the high-pressure rate constants
    for i in range(HIGH_PRESSURE_BLOCK_START, len(RATE_LINES)):
        if rxn in RATE_LINES[i]:
            rate_const_block_start = i
            rate_constants.append(grab_rate_constants(rate_const_block_start, rxn))
            break

    return rate_constants


def read_pdep_ks(output_string, reaction, run_temps, run_pressures):
    """ Read the pressure-dependent rate constants
    """
    
    for i in range(len(RATE_LINES)):
        if 'Temperature-Species Rate Tables:' in RATE_LINES[i]:
            TEMP_PRESSURE_BLOCK_START = i
            break
    
    # Get the rate constants at each listed pressure
    for i in range(TEMP_PRESSURE_BLOCK_START, len(RATE_LINES)):
        if 'Temperature-Pressure Rate Tables:' in RATE_LINES[i]:
            break
        elif rxn in RATE_LINES[i]:
            rate_const_block_start = i
            rate_constants.append(grab_rate_constants(rate_const_block_start, rxn))

    return rate_constants


def grab_rate_constants(block_start, rxn):
    """ Obtain the rate constants for each pressure """

    # Find the column corresponding to the reaction
    rxn_col = 0
    rxn_headers = RATE_LINES[block_start].strip().split()
    for i in range(len(rxn_headers)):
        if rxn == rxn_headers[i]:
            rxn_col = i
            break

    # Parse the following lines and store the constants in a list
    rate_constants = []
    for i in range(block_start+1, len(RATE_LINES)):
        if RATE_LINES[i].strip() == '':
            break
        else:
            rate_constants.append(RATE_LINES[i].strip().split()[rxn_col])

    return rate_constants


def get_temperatures(output_string):
    """ Determine the temperatures run in the reaction
    """

    RATE_LINES = output_string.splitlines()
    for i in range(len(RATE_LINES)):
        READ_TEMP_PRESSURE_BLOCK_START = i
    TEMPS = []
    for i in range(READ_TEMP_PRESSURE_BLOCK_START, len(RATE_LINES)):
        if 'Temperature =' in RATE_LINES[i]:
            tmp = RATE_LINES[i].strip().split()
            if tmp[2] not in TEMPS:
                TEMPS.append(tmp[2])
            else:
                TEMP_UNIT = tmp[3]
                break

    return TEMPS, TEMP_UNIT


def get_pressures(output_string):
    """ Determine the pressures run in the reaction
    """

    RATE_LINES = output_string.splitlines()
    for i in range(len(RATE_LINES)):
        READ_TEMP_PRESSURE_BLOCK_START = i
    PRESSURES = []
    for i in range(READ_TEMP_PRESSURE_BLOCK_START, len(RATE_LINES)):
        if 'P(' in RATE_LINES[i]:
            PRESSURE_UNIT = RATE_LINES[i].strip().split('(')[1].split(')')[0]
            press_start = i+1
            for j in range(press_start, len(RATE_LINES)):
                if 'O-O' in RATE_LINES[j]:
                    break
                else:
                    tmp = RATE_LINES[j].strip().split()
                    PRESSURES.append(tmp[0])
            break

    return PRESSURES, PRESSURE_UNIT
