"""
Grabs all of the rate constants out of a MESS output file
"""

import sys


def read_rates(rates_file_name, reaction, temps=[], pressures=[]):
    """ Reads the rates
    """

   


def print_rates(temps, pressures):
    """ Prints the rates
    """

    # Print header for rate constant table to the screen
    print('Rate Constants for the Reaction: {0}\n'.format(rxn))
    print_head = '{0:>5}'.format('T K')
    for i in range(len(pressures)):
        print_head = print_head + '{0:>16}'.format(pressures[i]+' atm')
    print_head = print_head + '{0:>16}'.format('High P')
    print(print_head)
    
    # Print values for rate constant table to the screen
    for i in range(len(temps)):
        print_str = '{0:>5}'.format(temps[i])
        for j in range(len(pressures)+1):
            print_str = print_str + '{0:>16}'.format(rate_constants[j][i])
        print(print_str)

    return None


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


# Read in the line for the rate constants that you want
rxn = sys.argv[1]

# Read in the contents of rate.out
with open('rate.out', 'r') as mess_rate_file:
    RATE_LINES = mess_rate_file.readlines()

# Find parts of file where info is read
READ_TEMP_PRESSURE_BLOCK_START = 0
HIGH_PRESSURE_BLOCK_START = 0
TEMP_PRESSURE_BLOCK_START = 0
for i in range(len(RATE_LINES)):
    if 'High Pressure Rate Coefficients (Temperature-Species Rate Tables):' in RATE_LINES[i]:
        HIGH_PRESSURE_BLOCK_START = i
    if 'Pressure-Species Rate Tables:' in RATE_LINES[i]:
        READ_TEMP_PRESSURE_BLOCK_START = i
    if 'Temperature-Species Rate Tables:' in RATE_LINES[i]:
        TEMP_PRESSURE_BLOCK_START = i
        break

# Get the temperatures and pressures
TEMPS = []
PRESSURES = []
TEMP_UNIT = ''
PRESSURE_UNIT = ''
for i in range(READ_TEMP_PRESSURE_BLOCK_START, len(RATE_LINES)):
    if 'Temperature =' in RATE_LINES[i]:
        tmp = RATE_LINES[i].strip().split()
        if tmp[2] not in TEMPS:
            TEMPS.append(tmp[2])
        else:
            TEMP_UNIT = tmp[3]
            break
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

# Initialize series of lists to store rate constants
rate_constants = []

# Get the rate constants at each listed pressure
for i in range(TEMP_PRESSURE_BLOCK_START, len(RATE_LINES)):
    if 'Temperature-Pressure Rate Tables:' in RATE_LINES[i]:
        break
    elif rxn in RATE_LINES[i]:
        rate_const_block_start = i
        rate_constants.append(grab_rate_constants(rate_const_block_start, rxn))

# Get the high-pressure rate constants
for i in range(HIGH_PRESSURE_BLOCK_START, len(RATE_LINES)):
    if rxn in RATE_LINES[i]:
        rate_const_block_start = i
        rate_constants.append(grab_rate_constants(rate_const_block_start, rxn))
        break
