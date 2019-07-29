""""
Run the make poly
"""

import os

FORMULA = 'CH4'
DELTAH = 15.1
ENTHALPYT = 0.0
BREAKT = 1000.0
THERMP_FILE_NAME = 'thermp.dat'
PF_FILE_NAME = 'pf.dat'
WORK_PATH = os.getcwd()

# Write thermp input file
write_thermp_input(FORMULA, DELTAH,
                   enthalpyT=ENTHALPYT, breakT=BREAKT,
                   thermp_file_name=THERMP_FILE_NAME)

# Run thermp
run_thermp(WORK_PATH,
           thermp_file_name=THERMP_FILE_NAME
           pf_file_name=PF_FILE_NAME)

# Run pac99
run_pac99(WORK_PATH, FORMULA)

# Read the data from thermp
with open() as thermp_file:
    thermp_str = thermp_file.read()
dHf298 = read_thermp_dHf(thermp_str)

# Read the data from pac99
with open() as pac99_file:
    pac99_str = pac99_file.read()





