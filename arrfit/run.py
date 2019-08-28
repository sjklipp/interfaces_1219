"""
  Routine to get a set of Arrhenius fits from MESS output
"""

import arrfit
import mess_io


# MESS info
MESS_PATH = 'rate.out'
REACTANT = ''
PRODUCT = ''

# Desired info
FIT_PS = ['0.1', '10.0', 'high']
TMIN = None
TMAX = None

# Dictionaries to store info; indexed by pressure (given in FIT_PS)
K_DCT = {}
FIT_PARAM_DCT = {}
FIT_K_DCT = {}
FIT_ERR_DCT = {}

# Fit info
T_REF = 1.0
FIT_TYPE = 'single'
FIT_METHOD = 'python'

# Read the MESS output file into a string
with open(MESS_PATH, 'r') as mess_file:
    OUTPUT_STRING = mess_file.read()

# Read the temperatures and pressures out of the MESS output
TEMPS, TUNIT = mess_io.reader.rates.get_temperatures(OUTPUT_STRING)
PRESSURES, PUNIT = mess_io.reader.rates.get_pressures(OUTPUT_STRING)

# Obtain the rate constants from the MESS output
for pressure in FIT_PS:

    # Read the rate constants
    if pressure == 'high':
        rate_ks = mess_io.reader.highp_ks(
            OUTPUT_STRING, REACTANT, PRODUCT)
    else:
        rate_ks = mess_io.reader.pdep_ks(
            OUTPUT_STRING, REACTANT, PRODUCT, PRESSURES, PUNIT)

    # Store in a the dictionary
    K_DCT[pressure] = rate_ks

# Filter any temperatures and rate_constants stored in the dictionary
for key, val in K_DCT.items():
    temps, calc_ks = arrfit.fit.get_valid_temps_rate_constants(
        TEMPS, val, tmin=TMIN, tmax=TMAX)
    K_DCT[key] = val

# Calculate the fitting parameters
for key, val in K_DCT.items():

    # Obtain the fitting parameters based on the desired fit and method
    if FIT_TYPE == 'single' and FIT_METHOD == 'python':
        fit_params = arrfit.fit.single_arrhenius_fit(
            TEMPS, val, T_REF)
    else:
        raise NotImplementedError

    # Store the fitting parameters in a dictionary
    FIT_PARAM_DCT[key] = fit_params

# Calculate fitted rate constants using the fitted parameters
for key, val in FIT_PARAM_DCT.items():

    # Calculate fitted rate constants, based on fit type
    if FIT_TYPE == 'single':
        fit_ks = arrfit.fit.single_arrhenius(
            val[0], val[1], val[2],
            T_REF, TEMPS)
    else:
        raise NotImplementedError

    # Store the fitting parameters in a dictionary
    FIT_K_DCT[key] = fit_params

# Calculate the errors for each fit k
for key, val in FIT_K_DCT.items():
    # Calculate the rror between the calc and fit ks
    sse, mean_avg_err, max_avg_err = arrfit.fit.calc_sse_and_mae(
        K_DCT[key], val)

    # Store in a dictionary
    FIT_ERR_DCT[key] = [sse, mean_avg_err, max_avg_err]
