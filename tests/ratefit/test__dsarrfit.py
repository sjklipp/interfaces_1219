"""
Test an Arrhenius fit of T, k(T,P) rates to both a
single and double Arrhenius function
where the fits are performed using the SJK dsarrfit code
"""

import ratefit


# Obtain list of temperatures and rate constants from initial pair list
PAIRS = [
    ['300', '***'],
    ['400', '***'],
    ['500', '-30.0'],
    ['600', '93.5678'],
    ['700', '2287.02'],
    ['800', '27103.4'],
    ['900', '193348'],
    ['1000', '955781'],
    ['1100', '3.60E+06'],
    ['1200', '1.10E+07'],
    ['1300', '2.85E+07'],
    ['1400', '6.51E+07'],
    ['1500', '1.34E+08'],
    ['1600', '2.52E+08'],
    ['1700', '4.43E+08'],
    ['1800', '7.34E+08'],
    ['1900', '1.16E+09'],
    ['2000', '1.74E+09'],
    ['2100', '2.53E+09'],
    ['2200', '3.56E+09'],
    ['2300', '4.86E+09'],
    ['2400', '6.48E+09'],
    ['2500', '8.45E+09'],
    ['2600', '1.08E+10'],
    ['2700', '1.36E+10'],
    ['2800', '1.68E+10'],
    ['2900', '2.06E+10'],
    ['3000', '2.48E+10']
]
TEMPS = [pair[0] for pair in PAIRS]
RATE_CONSTANTS = [pair[1] for pair in PAIRS]

# Set the T0 value in the (T/T0)^n term in the Arrhenius expr.
T_REF = 1.0

# Set the path to the dsarrfit executable
EXE_PATH = '../../ratefit/fit_arrhenius/external/dsarrfit'


def test__fit():
    """ fit test ratefit.fit_arrhenius.dsarrfit_io
    """

    # Print header
    print('\n\ndsarrfit Results:')

    # Filter the temperatures and rate constants to get valid values
    # k > 0 and k != *** and tmin <= T <= tmax
    tmin = 600
    tmax = 3000
    temps, calc_ks = ratefit.fit_arrhenius.util.get_valid_tk(
        TEMPS, RATE_CONSTANTS,
        tmin=None, tmax=None)
    print('Fit Range =', [tmin, tmax])

    # Write the input file for the ratefit code
    ratefit_inp_str = ratefit.fit_arrhenius.dsarrfit_io.write_input(
        temps, calc_ks)
    with open('arrfit.dat', 'w') as ratefit_infile:
        ratefit_infile.write(ratefit_inp_str)

    # Run the ratefit program
    ratefit.fit_arrhenius.dsarrfit_io.run_dsarrfit(EXE_PATH)

    # Read the output of the single and double fit
    with open('arrfit.out', 'r') as s_ratefit_outfile:
        sfit_out_str = s_ratefit_outfile.read()
    with open('darrfit.out', 'r') as d_ratefit_outfile:
        dfit_out_str = d_ratefit_outfile.read()

    # Parse the ratefit files for the Arrhenius fit parameters
    sfit_params = ratefit.fit_arrhenius.dsarrfit_io.read_params(
        sfit_out_str, 'single', 1.00)
    dfit_params = ratefit.fit_arrhenius.dsarrfit_io.read_params(
        dfit_out_str, 'double', 1.00)

    # Print the single Arrhenius fit parameters
    print('\nSingle Arrhenius Fit Results:')
    print('Fit Parameters:')
    print('A =', sfit_params[0][0])
    print('n =', sfit_params[0][1])
    print('Ea =', sfit_params[0][2])

    # Calculate fitted rate constants using the fitted parameters
    fit_ks1 = ratefit.fxns.single_arrhenius(
        sfit_params[0][0], sfit_params[0][1], sfit_params[0][2],
        T_REF, temps)

    # Print the fitted rate constants and errors
    print('\nComparison of Calculated vs. Fitted Rate Constants:')
    print('Temp (K)  Calc ks      Fit ks')
    for i, _ in enumerate(temps):
        print('{0:6.1f}    {1:1.5E}  {2:1.5E}'.format(
            temps[i], calc_ks[i], fit_ks1[i]))

    # Calculate the sum-of-square errors and mean-average-errors
    sse1, mean_err1, max_err1 = ratefit.err.calc_sse_and_mae(
        calc_ks, fit_ks1)
    print('\nSSE =', sse1)
    print('Mean Avg. Err = ', mean_err1)
    print('Max Avg. Err = ', max_err1)

    # Print the double Arrhenius fit parameters
    print('\n\nDouble Arrhenius Fit Results:')
    print('Fit Parameters:')
    print('A1 =', dfit_params[0][0])
    print('n1 =', dfit_params[0][1])
    print('Ea1 =', dfit_params[0][2])
    print('A2 =', dfit_params[1][0])
    print('n2 =', dfit_params[1][1])
    print('Ea2 =', dfit_params[1][2])

    # Calculate fitted rate constants using the fitted parameters
    fit_ks2 = ratefit.fxns.double_arrhenius(
        dfit_params[0][0], dfit_params[0][1], dfit_params[0][2],
        dfit_params[1][0], dfit_params[1][1], dfit_params[1][2],
        T_REF, temps)

    # Print the fitted rate constants and errors
    print('\nComparison of Calculated vs. Fitted Rate Constants:')
    print('Temp (K)  Calc ks      Fit ks')
    for i, _ in enumerate(temps):
        print('{0:6.1f}    {1:1.5E}  {2:1.5E}'.format(
            temps[i], calc_ks[i], fit_ks2[i]))

    # Calculate the sum-of-square errors and mean-average-errors
    sse2, mean_err2, max_err2 = ratefit.err.calc_sse_and_mae(
        calc_ks, fit_ks2)
    print('\nSSE =', sse2)
    print('Mean Avg. Err = ', mean_err2)
    print('Max Avg. Err = ', max_err2)


if __name__ == '__main__':
    test__fit()
