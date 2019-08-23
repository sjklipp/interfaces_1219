"""
test the Arrhenius fitting codes
"""

import numpy as np
import arrfit

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

TEMPS = []
RATE_CONSTANTS = []
for pair in PAIRS:
    TEMPS.append(pair[0])
    RATE_CONSTANTS.append(pair[1])

T_REF = 1.0
R = 8.314


def test__fit():
    """ fit test
    """

    # Run a single Arrhenius fit
    fit_params1, fit_range1 = arrfit.fit.single_arrhenius_fit(
        TEMPS, RATE_CONSTANTS, T_REF)
    #    tmin=600, tmax=max(TEMPS))
    print('\nSingle Arrhenius Fit:')
    print('A =', fit_params1[0])
    print('n =', fit_params1[1])
    print('Ea =', fit_params1[2])
    print('Fit Range =', fit_range1)

    # Calculate fitted rate constants using the fitted parameters
    fit_ks1 = arrfit.fit.single_arrhenius(
        fit_params1[0], fit_params1[1], fit_params1[2],
        T_REF, np.array(TEMPS[3:], dtype=np.float64))

    # Print the fitted rate constants and errors
    for i, _ in enumerate(TEMPS[3:]):
        print('{0}  {1}  {2}'.format(
            TEMPS[i], RATE_CONSTANTS[i], fit_ks1[i]))

    # Calculate the sum-of-square erros and mean-average-errors
    sse, mean_avg_err, max_avg_err = arrfit.fit.calc_sse_and_mae(
        np.array(RATE_CONSTANTS[3:], dtype=np.float64),
        np.array(fit_ks1, dtype=np.float64))
    print('SSE =', sse)
    print('Mean Avg. Err = ', mean_avg_err)
    print('Max Avg. Err = ', max_avg_err)


    # # Run a double Arrhenius fit using SJK code
    # fit_params3, fit_range3 = arrfit.fit.double_arrhenius_fit_dsarrfit(
    #      TEMPS, RATE_CONSTANTS, T_REF):
    #  print('\nDouble Arrhenius Fit (SJK arrfit):'
    #  print('A1 =', fit_params3[0])
    #  print('n1 =', fit_params3[1])
    #  print('Ea1 =', fit_params3[2])
    #  print('A2 =', fit_params3[3])
    #  print('n2 =', fit_params3[4])
    #  print('Ea2 =', fit_params3[5])
    #  print('Fit Range =', fit_range3)

    # # Run a double Arrhenius fit using scipy lsq fit code
    # fit_params4, fit_range4 = arrfit.fit.double_arrhenius_fit_scipy(
    #      TEMPS, RATE_CONSTANTS, T_REF):
    #  print('\nDouble Arrhenius Fit (SJK arrfit):'
    #  print('A1 =', fit_params4[0])
    #  print('n1 =', fit_params4[1])
    #  print('Ea1 =', fit_params4[2])
    #  print('A2 =', fit_params4[3])
    #  print('n2 =', fit_params4[4])
    #  print('Ea2 =', fit_params4[5])
    #  print('Fit Range =', fit_range4)


test__fit()
