"""
Writes strings containing the rate parameters
"""

def plog(reaction, rate_params_dct, err_dct):
    """ Write the string containing the fitting parameters
        formatted for CHEMKIN input files
    """

    # High-pressure rates required?
    assert 'high' in rate_params_dct.keys()

    # Find number of fitting parameters
    nparams = len(rate_params_dct['high'])

    # Make sure that you have 3 params (single fxn) or 6 params (double fxn)
    assert nparams in (3, 6)
    assert all(len(params) == nparams for params in rate_params_dct.values())

    # Obtain a list of the pressures and sort from low to high pressure
    pressures = [pressure for pressure in rate_params_dct.keys()
                 if pressure != 'high']
    pressures.sort()

    # Next is the reaction string and high-pressure params
    # Loop will build second ('DUPLICATE') section if double fit performed
    p_str = ''
    for i in range(nparams // 3):
        if i == 1:
            p_str += 'DUPLICATE\n'

        # Build the initial string with the reaction and high-pressure params
        high_a, high_n, high_ea = rate_params_dct['high'][3*i:3*i+3]
        p_str += '{0:<32s}{1:>10.4E}{2:>9.4f}{3:9.1f} /\n'.format(
            reaction, high_a, high_n, high_ea)

        # Build the PLOG string for each pressure
        for pressure in pressures:
            pdep_a, pdep_n, pdep_ea = rate_params_dct[pressure][3*i:3*i+3]
            p_str += '{0:>18s} /{1:>10.1f}  '.format(
                'PLOG', float(pressure))
            p_str += '{0:>10.4E}{1:>9.4f}{2:9.1f} /\n'.format(
                pdep_a, pdep_n, pdep_ea)

    # writing errors
    for key, val in err_dct.items():
        err_str = '{0:12s} {1:>6.3f}%, {2:8s} {3:>6.3f}%'.format(
            'MeanAbsErr =', val[0], 'MaxErr =', val[1])
        p_str += '! {0:<6s}: {1}\n'.format(str(key), err_str)

    return p_str
