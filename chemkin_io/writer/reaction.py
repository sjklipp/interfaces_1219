"""
Writes strings containing the rate parameters
"""


def plog(reaction, rate_params_dct, temp_dct=None, err_dct=None):
    """ Write the string containing the fitting parameters
        formatted for CHEMKIN input files
    """

    # Find nparams and ensure there are correct num in each dct entry
    nparams = len(next(iter(rate_params_dct.values())))
    assert nparams in (3, 6)
    assert all(len(params) == nparams for params in rate_params_dct.values())

    # Obtain a list of the pressures and sort from low to high pressure
    pressures = [pressure for pressure in rate_params_dct.keys()
                 if pressure != 'high']
    pressures.sort()

    # Add fake high pressure parameters if they are not in the dictionary
    if 'high' not in rate_params_dct:
        if nparams == 3:
            rate_params_dct['high'] = [1.00, 0.00, 0.00]
        elif nparams == 6:
            rate_params_dct['high'] = [1.00, 0.00, 0.00, 1.00, 0.00, 0.00]

    # Build the reaction string with high-pressure params and any plog params
    # Loop will build second ('DUPLICATE') section if double fit performed
    p_str = ''
    for i in range(nparams // 3):
        if i == 1:
            p_str += 'DUPLICATE\n'

        # Build the initial string with the reaction and high-pressure params
        high_a, high_n, high_ea = rate_params_dct['high'][3*i:3*i+3]
        p_str += '{0:<32s}{1:>10.3E}{2:>9.3f}{3:9.0f} /\n'.format(
            reaction, high_a, high_n, 1000*high_ea)

        # Build the PLOG string for each pressure, other than the HighP Limit
        for pressure in pressures:
            pdep_a, pdep_n, pdep_ea = rate_params_dct[pressure][3*i:3*i+3]
            p_str += '{0:>18s} /{1:>10.3f}  '.format(
                'PLOG', float(pressure))
            p_str += '{0:>10.3E}{1:>9.3f}{2:9.0f} /\n'.format(
                pdep_a, pdep_n, 1000*pdep_ea)

    # Write string showing the temp fit range and fit errors
    if temp_dct or err_dct:
        p_str += _fit_info_str(pressures, temp_dct, err_dct)

    return p_str


def _fit_info_str(pressures, temp_dct, err_dct):
    """ Write the string detailing the temperatures and errors associated
        with the rate constant fits at each pressure
    """

    # Make temp, err dcts empty if fxn receives None; add 'high' to pressures
    temp_dct = temp_dct if temp_dct else {}
    err_dct = err_dct if err_dct else {}
    if 'high' in temp_dct or 'high' in err_dct:
        pressures = ['high'] + pressures

    # Check the temp and err dcts have same presures as rate_dcts
    if temp_dct:
        assert set(pressures) == set(temp_dct.keys())
    err_dct = err_dct if err_dct else {}
    if err_dct:
        assert set(pressures) == set(err_dct.keys())

    # Write string showing the temp fit range and fit errors
    inf_str = '! Info Regarding Rate Constant Fits\n'
    for pressure in pressures:
        if err_dct:
            [mean_err, max_err] = err_dct[pressure]
            err_str = '{0:12s} {1:>5.1f}%,  {2:8s} {3:>5.1f}%'.format(
                'MeanAbsErr =', mean_err, 'MaxErr =', max_err)
        else:
            err_str = ''
        if temp_dct:
            [min_temp, max_temp] = temp_dct[pressure]
            temp_range_str = '{0:11s} {1:>.0f}-{2:<.0f} K, '.format(
                'TempRange =', min_temp, max_temp)
        else:
            temp_range_str = ''
        # Put together the who info string
        if pressure != 'high':
            pstr = '{0:<10.3f}'.format(pressure)
        else:
            pstr = '{0:<10s}'.format('High')
        inf_str += '! {0}: {1} {2}\n'.format(pstr, temp_range_str, err_str)

    return inf_str
