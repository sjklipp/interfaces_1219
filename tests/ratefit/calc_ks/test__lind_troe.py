"""
Test the ratefit rate constant calculators for Lindemann and Troe expressions
"""

import numpy as np
import pandas
import ratefit

TEMPS = np.array([
    300., 600., 900., 1200., 1500., 1800., 2100., 2400., 2700., 3000.])
PRESSURES = np.array([0.1, 2.0, 5.0, 10.0])
T_REF = 1.0

A_HIGH, N_HIGH, EA_HIGH = 2.000e+12, 0.900, 4.87490
A_LOW, N_LOW, EA_LOW = 2.490e24, -2.300, 4.87490
TROE_ALPHA, TROE_T3, TROE_T1, TROE_T2 = 6.0e-1, 1.0e3, 7.0, 1.7e3

np.set_printoptions(precision=15)


def _read_csv(filename):
    """ read csv values from file
    """
    csv_file = open(filename, 'r')
    data = pandas.read_csv(csv_file, comment='!', quotechar="'")
    csv_file.close()
    return data


def _calc_highp_ks():
    """ use Arrhenius calculator to get the high-pressure rate constants
    """
    highp_ks = ratefit.fxns.single_arrhenius(
        A_HIGH, N_HIGH, EA_HIGH,
        T_REF, TEMPS)
    return highp_ks


def _calc_lowp_ks():
    """ use Arrhenius calculator to get the low-pressure rate constants
    """
    lowp_ks = ratefit.fxns.single_arrhenius(
        A_LOW, N_LOW, EA_LOW,
        T_REF, TEMPS)
    return lowp_ks


def test__lindemann_and_troe():
    """ test ratefit.fxns.lindemann
    """
    highp_ks = _calc_highp_ks()
    lowp_ks = _calc_lowp_ks()

    lind_ktps = ratefit.fxns.lindemann(
        highp_ks, lowp_ks,
        PRESSURES, TEMPS)

    troe_ktps = ratefit.fxns.troe(
        highp_ks, lowp_ks,
        PRESSURES, TEMPS,
        TROE_ALPHA, TROE_T3, TROE_T1, TROE_T2)
    for k in troe_ktps.values():
        print(k)
    import sys 
    sys.exit()

    data_lind = _read_csv('lindemann.csv')
    data_troe = _read_csv('troe.csv')
    assert np.allclose(highp_ks, np.array(data_lind.kHighP), atol=0.01)
    assert np.allclose(lowp_ks, np.array(data_lind.kLowP), atol=0.01)
    assert np.allclose(lind_ktps[0.1], np.array(data_lind.kLind1), atol=0.01)
    assert np.allclose(lind_ktps[2.0], np.array(data_lind.kLind2), atol=0.01)
    assert np.allclose(lind_ktps[5.0], np.array(data_lind.kLind3), atol=0.01)
    assert np.allclose(lind_ktps[10.0], np.array(data_lind.kLind4), atol=0.01)
    assert np.allclose(troe_ktps[0.1], np.array(data_troe.kTroe1), atol=0.01)
    assert np.allclose(troe_ktps[2.0], np.array(data_troe.kTroe2), atol=0.01)
    assert np.allclose(troe_ktps[5.0], np.array(data_troe.kTroe3), atol=0.01)
    assert np.allclose(troe_ktps[10.0], np.array(data_troe.kTroe4), atol=0.01)


if __name__ == '__main__':
    test__lindemann_and_troe()
