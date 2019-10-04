"""
Test the ratefit rate constant calculators
"""

import numpy as np
import pandas
import ratefit

TEMPS = np.array([
    300., 400., 500., 600., 700., 800., 900., 1000.,
    1100., 1200., 1300., 1400., 1500., 1600., 1700., 1800., 1900., 2000.,
    2100., 2200., 2300., 2400., 2500., 2600., 2700., 2800., 2900., 3000.
])
T_REF = 1.0

A1, N1, EA1 = 7.015e+4, 2.053, -0.3557
A2, N2, EA2 = 5.757e+12, -0.664, 0.3318

np.set_printoptions(precision=15)


def _read_csv():
    """ read csv values from arrhenius.csv
    """
    csv_file = open('arrhenius.csv', 'r')
    data = pandas.read_csv(csv_file, comment='!', quotechar="'")
    csv_file.close()
    return data


def test__single_arrhenius():
    """ test ratefit.fxns.single_arrhenius
    """
    calc_ks = ratefit.fxns.single_arrhenius(
        A1, N1, EA1,
        T_REF, TEMPS)
    data = _read_csv()
    assert np.allclose(calc_ks, np.array(data.SingleArr), atol=0.01)


def test__double_arrhenius():
    """ test ratefit.fxns.double_arrhenius
    """
    calc_ks = ratefit.fxns.double_arrhenius(
        A1, N1, EA1,
        A2, N2, EA2,
        T_REF, TEMPS)
    data = _read_csv()
    print(calc_ks)
    print(np.array(data.DoubleArr))
    assert np.allclose(calc_ks, np.array(data.DoubleArr), atol=0.01)


if __name__ == '__main__':
    test__single_arrhenius()
    test__double_arrhenius()
