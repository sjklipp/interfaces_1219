"""
Functions to fit rate constants to Single or Double Arrhenius Functions
Performs fits either using SciPy or SJK's dsarrfit code
"""

from ratefit.fit_arrhenius.fit import single_arrhenius
from ratefit.fit_arrhenius.fit import double_arrhenius_dsarrfit
from ratefit.fit_arrhenius.fit import double_arrhenius_scipy
from ratefit.fit_arrhenius.util import get_valid_tk
from ratefit.fit_arrhenius.dsarrfit_io import write_input
from ratefit.fit_arrhenius.dsarrfit_io import run_dsarrfit
from ratefit.fit_arrhenius.dsarrfit_io import read_params


__all__ = [
    'single_arrhenius',
    'double_arrhenius_dsarrfit',
    'double_arrhenius_scipy',
    'write_input',
    'run_dsarrfit',
    'read_params'
]
