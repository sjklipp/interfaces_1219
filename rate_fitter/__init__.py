"""
Module deal with rate constant functional forms; Either
 (1) fits a set of rate constants [k(T, P)] to
     various functional forms
 (2) calculates k(T, P) using a functional form given the
     fitting parameters are known
"""

from arrfit.fit import get_valid_temps_rate_constants
from arrfit.fit import single_arrhenius_fit
from arrfit.dsarrfit_io import write_arrfit_inp
from arrfit.dsarrfit_io import run_dsarrfit
from arrfit.dsarrfit_io import parse_dsarrfit


__all__ = [
    'get_valid_temps_rate_constants',
    'single_arrhenius_fit',
    'write_arrfit_inp',
    'run_dsarrfit',
    'parse_dsarrfit'
]
