"""
Functions to fit rate constants to Single or Double Arrhenius Functions
Performs fits either using SciPy or SJK's dsarrfit code
"""

from rate_fitter.arrhenius.fit import to_single_arrhenius
from rate_fitter.arrhenius.fit import to_double_arrhenius


__all__ = [
    'single_arrhenius',
    'double_arrhenius'
]


