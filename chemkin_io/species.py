""" functions operating on the species block string
"""
import autoparse.find as apf


def names(block_str):
    """ species names
    """
    spc_names = apf.split_words(block_str)
    return spc_names
