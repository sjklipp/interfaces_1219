""" functions operating on the species block string
"""
import autoparse.find as apf


def names(block_str, exclude_names=('CH(6)',)):
    """ species names
    """
    spc_names = apf.split_words(block_str)
    spc_names = tuple(filter(lambda x: x not in exclude_names, spc_names))
    return spc_names
