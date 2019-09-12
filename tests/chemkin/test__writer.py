"""
 tests writers
"""

import chemkin_io
import automol

# def arhenisu_writer_plog

def test__transport_writer():
    """ test chemkin_io.lj.transport
    """

    # ichs = ['InChI=1S/He', 'InChI=1S/C', 'InChI=1S/N2/c1-2',
    #         'InChI=1S/CH/h1H', 'InChI=1S/CH2/h1H2', 'InChI=1S/CH3/h1H3',
    #         'InChI=1S/CH4/h1H4']

    # names = ['HE', 'C', 'N2', 'CH', 'CH2', 'CH3', 'CH4']
    # epsilons = [7.953, 1118.729, 68.001, 395.693,
    #             201.062, 240.615, 302.980, 627.229]
    # sigmas = [2.715, 1.849, 3.610, 2.676, 3.458, 3.535, 3.575, 3.655]
    # dipole_moments = [0.000, 0.000, 0.000, 1.489, 0.593, 0.000, 0.000, 0.000]
    # polarizabilities = [0.204, 1.748, 1.756, 2.271, 2.137, 2.354, 2.454, 3.397]

    geo = automol.inchi.smiles('InChI=1S/N2/c1-2')
    geo = automol.inchi.geometry('InChI=1S/N2/c1-2')

    ichs = ['InChI=1S/N2/c1-2',
            'InChI=1S/CH/h1H', 'InChI=1S/CH2/h1H2', 'InChI=1S/CH3/h1H3',
            'InChI=1S/CH4/h1H4']

    names = ['N2', 'CH', 'CH2', 'CH3', 'CH4']
    epsilons = [68.001, 395.693,
                201.062, 240.615, 302.980, 627.229]
    sigmas = [3.610, 2.676, 3.458, 3.535, 3.575, 3.655]
    dipole_moments = [0.000, 1.489, 0.593, 0.000, 0.000, 0.000]
    polarizabilities = [1.756, 2.271, 2.137, 2.354, 2.454, 3.397]

    transport_str = chemkin_io.transport.lj(
        ichs, names, epsilons, sigmas,
        dipole_moments, polarizabilities, z_rots=[])

    print(transport_str)


if __name__ == '__main__':
    test__transport_writer()
    # test__arrhenius_writer()
