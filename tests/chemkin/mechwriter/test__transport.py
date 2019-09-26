"""
 tests writers
"""

import chemkin_io


def test__transport_writer():
    """ test chemkin_io.mechwriter.transport.lj
    """

    geoms = [
        (('He', (0.0, 0.0, 0.0)),),
        (('N', (1.3794996697128696, 0.0, 0.0)),
         ('N', (-1.3794996697128696, 0.0, 0.0))),
        (('C', (0.4050179415, 1.1488934126, 0.0065627641)),
         ('H', (-0.2938067567, -0.2796488451, 1.3726046096)),
         ('H', (-0.2938162708, 0.6800975091, -1.9135987940))),
        (('C', (-1.877887794885e-07, -2.586267345369e-07, -2.8105369855e-06)),
         ('H', (1.9090421536565017, -0.5990348467223128, -0.426257958513)),
         ('H', (-1.4763355117155859, -1.3930299763066363, 0.254723919477)),
         ('H', (-0.4327064541521221, 1.9920650816556826, 0.171536849583))),
    ]
    names = ['HE', 'N2', 'CH2', 'CH4']
    epsilons = [7.953, 68.001, 240.615, 302.980]
    sigmas = [2.715, 3.610, 3.458, 3.575]
    dipole_moments = [0.000, 0.000, 0.593, 0.000]
    polarizabilities = [0.204, 1.756, 2.137, 2.454]

    transport_str = chemkin_io.mechwriter.transport.lennard_jones(
        names, geoms, epsilons, sigmas,
        dipole_moments, polarizabilities, z_rots=None)
    print('\ntransport_str')
    print(transport_str)


if __name__ == '__main__':
    test__transport_writer()
