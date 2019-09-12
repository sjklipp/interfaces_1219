"""
 writes the string for the chemkin
"""

import automol

# convert wavenumbers to Kelvin (check)
CM2K = 1.438776877


def lj(names, ichs,
       epsilons, sigmas,
       dipole_moments, polarizabilities,
       z_rots=[]):
    """ writes the string for the transport
    """

    # Initialize string with common header
    chemkin_string = """! THEORETICAL TRANSPORT PROPERTIES
!
! (1) Shape Idx, for 0 (atom), 1 (linear molecule), or 2 (nonlinear molecule);
! (2) Epsilon, the Lennard-Jones well depth, in K; and
! (3) Sigma, the Lennard-Jones collision diameter, in Angstrom.
! (4) Mu, total dipole moment, in Debye
! (5) Alpha, mean static polarizability, in Angstrom^3
! (6) Z_rot, Z rotational number
! Species     Shape Idx       Epsilon     Sigma      Mu      Alpha     Z_Rot"""

    # Find the length of the longest name string for formatting
    maxlen = 0
    for name in names:
        maxlen = max(maxlen, len(name))
    maxlen = str(maxlen)

    # get the shape index
    shape_idxs = []
    for ich in ichs:
        geo = automol.inchi.geometry(ich)
        if automol.geom.is_atom(geo):
            shape_idx = 0
        else:
            if automol.geom.is_linear(geo):
                shape_idx = 1
            else:
                shape_idx = 2
        shape_idxs.append(shape_idx)

    # Convert the cm-1 to K
    epsilons = [epsilon * CM2K for epsilon in epsilons]

    # Build the string
    mol_data = zip(names, shape_idxs, epsilons,
                   sigmas, dipole_moments, polarizabilities, z_rots)
    for name, shape, eps, sig, dmom, polr, zrot in mol_data:
        chemkin_string += (
            '{0:<'+maxlen+'}{1:>10.3f}{2:>10.3f}{3:>10.3f}' +
            '{4:>10.3f}{5:>10.3f}{6:>10.3f}').format(
            name, shape, eps, sig, dmom, polr, zrot)

    return chemkin_string
