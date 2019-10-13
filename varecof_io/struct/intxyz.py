"""
Sets up pivot point locations for a species
"""

import numpy as np
from qcelemental import constants as qcc

# Conversion factors
DEG2RAD = qcc.conversion_factor('degree', 'radian')
RAD2DEG = qcc.conversion_factor('radian', 'degree')


def find_xyzp(xyz1, xyz2, xyz3, pdist, pangle, pdihed):
    """ geometric approach for calculating the xyz coordinates of atom A
        when the xyz coordinates of the A B and C are known and
        the position is defined w/r to A B C with internal coordinates
    """
    # Set to numpy arrays
    xyz1 = np.array(xyz1)
    xyz2 = np.array(xyz2)
    xyz3 = np.array(xyz3)

    # Build initial coordinates
    pangle *= DEG2RAD
    pdihed *= DEG2RAD

    # if something:
    #     xyzp = _two_point(xyz1, xyz2, xyz3, pdist, pangle, pdihed)
    # else:
    #     xyzp = _three_point(xyz1, xyz2, pdist, pangle)
    xyzp = _three_point(xyz1, xyz2, xyz3, pdist, pangle, pdihed)

    return xyzp[0], xyzp[1], xyzp[2]


# def _two_point(xyz1, xyz2, pdist, pangle):
#     """ finds point P when based on two xyz points
#     """
#     return xyzp[0], xyzp[1], xyzp[2]


def _three_point(xyz1, xyz2, xyz3, pdist, pangle, pdihed):
    """ finds point P when based on three xyz points
    """

    # Set the coordinates of Point P in the RT system
    xyzp_rt = np.array([pdist * np.sin(pangle) * np.cos(pdihed),
                        pdist * np.cos(pangle),
                        -(pdist * np.sin(pangle) * np.sin(pdihed))])

    # Set the coordinates of the Point 2 and 3 in the RT system
    dist12 = np.linalg.norm(xyz1 - xyz2)
    dist13 = np.linalg.norm(xyz1 - xyz3)
    dist23 = np.linalg.norm(xyz2 - xyz3)
    xyz2_rt = np.array([0.0, dist12, 0.0])

    val = ((dist12**2 + dist13**2 - dist23**2) / 2.0 / dist12)
    valx3 = np.sqrt(dist13**2 - val**2)
    valy3 = ((dist12**2 + dist13**2 - dist23**2) / 2.0 / dist12)
    xyz3_rt = np.array([valx3, valy3, 0.0])

    # Translate original frame of ref coors so that xyz1 is at (0, 0, 0)
    xyz2_t = xyz2 - xyz1
    xyz3_t = xyz3 - xyz1

    # Rotation matrix to rotate back to the original ref system
    r12 = (xyz2[0] - xyz1[0]) / xyz2_rt[1]
    r22 = (xyz2[1] - xyz1[1]) / xyz2_rt[1]
    r32 = (xyz2[2] - xyz1[2]) / xyz2_rt[1]

    r11 = (xyz3[0] - xyz1[0] - xyz3_rt[1]*r12) / xyz3_rt[0]
    r21 = (xyz3[1] - xyz1[1] - xyz3_rt[1]*r22) / xyz3_rt[0]
    r31 = (xyz3[2] - xyz1[2] - xyz3_rt[1]*r32) / xyz3_rt[0]

    anum_aconst = xyz2_t[1] - (xyz3_t[1] / xyz3_t[0]) * xyz2_t[0]
    den_aconst = xyz2_t[2] - (xyz3_t[2] / xyz3_t[0]) * xyz2_t[0]

    if abs(anum_aconst) < 1.0e-6 and abs(den_aconst) < 1.0e-6:
        if anum_aconst < 0.0:
            aconst = -1.0e20
        else:
            aconst = 1.0e20
    elif abs(den_aconst) < 1.0e-6:
        if anum_aconst < 0.0:
            aconst = -1.0e20
        else:
            aconst = 1.0e20
    else:
        anum = xyz2_t[1] - (xyz3_t[1] / xyz3_t[0]) * xyz2_t[0]
        aden = xyz2_t[2] - (xyz3_t[2] / xyz3_t[0]) * xyz2_t[0]
        aconst = anum / aden

    den1 = (xyz3_t[1] / xyz3_t[0]) - aconst * (xyz3_t[2] / xyz3_t[0])
    if den1 == 0.0:
        den1 = 1.0e-20
    bconst = 1.0 / den1

    # Set vals for another point
    valx = -(1.0 / np.sqrt(1.0 + (bconst**2) * (1.0 + aconst**2)))
    valy = -(valx * bconst)
    xyz4_t = np.array([valx, valy, -(valy * aconst)])

    r13 = xyz4_t[0]
    r23 = xyz4_t[1]
    r33 = xyz4_t[2]
    r13n = -r13
    r23n = -r23
    r33n = -r33

    # Now rotate and translate back
    # Here I check  the (001) vector direction to decide whether
    # To take the positive of negative results of the
    # Square root taken above
    xap = (xyz1[0] + (r11 * xyzp_rt[0]) +
           (r12 * xyzp_rt[1]) + (r13 * xyzp_rt[2]))
    yap = (xyz1[1] + (r21 * xyzp_rt[0]) +
           (r22 * xyzp_rt[1]) + (r33 * xyzp_rt[2]))
    zap = (xyz1[2] + (r31 * xyzp_rt[0]) +
           (r32 * xyzp_rt[1]) + (r33 * xyzp_rt[2]))

    xan = (xyz1[0] + (r11 * xyzp_rt[0]) +
           (r12 * xyzp_rt[1]) + (r13n * xyzp_rt[2]))
    yan = (xyz1[1] + (r21 * xyzp_rt[0]) +
           (r22 * xyzp_rt[1]) + (r23n * xyzp_rt[2]))
    zan = (xyz1[2] + (r31 * xyzp_rt[0]) +
           (r32 * xyzp_rt[1]) + (r33n * xyzp_rt[2]))

    bvec = xyz1 - xyz2
    cvec = xyz2 - xyz3
    vec1 = (bvec[1] * cvec[2]) - (bvec[2] * cvec[1])
    vec2 = (bvec[2] * cvec[0]) - (bvec[0] * cvec[2])
    vec3 = (bvec[0] * cvec[1]) - (bvec[1] * cvec[0])

    if abs(xyz4_t[0]) > 1.0e-5:
        checkv = vec1 / xyz4_t[0]
    elif abs(xyz4_t[1]) > 1.0e-5:
        checkv = vec2 / xyz4_t[1]
    else:
        checkv = vec3 / xyz4_t[2]

    if checkv >= 0.0:
        xyzp = np.array([xap, yap, zap])
    else:
        xyzp = np.array([xan, yan, zan])

    return xyzp


if __name__ == '__main__':
    XYZ1 = (-4.2740380166, 0.2072902345, -1.7229419340)
    XYZ2 = (-5.7980622207, 0.1857236054, -1.7211160199)
    XYZ3 = (-6.1699775584, -0.1534668354, -2.7105155398)
    IN_DIST = 2.818
    IN_ANGLE = 109.1
    IN_DIHED = -167.3
    XP, YP, ZP = find_xyzp(
        XYZ1, XYZ2, XYZ3, IN_DIST, IN_ANGLE, IN_DIHED)
    print('\n\n\n\nactual:')
    print(-3.3521316720, 0.4963073589, 0.9239111161)
    print('\ncalculated:')
    print(XP, YP, ZP)
