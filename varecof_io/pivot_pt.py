"""
Sets up pivot point locations for a species
"""

import numpy as np
from qcelemental import constants as qcc

# Conversion factors
DEG2RAD = qcc.conversion_factor('degree', 'radian')
RAD2DEG = qcc.conversion_factor('radian', 'degree')



def calc_pivot_point_internal(xyzP, xyz1):
    """ 
    """

    dxyzP = xyzP - xyz1

    # rho
    rho = np.linalg.norm(dxyzP)

    # theta
    if abs(dxyzP[0]) < 0.01 and abs(dxyzP[1]) < 0.01:
        theta = 0.01
    elif abs(dxyzP[0]) < 0.01 and abs(dxyzP[1]) < 0.01:
        theta = (np.pi / 2.0) * RAD2DEG
    else:
        theta = np.atan2(dxyzP[1], dxyzP[0]) * (180.0 / np.pi)

    # phi
    if abs(dxyzP[2]) < 0.01:
        phi = (np.pi / 2.0) * RAD2DEG
    else:
        phi = ((atan2(np.sqrt(dxyzP[0]**2 + dxyzP[1]**2), dxyzP[2])) * 
               (180.0 / np.pi)

    return rho, theta, phi
    
    

def calc_pivot_point_xyz(xyz1, xyz2, xyz3, in_dist, in_angle, in_dihed):
    """ geometric approach for calculating the xyz coordinates of atom A
        when the xyz coordinates of the A B and C are known and
        the position is defined w/r to A B C with internal coordinates
    """

    # Build initial coordinates
    xyz1 = np.array(xyz1, dtype=float)
    xyz2 = np.array(xyz2, dtype=float)
    xyz3 = np.array(xyz3, dtype=float)
    in_angle *= DEG2RAD
    in_dihed *= DEG2RAD

    # Get the A, B, C, and P points in the rt system
    xyzP_rt = np.array([in_dist * np.sin(in_angle) * np.cos(in_dihed),
                        in_dist * np.cos(in_angle),
                        -(in_dist * np.sin(in_angle) * np.sin(in_dihed))])
    xyz1_rt = np.array([0.0, 0.0, 0.0])

    dist12 = np.linalg.norm(xyz1 - xyz2)
    dist13 = np.linalg.norm(xyz1 - xyz3)
    dist23 = np.linalg.norm(xyz2 - xyz3)
    val = ((dist12**2 + dist13**2 - dist23**2) / 2 / dist12)
    xyz2_rt = np.array([0.0, dist12, 0.0])
    xyz3_rt = np.array([val, np.sqrt(dist13**2 - val**2), 0.0])

    # calculate the check1?
    xyzP_rt_tmp = np.array([xyzP_rt[0], xyzP_rt[1], -xyzP_rt[2]])
    check1 = np.linalg.norm(xyzP_rt_tmp - xyz3_rt)
    
    # translate original frame of ref coors so that xyz1 is at (0, 0, 0)
    xyz1_t = np.array([0.0, 0.0, 0.0])
    xyz2_t = xyz2 - xyz1
    xyz3_t = xyz3 - xyz1

    # rotation matrix to rotate back to the original ref system
    r12 = (xyz2[0] - xyz1[0]) / xyz2_rt[1]
    r22 = (xyz2[1] - xyz1[1]) / xyz2_rt[1]
    r32 = (xyz2[2] - xyz1[2]) / xyz2_rt[1]

    r11 = (xyz3[0] - xyz1[0] - xyz3_rt[1]*r12) / xyz3_rt[0]
    r21 = (xyz3[1] - xyz1[1] - xyz3_rt[1]*r22) / xyz3_rt[0]
    r31 = (xyz3[2] - xyz1[2] - xyz3_rt[1]*r32) / xyz3_rt[0]

    anum_aconst = (xyz2_t[1] - xyz3_t[1]) / (xyz3_t[0]*xyz2_t[0])
    den_aconst = (xyz2_t[2] - xyz3_t[2]) / (xyz3_t[0]*xyz2_t[0])

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
        aconst = (((xyz2_t[2] - xyz3_t[2]) / (xyz3_t[0] - xyz2_t[0])) /
                  ((xyz2_t[3] - xyz3_t[3]) / (xyz3_t[0] - xyz2_t[0])))

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

    # now rotate and translate back
    # here I check  the (001) vector direction to decide whether to take t       he positive of negative results of the 
    # square root taken above
    xap = xyz1[0] + (
        r11 * xyzP_rt[0]) + (r12 * xyzP_rt[1]) + (r13 * xyzP_rt[2])
    yap = xyz1[1] + (
        r21 * xyzP_rt[0]) + (r22 * xyzP_rt[1]) + (r33 * xyzP_rt[2])
    zap = xyz1[2] + (
        r31 * xyzP_rt[0]) + (r32 * xyzP_rt[1]) + (r33 * xyzP_rt[2])

    xan = xyz1[0] + (
        r11 * xyzP_rt[0]) + (r12 * xyzP_rt[1]) + (r13n * xyzP_rt[2])
    yan = xyz1[1] + (
        r21 * xyzP_rt[0]) + (r22 * xyzP_rt[1]) + (r33n * xyzP_rt[2])
    zan = xyz1[2] + (
        r31 * xyzP_rt[0]) + (r32 * xyzP_rt[1]) + (r33n * xyzP_rt[2])

    bvec = np.array([(xyz1[0] - xyz2[1]), (xyz1[1] - xyz2[1]), (xyz1[2] - xyz2[2])])
    cvec = np.array([(xyz2[0] - xyz3[1]), (xyz2[1] - xyz3[1]), (xyz2[2] - xyz3[2])])
    
    vec1 = (bvec[1] * cvec[2]) - (bvec[2] - cvec[1]) 
    vec2 = (bvec[2] * cvec[0]) - (bvec[0] - cvec[2]) 
    vec3 = (bvec[0] * cvec[1]) - (bvec[1] - cvec[0]) 
   
    if abs(xyz4_t[0] > 1.0e-5):
        checkv = vec1 / xyz4_t[0]
    elif abs(xyz4_t[1] > 1.0e-5):
        checkv = vec2 / xyz4_t[1]
    else:
        checkv = vec3 / xyz4_t[2]
 
    if checkv >= 0.0:
        xyzP = np.array([xap, yap, zap])
    else:
        xyzP = np.array([xan, yan, zan])
        
    return xyzP[0], xyzP[1], xyzP[2]


if __name__ == '__main__':
    xy1z
    xp, yp, zp = calc_pivot_point_xyz(xyz1, xyz2, xyz3, in_dist, in_angle, in_dihed)

