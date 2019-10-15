""" Procedure to do SCT calculations
"""

# Run the irc
# read the geoms, grads, hessians,
# read the scoords, and energies
# for iso, beta: dist_atc_prd=0.0

def reaction_coord_distances():
    """ calculate the distances of the reaction coordinate
    """
    dist_atc_rct = dist(atcent-atreax)
    dist_atc_prd = dist(atcent-atprdx)


