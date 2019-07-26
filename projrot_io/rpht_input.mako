## Set the number of atoms 
Number_of_Atoms:                ${natoms}
## Keywords that will generally not change
Act_energy(kcal/mol):           0.0
Initial_Temperature:            200
Temperature_steps:              40
Temperature_increment:          50
Delta_Energy_rea:               0.
Delta_Energy_pro:               0.
Maxstep:                        1
Npointsint:                     5
Maxtdev:                        0.5
Rearrange(1=yes,0=no)           1
SaddlePoint                     1
ds(1=noexp,0=standard)          0
isct_vtst(1=vtst_sct,0=sct)     1
zerocurvature(1)                0
reduced_mass                    1.0
minimum_frequency               50
anim_freq(if_Maxstep=1)         2
onlyrotors(0=yes,1=no)          0
proj_rea_coo(0=yes(def),1=no)   1
## Set if projections to be done in internal or cartesian coordinates
% if coord_proj == 'cartesian':
internalcoord(1=yes)            0
% elif coord_proj == 'internal':
internalcoord(1=yes)            1
% endif
## Define all of the rotors
numrotors                       ${nrotors}
${rotors_str}
## Step Label
Step1
## Geometry Info
geometry
${geom_str}\
## Gradient Info
gradient
${grad_str}\
## Hessian Info
Hessian
${hess_str}\
