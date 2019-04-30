!+++++++++++++++++++++++++++++++++++++++++++++++++++
!  ENERGY TRANSFER SECTION
!+++++++++++++++++++++++++++++++++++++++++++++++++++
Model
  EnergyRelaxation
    Exponential
       Factor[1/cm]                     ${exp_factor}
       Power                            ${exp_power}
       ExponentCutoff                   ${exp_cutoff}
    End
  CollisionFrequency
    LennardJones
       Epsilons[K]                      ${epsilon1}    ${epsilon2}
       Sigmas[angstrom]                 ${sigma1}    ${sigma2}
       Masses[amu]                      ${mass1}    ${mass2}
    End\
