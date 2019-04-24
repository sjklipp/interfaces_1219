RRHO
## Core Section
${core}
## Frequencies Section
  Frequencies[1/cm]         ${nfreqs}
${freqs}
## Zero Energy Section
  ZeroEnergy[kcal/mol]      ${zero_energy}
## Electronic Levels Section
  ElectronicLevels[1/cm]    ${nlevels}
${levels}
## Hindered Rotor Section
% if hind_rot != '':
${hind_rot}
% endif 
## Anharmonicity Section
% if anharm != '':
${anharm}
% endif 
## Rovibrational Coupling Section
% if rovib_coup != '':
  RovibrationalCouplings[1/cm]
${rovib_coup}
% endif
## Rotational Distortion Section
% if rot_dist != '':
  RotationalDistortion[1/cm]
${rot_dist}
% endif
## Tunnel Section
% if tunnel != '':
${tunnel}
% endif
End
