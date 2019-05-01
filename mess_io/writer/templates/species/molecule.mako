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
${hind_rot}\
% endif 
## Anharmonicity Section
% if anharm != '':
  Anharmonicities[1/cm]
${anharm}\
% endif 
## Rovibrational Coupling Section
% if rovib_coups != '':
  RovibrationalCouplings[1/cm]
${rovib_coups}\
% endif
## Rotational Distortion Section
% if rot_dists != '':
  RotationalDistortion[1/cm]
${rot_dists}\
% endif
## Tunnel Section
% if tunnel != '':
${tunnel}\
% endif
