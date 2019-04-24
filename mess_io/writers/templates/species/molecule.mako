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
## Tunnel Section
% if tunnel != '':
${tunnel}
% endif
End
