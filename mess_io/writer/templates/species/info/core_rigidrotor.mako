  Geometry[angstrom]        ${natom1}
${geom1}
  Core RigidRotor
    SymmetryFactor          ${sym_factor}
% if interp_emax != '':     
    ZeroPointEnergy[1/cm]             0.0
    InterpolationEnergyMax[kcal/mol]  ${interp_emax}
% endif
  End\
