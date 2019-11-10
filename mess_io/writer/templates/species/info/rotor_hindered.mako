Rotor  Hindered
% if geom:
  Geometry[angstrom]   ${natom}
${geom}
% endif
% if use_quantum_weight:
  UseQuantumWeight
% endif
  Group                ${group}
  Axis                 ${axis}
  Symmetry             ${symmetry}
  Potential[kcal/mol]  ${npotential}
${potential} 
End
