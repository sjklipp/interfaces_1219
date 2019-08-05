Species
${formula}
  MonteCarlo
    MoleculeSpecification       ${natoms}
${atom_list}
${flux_mode_str}\
    DataFile                    ${data_file_name}
    GroundEnergy[kcal/mol]      ${ground_energy}
    ReferenceEnergy[kcal/mol]   ${reference_energy}
  End
