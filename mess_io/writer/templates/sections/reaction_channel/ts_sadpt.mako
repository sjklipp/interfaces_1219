Barrier ${ts_label} ${reac_label} ${prod_label}
${ts_data}\
## Zero Energy Section
    ZeroEnergy[kcal/mol]      ${zero_energy}
## Tunnel Section
% if tunnel != '':
${tunnel}\
% endif
