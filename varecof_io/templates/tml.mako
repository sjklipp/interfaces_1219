memory,${memory},m

GEOMETRY_HERE

basis = ${basis}
${wfn}
if (iterations.ge.0) then
   ${method}
   molpro_energy = energy + ${inf_sep_energy}
else
   molpro_energy = 10.0
endif
---
