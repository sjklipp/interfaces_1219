
      real*8 function ${species_name}_corr(natoms,x,rparm,iparm)
c
c     natoms = number of atoms (not used)
c     x = cartesian coordinates
c     rparm (not used)
c     iparm(1) = specifies which potential correction to use
c     rinp = ${asym}-${bsym} distance
% if pot_labels != '':
${pot_labels}
% endif
c
      implicit real*8 (a-h,o-z)
      dimension x(3,1)
      dimension iparm(${npot})
      dimension rinp(${npot_terms})
      dimension ${dv_defs}
      dimension dv20(${npot_terms})
      data dvp1,dvpn / 1.0d40,1.0d40 /
${rvals}
${dv_vals}
## Append further info that is needed for atom labeling
      data nrin / ${npot_terms} /
      data r${asym}${bsym}min,r${asym}${bsym}max / ${rmin},${rmax} /

      ipot = iparm(1)

      ${bond_dist_string}

## Set strings if there are distance comparisons to adjust potential
% for i in range(comp_distance_strings):
comp_distance_strings[i]
% endfor

## Append the lines for setting all the values
      delmlt = 1.0d0
      if(r${asym}${bsym}.le.r${asym}${bsym}min) r${asym}${bsym} = r${asym}${bsym}min
      if(r${asym}${bsym}.ge.r${asym}${bsym}max) then
        delmlt = exp(-2.d0*(r${asym}${bsym}-r${asym}${bsym}max))
        r${asym}${bsym}=r${asym}${bsym}max
      endif

## Append the lines for declaring the spline functions
${spline_strings}

      ${species_name}_corr = ${species_name}_corr*delmlt/627.5095

      return
      end

