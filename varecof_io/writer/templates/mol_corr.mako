
      real*8 function ${species_name}_corr(natoms,x,rparm,iparm)
c
c     natoms = number of atoms (not used)
c     x = cartesian coordinates
c     iparm(1) = specifiies which potential correction to use
c     rparm (not used)
c     rinp = ${asym}-${bsym} distance
% if corr_labels_str != '':
${corr_labels_str}
% endif
c
      implicit real*8 (a-h,o-z)
      dimension x(3,1)
      dimension iparm({$npot})
      dimension rinp(${npot_terms})
      dimension ${dv_defs}
      dimension dv20(${npot_terms})
      data dvp1,dvpn / 1.0d40,1.0d40 /
${rvals}
${dv_vals}

## Append further info that is needed for atom labeling
      data nrin / ${npot_terms} /'
      data r${asym}${bsym}min,r${asym}${bsym}max / ${rmin},${rmax} /

      ipot = iparm(1)

      na = ${aidx}
      nb = ${bidx}

      r${asym}{bsym} = dsqrt( (x(1,nb)-x(1,na))**2 + (x(2,nb)-x(2,na))**2
     x             (x(3,nb)-x(3,na))**2 )')

      r${asym}{bsym} = r${asym}${bsym}*0.52917


## Append the lines for setting all the values
      delmlt = 1.0d0
      if(r${asym}${bsym}.le.r${asym}${bsym}min) r${asym}${bsym} = r${asym}${bsym}min
      if(r${asym}${bsym}.ge.r${asym}${bsym}max) then
        delmlt = exp(-2.d0*(r${asym}${bsym}-r${asym}${bsym}max))
        r${asym}${bsym}=r${asym}${bsym}max'.format(asym,bsym))
      endif

## Append the lines for declaring the spline functions
${spline_str}

      ${species_name}_corr = ${species_name}_corr*delmlt/627.5095

      return
      end

