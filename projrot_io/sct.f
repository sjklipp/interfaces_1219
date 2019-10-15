cc first create the at_geom vector as necessary for freq projection

      do j=1,numpointstot
         if(j.ne.numpointsf+1)then
            do  iatom = 1, natom
               open (unit=99, status='unknown')
               write(99,*)grad(iatom,j)
               rewind(99)
               read(99,*)ind1,ind2
               ind3=0
               rewind(99)
               write(99,*)atgeom_me(iatom,j)
               rewind(99)
               read(99,*)cjunk,coox,cooy,cooz
               rewind(99)
               write(99,1199)ind1,ind2,ind3,coox,cooy,cooz
               rewind(99)
               read (99,'(A100)') atgeom_pr(iatom,j)
               close(99)
            enddo
         else
            open(unit=108,file='geoms/tsgta_l1.xyz',status='unknown')
            read(108,*)
            read(108,*)
            do iatom = 1, natom
               read(108,*)cjunk,coox,cooy,cooz
               open (unit=99, status='unknown')
               write(99,*)grad(iatom,numpointsf)
               rewind(99)
               read(99,*)ind1,ind2
               ind3=0
               close(99)
               open (unit=99, status='unknown')
               rewind(99)
               write(99,1199)ind1,ind2,ind3,coox,cooy,cooz
               rewind(99)
               read (99,'(A100)') atgeom_pr(iatom,numpointsf+1)
               close(99)
            enddo
            close(108)
         endif
      enddo
 1199 format (I3,1x,I3,1x,I3,1x,F9.5,1x,F9.5,1x,F9.5)


cc now write coordinates and energies file
c      
      open (unit=333,file='RPHt_coord_en.dat',status='unknown')

      write(333,*)'Point Coordinate Energy Bond1 Bond2'

      call read_zmat(atomlabel,natom,natomt,intcoor,bislab,ibconn,
     $ iaconn,idconn,bname,anname,dname,atname,idummy,isited,jsited,
     $ ksited,bconnt,aconnt,dconnt)


      do i=1,numpointsf+numpointsb+1
            open (unit=99, status='unknown')
c            write(99,*)'isited is',isited
c            write(99,*)'natomt is',natomt
c            write(99,*)'natom is',natom
            write(99,*)atgeom_me(isited,i)
            write(99,*)atgeom_me(jsited,i)
            if(iadd.eq.1.or.iabs.eq.1)then
               write(99,*)atgeom_me(natom1+1,i)
            else if (iiso.eq.1.or.ibeta.eq.1) then
               write(99,*)'X', 0.0, 0.0, 0.0
            endif
            rewind(99)
            read(99,*)cjunk,atcentx,atcenty,atcentz
            read(99,*)cjunk,atreax,atreay,atreaz
            read(99,*)cjunk,atprodx,atprody,atprodz
            close(99)
            dist_atc_rea=sqrt((atcentx-atreax)**2.
     $       +(atcenty-atreay)**2.
     $       +(atcentz-atreaz)**2.0)
            dist_atc_pro=sqrt((atcentx-atprodx)**2.0
     $       +(atcenty-atprody)**2.0
     $       +(atcentz-atprodz)**2.0)
 

c            write(*,*)'energy = ',rc_ene(i)
         if(iadd.eq.1.or.iabs.eq.1)then
            write(333,1333)i,rc_coord(i),rc_ene(i),dist_atc_rea,
     $                  dist_atc_pro
            else if (iiso.eq.1.or.ibeta.eq.1) then
               write(333,1333)i,rc_coord(i),rc_ene(i),dist_atc_rea,
     $                  0.0
            endif
      enddo
c      do i=1,numpointsf+numpointsb+1
c         write(*,*)'FC = ',force_con(2,i)
c      enddo
      close(333)
 
 1333 format (I3,1X,F8.5,1X,F10.7,1X,F8.5,1X,F8.5)
 1032 format (A120)

cc
cc set to 0 the gradient at the saddle point

      do j=1,natom
         gradts(j)=grad(j,numpointsf+1)
         grad(j,numpointsf+1)=' 1 1 0. 0. 0.'
      enddo


cc at this point, all vectors are filled in 
cc now we can compute the projected frequencies

      open (unit=16,file='freqout.dat'
     $         ,status='unknown')
      open (unit=17,file='fresub.dat'
     $         ,status='unknown')
      open (unit=18,file='freqRTout.dat'
     $         ,status='unknown')

cc write input for multidimensional tunneling calculations

      open (unit=333,file='RPHt_all_data.dat',status='unknown')
      write(333,*)'Number_of_Atoms:            ',natom
      write(333,*)'Act_energy(kcal/mol):       0. '
      write(333,*)'Initial_Temperature:        50'
      write(333,*)'Temperature_steps:          40'
      write(333,*)'Temperature_increment:      40'
      write(333,*)'Delta_Energy_rea:           0.'
      write(333,*)'Delta_Energy_pro:           0.'
      write(333,*)'Maxstep:                   ',numpointstot
      write(333,*)'Npointsint:                 5 '
      write(333,*)'Maxtdev:                    0.5'
      write(333,*)'Rearrange(1=yes,0=no)       1'
      write(333,*)'SaddlePoint                ',numpointsf
         if(intfreq.eq.0)then
            write(333,*)'internalcoord(1=yes)     0'
         else if (intfreq.eq.1)then
            write(333,*)'internalcoord(1=yes)     1'            
         endif
      write(333,*)'isct_vtst(1=vtst_sct,0=sct) 1'
      write(333,*)'zerocurvature(1)            0'
      write(333,*)'reduced_mass                1.0'
      write(333,*)'minimum_frequency            50'
      write(333,*)'anim_freq(if_Maxstep=1)          2'
      write(333,*)'onlyrotors(0=yes,1=no)        1'
      write(333,*)'proj_rea_coo(0=yes(def),1=no) ',iprojrcoo
      write(333,*)'numrotors                     ',0

cc here starts the projection cycle over the total number of IRC points

      
      do inumpoints = 1, numpointstot
         open (unit=133,file='RPHt_input_data.dat',status='unknown')
c         open (unit=134,file='./data/hind_rot_head.dat',
c     +         status='unknown')

         write(133,*)'Number_of_Atoms: ',natom
         write(133,*)'Act_energy(kcal/mol):       0. '
         write(133,*)'Initial_Temperature:        200'
         write(133,*)'Temperature_steps:          40'
         write(133,*)'Temperature_increment:      40'
         write(133,*)'Delta_Energy_rea:           0.'
         write(133,*)'Delta_Energy_pro:           0.'
         write(133,*)'Maxstep:                    1'
         write(133,*)'Npointsint:                 5 '
         write(133,*)'Maxtdev:                    0.5'
         write(133,*)'Rearrange(1=yes,0=no)       1'
         write(133,*)'SaddlePoint                 1'
         if(intfreq.eq.0)then
            write(133,*)'internalcoord(1=yes)     0'
         else if (intfreq.eq.1)then
            write(133,*)'internalcoord(1=yes)     1'            
         endif
         write(133,*)'isct_vtst(1=vtst_sct,0=sct) 1'
         write(133,*)'zerocurvature(1)            0'
         write(133,*)'reduced_mass                1.0'
         write(133,*)'minimum_frequency            50'
         write(133,*)'anim_freq(if_Maxstep=1)       2'
         write(133,*)'onlyrotors(0=yes,1=no)        0'


c         do iatom = 1, 17
c            read (134,'(A60)') atomlabel(iatom) 
c            write (133,'(A60)') atomlabel(iatom) 
c         enddo
c         close(134)
         
         if(nhind.ne.0) then
            open (unit=15,file='./output/hrdata4proj_ts.dat'
     $           ,status='unknown')
            read (15,*)cjunk
            if(inumpoints.eq.numpointsf+1)then
               write (133,*)'proj_rea_coo(0=yes(def),1=no) ',1
            else
               write (133,*)'proj_rea_coo(0=yes(def),1=no) ',iprojrcoo
            endif
            read (15,*)cjunk,nhind
            write (133,*)cjunk,nhind

            do ir=1,nhind
               read (15,*)cjunk,ipivotA(ir)
               write (133,*)cjunk,ipivotA(ir)

