      program nst

c     A very simple minimum-on-the-seam-of-crossings (MSX) optimizer 
c     and nonadiabatic statistical theory (NST) flux calculator
c     Jasper, 2016
c
c     The MSX optimization strategy is from
c     1. M. J. Bearpark, M. A. Robb, and H. B. Schlegel, Chem. Phys. Lett. 223, 269 (1994)
c
c     The effective Hessian used in the NST calculation is from
c     2. J. N. Harvey and M. Aschi, Phys. Chem. Chem. Phys. 1, 5555-5563 (1999)
c
c     NST (often called NA TST) has been described
c     3. J. N. Harvey, Phys. Chem. Chem. Phys. 9, 331−343 (2007)
c     4. A. W. Jasper, J. Phys. Chem. A 119, 7339-7351 (2015)

c      implicit double precision(a-h,o-z)
      implicit none
      integer::i,j,n1,n2,nargs,maxstep,mnat,mnsurf,nclu,is,nsurft,
     &  nupper,nlower,qcprog
      double precision::etol,gtol,autoang,amutoau,autokcal,
     &  alpha_gf,alpha_gf_last,alpha_gg,alpha_gg_last,
     &  delg_dot_delg_gf,delg_dot_delg_gg,delx_dot_delg_gf,
     &  delx_dot_delg_gg,dot1,eavg,ediff1,ediff2,egap,elast1,elast2,
     &  etmp1,etmp2,stepscale,x1mag,zero,dot2,glownorm,ggnorm,
     &  glnlast,ggnlast,glndiff,ggndiff,exctole,exctolg,excscal,
     &  egaplast,egapdiff,ediffavg
      parameter(mnat=13)    ! maximum number of atoms
      parameter(mnsurf=3)   ! maximum number of electronic surfaces in the PES call
      parameter(autoang=0.52917706d0) ! bohr to A
      parameter(amutoau=1822.844987d0) ! amu to au
      parameter(autokcal=627.509d0) ! hartree to kcal/mol

      logical ldonej,firstcycle,repeatcyc,barbor
      character*2 symb(mnat)
      character*4,allocatable::args(:)
      double precision mm(mnat)
      double precision xx(3,mnat),x(mnat),y(mnat),z(mnat),
     & gpemd(3,mnat,mnsurf,mnsurf),pemd(mnsurf,mnsurf),grad(3,mnat),
     & x1(3,mnat),gupp(3,mnat),glow(3,mnat),gg(3,mnat),gf(3,mnat),
     & diffx_gf(mnat),diffy_gf(mnat),diffz_gf(mnat),
     & diffx_gg(mnat),diffy_gg(mnat),diffz_gg(mnat),
     & gglast(3,mnat),gflast(3,mnat),diffgf(3,mnat),diffgg(3,mnat)

      firstcycle=.true.
      repeatcyc=.false.
      barbor=.true.  ! Turn on Barzilai-Borwein step-size control by default
      qcprog=1       ! Gaussian as default program package
      
      nargs=COMMAND_ARGUMENT_COUNT()
      if (nargs .gt. 0) then
      allocate(args(nargs))
      do i=1,nargs
        call GET_COMMAND_ARGUMENT (i, args(i))
        args(i)=adjustl(args(i))
        if (trim(args(i))=='nobb' .or. trim(args(i))=='NOBB' .or.
     &    trim(args(i))=='noBB' .or. trim(args(i))=='NoBB') then
          write(*,*) 'Barzilai-Borwein method switched off.'
          barbor=.false.
        elseif (trim(args(i))=='m' .or. trim(args(i))=='M') then
          write(*,*) 'Using Molpro interface instead of Gaussian.'
          qcprog=2
        endif
      enddo
      deallocate(args)
      endif

c set these
      etol = 5.d-5   ! energy gap convergence tolerance in kcal/mol
      gtol = 1.d-5   ! gradient norm convergence tolerance in Hartree/bohr
      maxstep = 1000 ! maximum number of steps
      stepscale = 1. ! scales stepsize
      n1=1           ! PES call index for surface 1
      n2=2           ! PES call index for surface 2
      ldonej=.true.  ! Calculate and write N_ISC(E,J)?
      exctole = 10.0d0  ! energy tolerance for excessive step identification (kcal/mol)
      exctolg = 0.3d0  ! gradient norm tolerance for excessive step id. (Hartree/bohr)
      excscal = 0.5d0   ! scale down step by this amount if excessive step was detected

c read initial guess
      read(5,*)
      read(5,*)nclu,zero
      if (nclu.gt.mnat) then
        print *,"nclu (",nclu,") > mnat (",mnat,")"
        stop
      endif
      do i=1,nclu
      read(5,*)symb(i),mm(i),(xx(j,i),j=1,3)
      x(i)=xx(1,i)/autoang
      y(i)=xx(2,i)/autoang
      z(i)=xx(3,i)/autoang
      mm(i)=mm(i)*amutoau
      enddo
      read(5,*)ldonej

c prepare unit 80 for molden style optimization movie
      write(80,*)"[Molden Format]"
      write(80,*)"[GEOMETRIES] XYZ"

c dont change
      nsurft=2

c initialize
      elast1 = 0.d0
      elast2 = 0.d0
      ggnorm = 0.d0
      glownorm=0.d0
      egap =0.d0

c header
      write(6,*)"MSX optimization"
      write(6,*)
      write(6,*)"Energies in kcal/mol"
      write(6,'(a5,12a15)')"step","E1","E2","<E>","E1-E2","dE1","dE2",
     &  "gr.norm(lo)","gr.norm(up)"

c main loop
      do is=1,maxstep
 101  continue
c call pes
      call pot(symb,x,y,z,pemd,gpemd,nclu,mnat,nsurft,mnsurf,qcprog)
      pemd(n1,n1)=pemd(n1,n1)-zero
      pemd(n2,n2)=pemd(n2,n2)-zero
      
      egaplast=egap
      egap=pemd(n1,n1)-pemd(n2,n2)

c calculate energy differences and average
      etmp1 = pemd(n1,n1)
      etmp2 = pemd(n2,n2)
      ediff1 = (elast1-etmp1)*autokcal
      ediff2 = (elast2-etmp2)*autokcal
      eavg=(etmp1+etmp2)/2.d0*autokcal
      egapdiff=(dabs(egap)-dabs(egaplast))*autokcal
      ediffavg=(ediff1+ediff2)/2.d0

c identify higher E state
      nlower=n1
      nupper=n2
      if (pemd(n1,n1).gt.pemd(n2,n2)) then
        nupper=n1
        nlower=n2
      endif

c get x1 and gupper
      x1mag = 0.d0
      do i=1,nclu
      do j=1,3
      x1(j,i)   = gpemd(j,i,n1,n1)-gpemd(j,i,n2,n2)    ! grad(E1-E2)
      x1mag     = x1mag+x1(j,i)**2
      gupp(j,i) = gpemd(j,i,nupper,nupper)         ! grad of upper surface
      glow(j,i) = gpemd(j,i,nlower,nlower)         ! grad of lower surface
      enddo
      enddo
      x1mag = dsqrt(x1mag)                         ! magnitude of grad(E1-E2)

c ****
c for future improvement to treat nonzero nonadiabatic coupling...
c     get or compute d
c     compute x2(j,i)=d(j,i)/dmag  ! unit vec in direction of d
c     compute cross product of x2 and x1/x1mag (cross prod of the two unit vectors)
c     the resulting unit vector (xp) points in the direction we want from gupp
c     project the rest out from gupp (take the dot product of gupp and xd (dot2))
c     gg(j,i) = dot2 * xp(j,i)
c     gf(j,i) is the same as belose
c     add gg and gf to get the total gradient
c ****

c for zero nonadiabatic coupling...
c compute dot prod
      dot1=0.d0
      dot2=0.d0
      do i=1,nclu
      do j=1,3
        dot1=dot1+x1(j,i)*gupp(j,i)/x1mag           ! grad upper . x1/x1mag
        dot2=dot2+x1(j,i)*glow(j,i)/x1mag           ! grad lower . x1/x1mag
      enddo
      enddo
      
c project out componenets in the direction of x1
      do i=1,nclu
      do j=1,3
      gflast(j,i)=gf(j,i)
      gf(j,i) = 2.d0*egap*x1(j,i)/x1mag            ! f = 2*(E1-E2)*x1/x1mag
      gglast(j,i)=gg(j,i)
      gg(j,i)   = gupp(j,i)-dot1*x1(j,i)/x1mag       ! g = grad upper projected onto plane perp to f
      glow(j,i) = glow(j,i)-dot2*x1(j,i)/x1mag
cccc      grad(j,i) = gg(j,i)+gf(j,i)                  ! total gradient (now computed below!)
      enddo
      enddo

c calculate gradient norms of gg and glow (for convergence check)

      ggnlast = ggnorm
      glnlast = glownorm
      ggnorm    = dsqrt(sum(gg(1:3,1:nclu)*gg(1:3,1:nclu)))
      glownorm  = dsqrt(sum(glow(1:3,1:nclu)*glow(1:3,1:nclu)))
      ggndiff = ggnorm-ggnlast
      glndiff = glownorm-glnlast
      
c check if the last last step caused an excessive increase in energy (gap) and/or gradient norms
c reminder: ediff was defined as elast-ecurrent (postive sign is desired)
c If an excessive step is detected, it is rejected and replaced by one that is scaled
c down significantly (defined by excscal parameter)
      if (barbor .and. .not.(firstcycle.or.repeatcyc)) then
      if ( (egapdiff .gt. exctole .and. -ediffavg.gt.exctole)
     &    .or. ggndiff.gt.exctolg .or. glndiff .gt.exctolg ) then
        write(*,*) 'Excessive step detected. Scaling down & goin back'
        do i=1,nclu
          x(i) = x(i) + (1.d0-excscal)*grad(1,i)
          y(i) = y(i) + (1.d0-excscal)*grad(2,i)
          z(i) = z(i) + (1.d0-excscal)*grad(3,i)
        enddo
        ggnorm=ggnlast
        glownorm=glnlast
        gf = gflast
        gg = gglast
        repeatcyc=.true.
c        firstcycle=.true.
        goto 101
      endif
      endif
c      repeatcyc=.false.
      
      if (firstcycle .or. .not. barbor) then
c  do an ordinary step in the first iteration (and subsequent if Barzilai-Borwein method is switched off)
        alpha_gg=stepscale
        alpha_gf=stepscale
        firstcycle=.false.
      else
c  Use Barzilai-Borwein step-size control in the subsequent steps
c  form difference vectors delta_x and delta_g
        alpha_gf_last=alpha_gf
        alpha_gg_last=alpha_gg
        do i=1,nclu
          diffx_gf(i) = -alpha_gf_last*gflast(1,i)
          diffx_gg(i) = -alpha_gg_last*gglast(1,i)
          diffy_gf(i) = -alpha_gf_last*gflast(2,i)
          diffy_gg(i) = -alpha_gg_last*gglast(2,i)
          diffz_gf(i) = -alpha_gf_last*gflast(3,i)
          diffz_gg(i) = -alpha_gg_last*gglast(3,i)
        enddo
        if (repeatcyc) then
          diffx_gf=excscal*diffx_gf
          diffx_gg=excscal*diffx_gg
          diffy_gf=excscal*diffy_gf
          diffy_gg=excscal*diffy_gg
          diffz_gf=excscal*diffz_gf
          diffz_gg=excscal*diffz_gg
        endif
        diffgf=gf-gflast
        diffgg=gg-gglast
c  compute step sizes alpha_g=<delta_x|delta_g>g / <delta_g|delta_g>g
c  and                alpha_f=<delta_x|delta_g>f / <delta_g|delta_g>f
        delx_dot_delg_gg=0.0d0
        delx_dot_delg_gf=0.0d0
        delg_dot_delg_gg=0.0d0
        delg_dot_delg_gf=0.0d0
        do i=1,nclu
          do j=1,3
            if     (j.eq.1) then 
              delx_dot_delg_gg=delx_dot_delg_gg+
     &                         diffx_gg(i)*diffgg(j,i)
              delx_dot_delg_gf=delx_dot_delg_gf+
     &                         diffx_gf(i)*diffgf(j,i)
            elseif (j.eq.2) then
              delx_dot_delg_gg=delx_dot_delg_gg+
     &                         diffy_gg(i)*diffgg(j,i)
              delx_dot_delg_gf=delx_dot_delg_gf+
     &                         diffy_gf(i)*diffgf(j,i)
            else
              delx_dot_delg_gg=delx_dot_delg_gg+
     &                         diffz_gg(i)*diffgg(j,i)
              delx_dot_delg_gf=delx_dot_delg_gf+
     &                         diffz_gf(i)*diffgf(j,i)
            endif
            delg_dot_delg_gg=delg_dot_delg_gg + diffgg(j,i)**2
            delg_dot_delg_gf=delg_dot_delg_gf + diffgf(j,i)**2
          enddo
        enddo
        alpha_gg= delx_dot_delg_gg / delg_dot_delg_gg
        alpha_gf= delx_dot_delg_gf / delg_dot_delg_gf
        firstcycle=.false.
      endif
      
      grad = alpha_gg * gg + alpha_gf * gf  ! total gradient
      
      repeatcyc=.false.
      
c take step
      do i=1,nclu
      x(i) = x(i) - grad(1,i)
      y(i) = y(i) - grad(2,i)
      z(i) = z(i) - grad(3,i)
      enddo

c save the current energies (for convergence check in next cycle)
      elast1 = etmp1
      elast2 = etmp2
      
      write(6,100)is,pemd(n1,n1)*autokcal,pemd(n2,n2)*autokcal,eavg,
     & egap*autokcal,ediff1,ediff2,glownorm,ggnorm

      write(80,*)nclu 
      write(80,*)is,egap*autokcal
      do i=1,nclu
      write(80,'(a,3f15.5)')symb(i),
     &      x(i)*autoang,y(i)*autoang,z(i)*autoang
      enddo

      if (ggnorm.lt.gtol.and.glownorm.lt.gtol
     &               .and.(dabs(egap)*autokcal).lt.etol) go to 999
     
c check if energies haven't changed even though convergence is not
c achieved yet. If yes, do a regular gradient descent step in the next
c cycle
      if (dabs(ediff1).lt.etol .and. dabs(ediff2).lt.etol
     &     .and. barbor) then
        firstcycle=.true.
      endif

      enddo

 999  continue

      write(6,*)
      write(6,*)"Optimized geometry (A)"
      do i=1,nclu
      write(6,'(a,3f22.15)')symb(i),
     &      x(i)*autoang,y(i)*autoang,z(i)*autoang
      enddo
      write(6,*)

 100  format(i5,3f15.8,3f15.8,2f15.8)

      do i=1,nclu
       xx(1,i)=x(i)
       xx(2,i)=y(i)
       xx(3,i)=z(i)
      enddo

c     calculate frequencies on the seam
      if (ldonej) call msxfreq(xx,mm,nclu,symb,zero,n1,n2,qcprog)

      print *,"Well done!"
      end


