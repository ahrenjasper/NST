      subroutine msxfreq(xx0,mm,nclu,symb,el_zero,jobtype,idebug,icut,
     &                    es,emax,js,jmax,hso12,sc_qelec)

      implicit none

c     Passed into the subroutine
      double precision xx0(3,nclu),mm(nclu)
      integer nclu
      character*2 symb(nclu)
      character*5 jobtype
      logical idebug,icut

c     from input file
      double precision el_zero
      double precision es,emax
      integer js,jmax
      double precision hso12
      double precision sc_qelec

c     Phys constants and conversions
      double precision :: pi=dacos(-1.d0)
      double precision autoang,autokcal,autos
      double precision autoev,amutoau,autocmi
      double precision mu,kb
      parameter(autoang=0.52917706d0)
      parameter(autocmi=219474.63067d0)
      parameter(autokcal=627.509d0)
      parameter(autos=2.4188843d-17)
      parameter(amutoau=1822.844987d0) 
      parameter(autoev=27.2113961d0)
      parameter(kb=3.166829d-6)     
      parameter(mu=1.d0*amutoau)

c     Grid and surf parameters
      integer mnsurf,mgrid
      parameter(mnsurf=3)
      parameter(mgrid=100000)
      integer nsurf
      integer mnclu

c     Variables Needed for msxfreq
      integer :: n1 = 1
      integer :: n2 = 2

c     Values initialized to Zero
      double precision :: ezero = 0.0
      double precision :: emin = 0.0
      integer :: jmin = 0
  
c     For linear checks      
      logical linear
      double precision linearinfinity

c     mass variables 
      double precision mtot,mx,my,mz

c     coordinates variables
      double precision xx(3,nclu)
      double precision x(nclu),y(nclu),z(nclu)
      double precision xxx
      double precision be

c     moment-of-intertia variables
      double precision rot(3,3),mom(3,3),momi(3),mom1,mom2
      double precision eig(3),ap(6)

c     gradients variables     
      double precision pemd(mnsurf,mnsurf)
      double precision gpemd(3,nclu,mnsurf,mnsurf)
      double precision gv1a(3,nclu),gv2a(3,nclu)
      double precision gv1b(3,nclu),gv2b(3,nclu)
      double precision gperp1(3*nclu),gperp2(3*nclu),gperp3(3*nclu)
      double precision hh
      double precision :: gdot = 0.0

c     Hessian variables
      integer nmax,ndim
      double precision hessa(3*nclu,3*nclu),hessb(3*nclu,3*nclu)
      double precision mwhessa(3*nclu,3*nclu),mwhessb(3*nclu,3*nclu)
      double precision hessax(3*nclu,3*nclu)
      double precision hessval

c     Frequency variables
      integer nfreq
      double precision freq(3*nclu)

c     State-count variables
      double precision ee(50000)
      integer ne
      double precision fx(3*nclu)
      integer is,im
      double precision t(mgrid),at(mgrid),estep,ejk,djk,eee
      integer imax,iejk,kmax,ifreq
      double precision etmp,ss
      double precision plz,rho,rholz,p2pass
     
c     Variables for matrix diagonalization routines
      integer info
      integer lwork
      double precision work(9*nclu-1),work2(9)

c     Loop Variables
      integer i,j,ii,ij,jj,k,l,kl

c     Variables for Temp Storage
      double precision tmp,tmp1,tmp2,tmp3
      double precision temp1,temp2,temp3

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
c     Initialize certain variables as necessary

      lwork=9*nclu-1

      do i=1,3*nclu
      do j=1,3*nclu
      hessa(i,j) = 0.d0
      hessb(i,j) = 0.d0
      enddo
      enddo

      mtot=0.d0
      do i=1,nclu
      mtot=mtot+mm(i)
      enddo
      
      do i=1,3
        do j=1,nclu
          xx0(i,j)=xx0(i,j)/autoang
        enddo
      enddo

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     Calculate Moments-of-Inertia; Reorient Geometry

c     center & reorient
      mx=0.d0
      my=0.d0
      mz=0.d0
      do i=1,nclu
      mx=mx+mm(i)*xx0(1,i)/mtot
      my=my+mm(i)*xx0(2,i)/mtot
      mz=mz+mm(i)*xx0(3,i)/mtot
      enddo
      do i=1,nclu
       xx0(1,i)=xx0(1,i)-mx
       xx0(2,i)=xx0(2,i)-my
       xx0(3,i)=xx0(3,i)-mz
      enddo

c     compute moment of intertia matrix (mom)
      do i=1,3
      do j=1,3
          mom(i,j) = 0.d0
      enddo
      enddo

      do i=1,nclu
         mom(1,1)=mom(1,1)+mm(i)*(xx0(2,i)**2+xx0(3,i)**2)
         mom(2,2)=mom(2,2)+mm(i)*(xx0(1,i)**2+xx0(3,i)**2)
         mom(3,3)=mom(3,3)+mm(i)*(xx0(1,i)**2+xx0(2,i)**2)
         mom(1,2)=mom(1,2)-mm(i)*(xx0(1,i)*xx0(2,i))
         mom(1,3)=mom(1,3)-mm(i)*(xx0(1,i)*xx0(3,i))
         mom(2,3)=mom(2,3)-mm(i)*(xx0(2,i)*xx0(3,i))
      enddo
      mom(2,1)=mom(1,2)
      mom(3,1)=mom(1,3)
      mom(3,2)=mom(2,3)

c     diagonalize the mom matrix
      do i=1,3
      do j=i,3
        ap(i+(j-1)*j/2)=mom(i,j)
      enddo
      enddo
      call dspev('v','u',3,ap,eig,rot,3,work2,info )

      do i=1,3
      momi(i)=0.5d0/eig(i)
      enddo
      write(6,'(a,3f15.5)')
     & " Moments of intertia (cm-1)",(autocmi*momi(i),i=1,3)
      tmp1=dabs(momi(1)-momi(2))
      tmp2=dabs(momi(2)-momi(3))
      tmp3=dabs(momi(1)-momi(3))
      if (tmp1.lt.tmp2.and.tmp1.lt.tmp3) then
      mom2=(momi(1)+momi(2))/2.d0
      mom1=momi(3)
      elseif (tmp2.lt.tmp3.and.tmp2.lt.tmp1) then
      mom2=(momi(2)+momi(3))/2.d0
      mom1=momi(1)
      else
      mom2=(momi(1)+momi(3))/2.d0
      mom1=momi(2)
      endif

c     check if the molecule is linear using mom eigenvalues 
      linearinfinity=1.d3
      if(mom1.gt.linearinfinity) then
      write(6,*)" Linear species found!"
      linear=.true.
      write(6,'(a,1f15.5,a)')" Symmetrized to (cm-1)     ",
     &    mom2*autocmi," (x2)"
      else
      linear=.false.
      write(6,'(a,2f15.5,a)')" Symmetrized to (cm-1)     ",
     &    mom1*autocmi,mom2*autocmi," (x2)"
      endif
      write(6,*)

c     rotate to diagonalize mom
      do i=1,nclu
         temp1 = xx0(1,i)
         temp2 = xx0(2,i)
         temp3 = xx0(3,i)
         xx0(1,i)=temp1*rot(1,1)+temp2*rot(2,1)+temp3*rot(3,1)
         xx0(2,i)=temp1*rot(1,2)+temp2*rot(2,2)+temp3*rot(3,2)
         xx0(3,i)=temp1*rot(1,3)+temp2*rot(2,3)+temp3*rot(3,3)
      enddo

      do i=1,nclu
      x(i)=xx0(1,i)
      y(i)=xx0(2,i)
      z(i)=xx0(3,i)
      enddo

c      code does something???
c      nsurf=mnsurf
c      mnclu=nclu
c      call pot(symb,x,y,z,pemd,gpemd,nclu,mnclu,mnsurf,nsurf)
c      pemd(n1,n1)=pemd(n1,n1)-el_zero
c      pemd(n2,n2)=pemd(n2,n2)-el_zero
c      if (idebug) write(6,*)"angmom energy at MSX = ",
c     & pemd(n1,n1)*627.509,pemd(n2,n2)*627.509
c      ethresh=pemd(n1,n1)
c      hso12=pemd(n1,n2)

c     rotate geometry to align with principal axes of rotation
      write(6,*)"rotated geometry"
      do i=1,nclu
      write(6,'(a,3f22.15)')symb(i),
     &      x(i)*autoang,y(i)*autoang,z(i)*autoang
      enddo

c     kill the calculation if only rotation is requested      
      if (jobtype.eq."ROTGM") then
        stop
      endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     Calculate the gradients

c     gradient of the gap
      write(6,*)
      write(6,*)
      do i=1,nclu
      x(i)=xx0(1,i)
      y(i)=xx0(2,i)
      z(i)=xx0(3,i)
      enddo
      nsurf=mnsurf
      mnclu=nclu
      call pot(symb,x,y,z,pemd,gpemd,nclu,mnclu,mnsurf,nsurf)
      pemd(n1,n1)=pemd(n1,n1)-el_zero
      pemd(n2,n2)=pemd(n2,n2)-el_zero
      write(6,*)"energy at MSX = ",pemd(n1,n1)*627.509,pemd(n2,n2)*627.509
      write(6,*)
      tmp1=0.d0
      tmp2=0.d0
      tmp3=0.d0
      do i=1,3
      do j=1,nclu
        tmp1=tmp1+(gpemd(i,j,n1,n1)**2)*mu/mm(j)
        tmp2=tmp2+(gpemd(i,j,n2,n2)**2)*mu/mm(j)
        tmp3=tmp3+((gpemd(i,j,n1,n1)-gpemd(i,j,n2,n2))**2)*mu/mm(j)
      enddo
      enddo
      tmp1=dsqrt(tmp1)
      tmp2=dsqrt(tmp2)
      tmp3=dsqrt(tmp3)
 144  format(i3,7e15.5)
      write(6,*)"grad(E1) grad(E2) grad(E1-E2)"
      do i=1,3
      do j=1,nclu
        ij = (i-1)*nclu + j
        gperp1(ij)=gpemd(i,j,n1,n1)*dsqrt(mu/mm(j))/tmp1
        gperp2(ij)=gpemd(i,j,n2,n2)*dsqrt(mu/mm(j))/tmp2
        gperp3(ij)=(gpemd(i,j,n1,n1)-gpemd(i,j,n2,n2))
     &    *dsqrt(mu/mm(j))/tmp3
      write(6,'(3i5,3f15.5)')ij,i,j,gperp1(ij),gperp2(ij),gperp3(ij)
      enddo
      enddo
      write(6,*)

c     calculate the dot products
      do i=1,3
        do j=1,nclu
          gdot = gdot + gpemd(i,j,n1,n1) * gpemd(i,j,n2,n2) 
        enddo
      enddo

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c       Obtain the Hessians

c     either calculates Hess with NST or reads them from hess.x file
c     procedure is determined based on user input
      write(6,*)"Obtaining the Hessians"
      write(6,*)jobtype
      if (jobtype.eq."HESSC") then
        write(6,*)"Calculating Hessians..."
        hh = 0.01
        do i=1,3
        do j=1,nclu
          write(6,'(4(a,i5))')" step (",i,",",j,") of ( 3,",nclu,")"
          ij = (i-1)*nclu + j
          do k=1,3
          do l=1,nclu
            xx(k,l) = xx0(k,l)
          enddo
          enddo
          xx(i,j) = xx0(i,j) + hh
        do k=1,nclu
        x(k)=xx(1,k)
        y(k)=xx(2,k)
        z(k)=xx(3,k)
        enddo
        call pot(symb,x,y,z,pemd,gpemd,nclu,mnclu,mnsurf,nsurf)
          do k=1,3
          do l=1,nclu
            gv1a(k,l) = gpemd(k,l,n1,n1)
            gv1b(k,l) = gpemd(k,l,n2,n2)
          enddo
          enddo
          xx(i,j) = xx0(i,j) - hh
        do k=1,nclu
        x(k)=xx(1,k)
        y(k)=xx(2,k)
        z(k)=xx(3,k)
        enddo
        call pot(symb,x,y,z,pemd,gpemd,nclu,mnclu,mnsurf,nsurf)
            do k=1,3
            do l=1,nclu
              gv2a(k,l) = gpemd(k,l,n1,n1)
              gv2b(k,l) = gpemd(k,l,n2,n2)
            enddo
            enddo
          do k=1,3
          do l=1,nclu
            kl = (k-1)*nclu + l
c           Hessian matrix, 3NCLU X 3NCLU matrix
c           data ordered (x1,x2,...,y1,...,z1,...,zNCLU)
            hessa(ij,kl) = (gv1a(k,l) - gv2a(k,l))/(2.d0*hh)
            hessb(ij,kl) = (gv1b(k,l) - gv2b(k,l))/(2.d0*hh)
c           mass-scale
            mwhessa(ij,kl) = hessa(ij,kl)*mu/dsqrt(mm(j)*mm(l))
            mwhessb(ij,kl) = hessb(ij,kl)*mu/dsqrt(mm(j)*mm(l))
          enddo
          enddo
        enddo
        enddo
      elseif (jobtype.eq."HESSR") then
        write(6,*) "Now reading the Hessians"      
        open(30, file="hess.1")
        do i=1,3
          do j=1,nclu
            do k=1,3
              do l=1,nclu
                ij = nclu*(j-1)+i
                kl = nclu*(l-1)+k
                read(30,*) hessval
                mwhessa(ij,kl)=hessval*mu/dsqrt(mm(i)*mm(k))
              enddo
            enddo
          enddo
        enddo
        open(31, file="hess.3")
        do i=1,3
          do j=1,nclu
            do k=1,3
              do l=1,nclu
                ij=i*j
                kl=k*l
                ij = nclu*(j-1)+i
                kl = nclu*(l-1)+k
                read(31,*) hessval
                mwhessb(ij,kl)=hessval*mu/dsqrt(mm(i)*mm(k))
              enddo
            enddo
          enddo
        enddo
      endif

      write(6,*)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     Calculate a bunch of frequencies from Hessians for debugging

      if (idebug.eqv..true.) then
c       calc freqs
        write (6,*)"SURFACE 1"
        ndim = 3*nclu
        nmax = 3*nclu
        do i=1,ndim
        do j=1,ndim
          hessax(i,j)=mwhessa(i,j)
        enddo
        enddo
        call dsyev( 'v','u',ndim,hessax,nmax,freq,work,lwork,info)
        write(6,*)"  Index  Force Const(mass-scaled Eh/a0^2) Freq(cm-1)"
        do k=1,ndim
          if (freq(k).gt.0.d0) then
            tmp=dsqrt(freq(k)/mu)
          else
            tmp=-dsqrt(-freq(k)/mu)
          endif
          write(6,150)k,freq(k),tmp*autocmi
        enddo
        write(6,*)
c       calc freqs
        write (6,*)"SURFACE 1 projected"
        ndim = 3*nclu
        nmax = 3*nclu
        do i=1,ndim
        do j=1,ndim
        hessax(i,j)=mwhessa(i,j)
        enddo
        enddo
        call proj(symb,xx0,mm,mu,nclu,nsurf,el_zero,gperp1,hessax)
        call dsyev('v','u',ndim,hessax,nmax,freq,work,lwork,info)
        write(6,*)"  Index  Force Const(mass-scaled Eh/a0^2) Freq(cm-1)"
        do k=1,ndim
          if (freq(k).gt.0.d0) then
              tmp=dsqrt(freq(k)/mu)
          else
              tmp=-dsqrt(-freq(k)/mu)
          endif
          write(6,150)k,freq(k),tmp*autocmi
        enddo
        write(6,*)
c       calc freqs
        write (6,*)"SURFACE 2"
        ndim = 3*nclu
        nmax = 3*nclu
        do i=1,ndim
        do j=1,ndim
        hessax(i,j)=mwhessb(i,j)
        enddo
        enddo
        call dsyev('v','u',ndim,hessax,nmax,freq,work,lwork,info)
        write(6,*)"  Index  Force Const(mass-scaled Eh/a0^2) Freq(cm-1)"
        do k=1,ndim
          if (freq(k).gt.0.d0) then
              tmp=dsqrt(freq(k)/mu)
          else
              tmp=-dsqrt(-freq(k)/mu)
          endif
          write(6,150)k,freq(k),tmp*autocmi
        enddo
        write(6,*)
c       calc freqs
        write (6,*)"SURFACE 2 projected"
        ndim = 3*nclu
        nmax = 3*nclu
        do i=1,ndim
        do j=1,ndim
        hessax(i,j)=mwhessb(i,j)
        enddo
        enddo
        call proj(symb,xx0,mm,mu,nclu,nsurf,el_zero,gperp2,hessax)
        call dsyev('v','u',ndim,hessax,nmax,freq,work,lwork,info)
        write(6,*)"  Index  Force Const(mass-scaled Eh/a0^2) Freq(cm-1)"
        do k=1,ndim
          if (freq(k).gt.0.d0) then
              tmp=dsqrt(freq(k)/mu)
          else
              tmp=-dsqrt(-freq(k)/mu)
          endif
          if (idebug) write(6,150)k,freq(k),tmp*autocmi
        enddo
        write(6,*)
c       calc freqs
        write (6,*)"Surface 1 grad(E1-E2)-projected"
        ndim = 3*nclu
        nmax = 3*nclu
        do i=1,ndim
        do j=1,ndim
        hessax(i,j)=mwhessa(i,j)
        enddo
        enddo
        call proj(symb,xx0,mm,mu,nclu,nsurf,el_zero,gperp3,hessax)
        call dsyev('v','u',ndim,hessax,nmax,freq,work,lwork,info)
        write(6,*)"  Index  Force Const(mass-scaled Eh/a0^2) Freq(cm-1)"
        do k=1,ndim
          if (freq(k).gt.0.d0) then
              tmp=dsqrt(freq(k)/mu)
          else
              tmp=-dsqrt(-freq(k)/mu)
          endif
          write(6,150)k,freq(k),tmp*autocmi
        enddo
        write(6,*)
c       calc freqs
        write (6,*)"Surface 2 grad(E1-E2)-projected"
        ndim = 3*nclu
        nmax = 3*nclu
        do i=1,ndim
        do j=1,ndim
        hessax(i,j)=mwhessb(i,j)
        enddo
        enddo
        call proj(symb,xx0,mm,mu,nclu,nsurf,el_zero,gperp3,hessax)
        call dsyev('v','u',ndim,hessax,nmax,freq,work,lwork,info)
        write(6,*)"  Index  Force Const(mass-scaled Eh/a0^2) Freq(cm-1)"
        do k=1,ndim
          if (freq(k).gt.0.d0) then
              tmp=dsqrt(freq(k)/mu)
          else
              tmp=-dsqrt(-freq(k)/mu)
          endif
          write(6,150)k,freq(k),tmp*autocmi
        enddo
        write(6,*)
c       calc freqs
        write (6,*)"Effective two-state"
        ndim = 3*nclu
        nmax = 3*nclu
        do i=1,ndim
        do j=1,ndim
        hessax(i,j)=(mwhessb(i,j)+mwhessa(i,j))/2.d0
        enddo
        enddo
c        call proj(symb,xx0,mm,mu,nclu,nsurf,zero,gperp3,hessax)
        call dsyev('v','u',ndim,hessax,nmax,freq,work,lwork,info)
        write(6,*)"  Index  Force Const(mass-scaled Eh/a0^2) Freq(cm-1)"
        do k=1,ndim
          if (freq(k).gt.0.d0) then
              tmp=dsqrt(freq(k)/mu)
          else
              tmp=-dsqrt(-freq(k)/mu)
          endif
          if (idebug) write(6,150)k,freq(k),tmp*autocmi
        enddo
        write(6,*)
      endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     Calculate the Hessians used to determine the state counts

      write (6,*)"Effective two-state grad(E1-E2)-projected"
      write (6,*)
      write (6,*)


c     calculate the effective two-state Hessian from paper      
      write(6,*)"Effective Hessian from Harvey"
      ndim = 3*nclu
      nmax = 3*nclu
      do i=1,ndim
      do j=1,ndim
      if (gdot.gt.0.0) then
        hessax(i,j)=((tmp2*mwhessa(i,j)+tmp1*mwhessb(i,j))/(abs(tmp3)))
      elseif (gdot.lt.0.0) then
        hessax(i,j)=((tmp2*mwhessa(i,j)-tmp1*mwhessb(i,j))/(abs(tmp3)))
      endif
      enddo
      enddo
      call proj(symb,xx0,mm,mu,nclu,nsurf,el_zero,gperp3,hessax)
      call dsyev('v','u',ndim,hessax,nmax,freq,work,lwork,info)

      write(6,*)"   Index  Force Const (mass-scaled Eh/a0^2) Freq(cm-1)"
      do k=1,ndim
        if (freq(k).gt.0.d0) then
            tmp=dsqrt(freq(k)/mu)
        else
            tmp=-dsqrt(-freq(k)/mu)
        endif
        write(6,150)k,freq(k),tmp*autocmi
      enddo
      write(6,*)

c     simply average the two elements of each Hessain     
      write (6,*)"Simple Average of Hessian Elements of Two States"
      ndim = 3*nclu
      nmax = 3*nclu
      do i=1,ndim
        do j=1,ndim
         hessax(i,j)=(mwhessb(i,j)+mwhessa(i,j))/2.d0
        enddo
      enddo

      call proj(symb,xx0,mm,mu,nclu,nsurf,el_zero,gperp3,hessax)
      call dsyev( 'v','u',ndim,hessax,nmax,freq,work,lwork,info)

      write(6,*)"   Index  Force Const (mass-scaled Eh/a0^2) Freq(cm-1)"
      do k=1,ndim
        if (freq(k).gt.0.d0) then
            tmp=dsqrt(freq(k)/mu)
        else
            tmp=-dsqrt(-freq(k)/mu)
        endif
        write(6,150)k,freq(k),tmp*autocmi
      enddo
      write(6,*)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     Calculate cuts along each normal mode if requested

      if (icut.eqv..true.) then
      do im=0,ndim
        if (im.gt.0) print *,"mode=",im,dsqrt(dabs(freq(k))/mu)*autocmi
        if (im.eq.0) print *,"gradient"
        do is=-10,10
          ss=dble(is)*.2d0
          do i=1,3
          do j=1,nclu
            ij = (i-1)*nclu + j
            if (im.gt.0) xx(i,j) = xx0(i,j) + hessax(ij,im)*ss
            if (im.eq.0) xx(i,j) = xx0(i,j) + gperp3(ij)*ss
          enddo
          enddo
          do k=1,nclu
          x(k)=xx(1,k)
          y(k)=xx(2,k)
          z(k)=xx(3,k)
          enddo
          call pot(symb,x,y,z,pemd,gpemd,nclu,mnclu,mnsurf,nsurf)
          pemd(n1,n1)=pemd(n1,n1)-el_zero
          pemd(n2,n2)=pemd(n2,n2)-el_zero
          print *,ss,pemd(n1,n1)*autoev,pemd(n2,n2)*autoev
        enddo
      enddo
      ENDIF

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     Calculate the state counts and write the ne* files for MESS

c numerical grid for nej.dat
      ne=int(emax/es)
      emin=0.  
      emax=emax/autocmi
      emin=emin/autocmi
      jmin=0.  
      hso12=hso12/autocmi

c MSX properties
      nfreq=(3*nclu-7)
      if (linear) nfreq=nfreq+1
      print *,"MSX frequencies"
      ii=0
      rho=1.d0  ! calc classical harmonic rho prefactor
      do k=ndim-nfreq+1,ndim
       tmp=dsqrt(freq(k)/mu)  
       ii=ii+1
       fx(ii)=tmp
       rho=rho/tmp
       print *,tmp*autocmi
      enddo
      do i=2,nfreq-1
      rho=rho/dble(i)
      enddo
      print *

c NEJ.DAT
      open(33,file="nej.dat")
c     initialize density of state arrays t(i) (and at(i))
c     t(i) corresponds to energy bin from E = ESTEP*(i-1) to ESTEP*i
      estep=es/autocmi
      ezero=0.d0
      imax = int(emax/estep)+1
      write(33,*)imax,jmax/js+1,1.d0

      do jj=jmin,jmax,js

      do i=1,imax
      at(i) = 0.d0
      enddo
      if (linear) then
        djk = dble(2*jj+1)         ! degeneracy for tau
        ejk = mom2*dble(jj*(jj+1))
        iejk = int(ejk/estep)+1
        if (iejk.le.0) then
          write(6,*)"Rotational energy level found below energy grid"
          write(6,*)jj,i,iejk,ejk,estep
          stop
        endif
        if (iejk.le.imax) at(iejk) = at(iejk)+djk
      else
      do k=0,jj
        djk = dble(2*jj+1)        ! degeneracy for tau
        if (k.ne.0) djk=djk*2.d0  ! degeneracy for k
        ejk = mom2*dble(jj*(jj+1))+(mom1-mom2)*dble(k**2)
        iejk = int(ejk/estep)+1
        if (iejk.le.0) then
          write(6,*)"Rotational energy level found below energy grid"
          write(6,*)jj,i,iejk,ejk,estep
          stop
        endif
        if (iejk.le.imax) at(iejk) = at(iejk)+djk
      enddo
      endif
      do i=1,imax
        t(i)=at(i)
      enddo
      do j=1,nfreq
        kmax = int(emax/fx(j))  ! maximum quanta w/ E < EMAX
        do k=1,kmax               ! loop over allowed quanta
          ifreq = int(dble(k)*fx(j)/estep) ! vib level spacing in grid units
          do i=1,imax-ifreq
          at(i+ifreq)=at(i+ifreq)+t(i)
          enddo
        enddo
        do i=1,imax
          t(i)=at(i)
        enddo
      enddo

      do i=1,imax
        t(i)=at(i)
      enddo
      do i=1,imax
       ee(i)=dble(i)*estep
       rholz=0.d0
       do j=1,i ! convolute
        if (i.eq.j) then
          etmp = estep/2.d0
        else
          etmp = ee(i-j)+estep/2.d0
        endif
        plz=1.d0-dexp(-2.d0*pi*hso12**2/tmp3*dsqrt(0.5d0*mu/etmp)) 
        p2pass=plz+(1.d0-plz)*plz
c        p2pass=1.d0  ! test
        rholz=rholz+at(j)*p2pass
       enddo
       t(i)=rholz
      enddo

      eee=0.d0
      do while(eee.lt.ezero)
      write(33,133)eee*autocmi,jj,0.d0
      eee=eee+estep
      enddo
      write(33,133)eee*autocmi,jj,0.d0
      do i=1,imax
      eee=dble(i)*estep+ezero
      if (eee.le.emax) write(33,133)eee*autocmi,jj,t(i)*sc_qelec
      enddo

      enddo

c create NE.DAT from NEJ.DAT by just summing over J
      close(33)
      rewind(33)
      open(33,file="nej.dat")
      read(33,*)
      open(34,file="ne.dat")
      do i=1,imax
         t(i)=0.d0
      enddo
      do j=0,jmax,js
      do i=1,imax
      read (33,*)eee,jj,xxx
      t(i)=t(i)+xxx
      if (j.eq.jmax) write(34,134)eee,t(i)*js
      enddo
      enddo

      write(6,*)"Finished writing the nej.dat and ne.dat files."
      write(6,*)

 133  format(f12.2,i12,1pe20.8)
 134  format(f12.2,1pe20.8)

 101  format(i10,f15.5,10e15.5)
 150  format(i10,e15.5,f15.5)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      end subroutine msxfreq
