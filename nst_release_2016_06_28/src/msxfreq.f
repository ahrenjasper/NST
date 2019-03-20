      subroutine msxfreq(xx0,mm,nclu,symb,zero,n1,n2,qcprog)

c adapted from the normod subroutine in DiNT

      implicit none

      double precision autoang,autokcal,autos
      parameter(autoang=0.52917706d0)
      integer mnsurf,mnclu,nt,it,mgrid,jmax,js,jmin,qcprog
      parameter(mgrid=100000)
      double precision t(mgrid),at(mgrid),estep,ejk,djk,mom1,mom2,eee,es
      double precision t_lz(mgrid),t_ai(mgrid)
      integer jj,imax,iejk,kmax,ifreq
      double precision autoev,mu,amutoau,zero,autocmi,kb,ttt,qq,qqq,
     & qqlz,tt(100),qelec,etmp,xe,xmd
      parameter(kb=3.166829d-6)     ! hartree/K
      parameter(autocmi=219474.63067d0)
      parameter(autokcal=627.509d0)
      parameter(autos=2.4188843d-17)
      parameter(amutoau=1822.844987d0) ! best
      parameter(autoev=27.2113961d0)
      parameter(mnsurf=3)
      parameter(mu=1.d0*amutoau)  ! mass-scaling mass
c      parameter(mu=1.d0)  ! mass-scaling mass
      integer nclu,repflag,nsurf,nmtype,n1,n2
      double precision xx0(3,nclu),mm(nclu),tmp1,tmp2,tmp3,
     & ethresh
      character*2 symb(nclu)

      logical ldebug
      integer lwork,nbound,ne
      integer i,j,ij,k,i1,i2,l,kl,nmax,ndim,info,j1,j2,nfreq,nfw,ii
      double precision freq(3*nclu),ee(50000),emin,emax,h12,plz,pi,
     & rho,rholz_lz,rholz_ai,p2pass_lz,p2pass_ai,rhox,rhox0,ezero
      double precision xx(3,nclu),hh,ewell,evib
      double precision pemd(mnsurf,mnsurf),
     & gpemd(3,nclu,mnsurf,mnsurf),
     & gv1a(3,nclu),gv2a(3,nclu),
     & work(9*nclu-1),hessa(3*nclu,3*nclu),
     & tmp,klz,fw(3*nclu),wrho,fx(3*nclu),
     & x(nclu),y(nclu),z(nclu),gv1b(3,nclu),
     & gv2b(3,nclu),hessb(3*nclu,3*nclu),aaa1(3*nclu,3*nclu),
     & hessax(3*nclu,3*nclu),hessbx(3*nclu,3*nclu),
     & gperp1(3*nclu),gperp2(3*nclu),gperp3(3*nclu),
     & mtot,mx,my,mz,mom(3,3),ap(6),momi(3),
     * rot(3,3),eig(3),work2(9),temp1,temp2,temp3

      double precision a1,a2,a3,a4,xxx,be,e0,ex,tmpf,pairy,
     & airyarg,airypre
      
      double precision :: gdot = 0.0

      double precision tmpe,mnstEX,mnstEE,mnstA0,mnstA1,mnstE0,
     & mnstEC,pmnst

      lwork=9*nclu-1
      pi=dacos(-1.d0)
      ldebug=.true.
      ldebug=.false.

c Initialize
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
c center & reorient
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

c compute moment of intertia matrix mom
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
      call dspev( 'v','u',3,ap,eig,rot,3,work2,info )

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
      write(6,'(a,2f15.5,a)')" Symmetrized to (cm-1)     ",
     &    mom1*autocmi,mom2*autocmi," (x2)"
      write(6,*)

c rotate to diagonalize mom
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
      nsurf=mnsurf
      mnclu=nclu
      call pot(symb,x,y,z,pemd,gpemd,nclu,mnclu,mnsurf,nsurf,qcprog)
      pemd(n1,n1)=pemd(n1,n1)-zero
      pemd(n2,n2)=pemd(n2,n2)-zero
      if (ldebug) write(6,*)"angmom energy at MSX = ",
     & pemd(n1,n1)*627.509,pemd(n2,n2)*627.509
      ethresh=pemd(n1,n1)
      h12=pemd(n1,n2)

      if(ldebug) then
      write(6,*)"rotated geometry"
      do i=1,nclu
      write(6,'(a,3f22.15)')symb(i),
     &      x(i)*autoang,y(i)*autoang,z(i)*autoang
      enddo
      endif


c gradient of the gap
      do i=1,nclu
      x(i)=xx0(1,i)
      y(i)=xx0(2,i)
      z(i)=xx0(3,i)
      enddo
      nsurf=mnsurf
      mnclu=nclu
      call pot(symb,x,y,z,pemd,gpemd,nclu,mnclu,mnsurf,nsurf,qcprog)
      pemd(n1,n1)=pemd(n1,n1)-zero
      pemd(n2,n2)=pemd(n2,n2)-zero
      if (ldebug) write(6,*)"energy at MSX = ",
     &   pemd(n1,n1)*627.509,pemd(n2,n2)*627.509
      if (ldebug) write(6,*)
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
      if (ldebug) write(6,*)"grad(E1) grad(E2) grad(E1-E2)"
      do i=1,3
      do j=1,nclu
        ij = (i-1)*nclu + j
        gperp1(ij)=gpemd(i,j,n1,n1)*dsqrt(mu/mm(j))/tmp1
        gperp2(ij)=gpemd(i,j,n2,n2)*dsqrt(mu/mm(j))/tmp2
        gperp3(ij)=(gpemd(i,j,n1,n1)-gpemd(i,j,n2,n2))
     &    *dsqrt(mu/mm(j))/tmp3
        if (ldebug) 
     &  write(6,'(3i5,3f15.5)')ij,i,j,gperp1(ij),gperp2(ij),gperp3(ij)
      enddo
      enddo
     
c     calculate the dot products
      do i=1,3
        do j=1,nclu
          gdot = gdot + gpemd(i,j,n1,n1) * gpemd(i,j,n2,n2) 
        enddo
      enddo

c     write gradient information to output file
      write(6, *) "|grad(E1)| =  ",tmp1
      write(6, *) "|grad(E2)| =  ",tmp2
      write(6, *) "|grad(E1-E2)| =  ",tmp3
      write(6, *) "|grad(E1)| / |grad(E1-E2)| =  ",tmp1/tmp3
      write(6, *) "|grad(E2)| / |grad(E1-E2)| =  ",tmp2/tmp3
      write(6, *) "grad(E1) .dot. grad(E2) =  ",gdot
      write(6,*)

      if (ldebug) write(6,*)
      if (ldebug) write(6,*)

c calculate hessians of both states
c     stepsize
      write(6,*)"Calculating Hessians..."
      hh = 0.001d0
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
      call pot(symb,x,y,z,pemd,gpemd,nclu,mnclu,mnsurf,nsurf,qcprog)
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
      call pot(symb,x,y,z,pemd,gpemd,nclu,mnclu,mnsurf,nsurf,qcprog)
          do k=1,3
          do l=1,nclu
            gv2a(k,l) = gpemd(k,l,n1,n1)
            gv2b(k,l) = gpemd(k,l,n2,n2)
          enddo
          enddo
        do k=1,3
        do l=1,nclu
          kl = (k-1)*nclu + l
c         Hessian matrix, 3NCLU X 3NCLU matrix
c         data ordered (x1,x2,...,y1,...,z1,...,zNCLU)
          hessa(ij,kl) = (gv1a(k,l) - gv2a(k,l))/(2.d0*hh)
          hessb(ij,kl) = (gv1b(k,l) - gv2b(k,l))/(2.d0*hh)
c         mass-scale
          hessa(ij,kl) = hessa(ij,kl)*mu/dsqrt(mm(j)*mm(l))
          hessb(ij,kl) = hessb(ij,kl)*mu/dsqrt(mm(j)*mm(l))
        enddo
        enddo
      enddo
      enddo
      write(6,*)






      if (ldebug) write (6,*)"SURFACE 1"

      ndim = 3*nclu
      nmax = 3*nclu
      do i=1,ndim
      do j=1,ndim
      hessax(i,j)=hessa(i,j)
      enddo
      enddo

      call dsyev( 'v','u',ndim,hessax,nmax,freq,work,lwork,info )

      if (ldebug) 
     & write(6,*)"   Index  Force Const(mass-scaled Eh/a0^2) Freq(cm-1)"
      do k=1,ndim
        if (freq(k).gt.0.d0) then
            tmp=dsqrt(freq(k)/mu)
        else
            tmp=-dsqrt(-freq(k)/mu)
        endif
        if (ldebug) write(6,150)k,freq(k),tmp*autocmi
      enddo
      if (ldebug) write(6,*)






      if (ldebug) write (6,*)"SURFACE 1 projected"

      ndim = 3*nclu
      nmax = 3*nclu
      do i=1,ndim
      do j=1,ndim
      hessax(i,j)=hessa(i,j)
      enddo
      enddo

      call proj(symb,xx0,mm,mu,nclu,nsurf,zero,gperp1,hessax)
      call dsyev( 'v','u',ndim,hessax,nmax,freq,work,lwork,info )

      if (ldebug) 
     & write(6,*)"   Index  Force Const(mass-scaled Eh/a0^2) Freq(cm-1)"
      do k=1,ndim
        if (freq(k).gt.0.d0) then
            tmp=dsqrt(freq(k)/mu)
        else
            tmp=-dsqrt(-freq(k)/mu)
        endif
        if (ldebug) write(6,150)k,freq(k),tmp*autocmi
      enddo
      if (ldebug) write(6,*)




      if (ldebug) write (6,*)"SURFACE 2"

      ndim = 3*nclu
      nmax = 3*nclu
      do i=1,ndim
      do j=1,ndim
      hessax(i,j)=hessb(i,j)
      enddo
      enddo

      call dsyev( 'v','u',ndim,hessax,nmax,freq,work,lwork,info )

      if (ldebug) 
     & write(6,*)"   Index  Force Const(mass-scaled Eh/a0^2) Freq(cm-1)"
      do k=1,ndim
        if (freq(k).gt.0.d0) then
            tmp=dsqrt(freq(k)/mu)
        else
            tmp=-dsqrt(-freq(k)/mu)
        endif
        if (ldebug) write(6,150)k,freq(k),tmp*autocmi
      enddo
      if (ldebug) write(6,*)



      if (ldebug) write (6,*)"SURFACE 2 projected"

      ndim = 3*nclu
      nmax = 3*nclu
      do i=1,ndim
      do j=1,ndim
      hessax(i,j)=hessb(i,j)
      enddo
      enddo

      call proj(symb,xx0,mm,mu,nclu,nsurf,zero,gperp2,hessax)
      call dsyev( 'v','u',ndim,hessax,nmax,freq,work,lwork,info )

      if (ldebug) 
     & write(6,*)"   Index  Force Const(mass-scaled Eh/a0^2) Freq(cm-1)"
      do k=1,ndim
        if (freq(k).gt.0.d0) then
            tmp=dsqrt(freq(k)/mu)
        else
            tmp=-dsqrt(-freq(k)/mu)
        endif
        if (ldebug) write(6,150)k,freq(k),tmp*autocmi
      enddo
      if (ldebug) write(6,*)





      write (6,*)"Surface 1 grad(E1-E2)-projected"

      ndim = 3*nclu
      nmax = 3*nclu
      do i=1,ndim
      do j=1,ndim
      hessax(i,j)=hessa(i,j)
      enddo
      enddo

      call proj(symb,xx0,mm,mu,nclu,nsurf,zero,gperp3,hessax)
      call dsyev( 'v','u',ndim,hessax,nmax,freq,work,lwork,info )

      write(6,*)"   Index  Force Const(mass-scaled Eh/a0^2) Freq(cm-1)"
      do k=1,ndim
        if (freq(k).gt.0.d0) then
            tmp=dsqrt(freq(k)/mu)
        else
            tmp=-dsqrt(-freq(k)/mu)
        endif
        write(6,150)k,freq(k),tmp*autocmi
      enddo
      write(6,*)




      write (6,*)"Surface 2 grad(E1-E2)-projected"

      ndim = 3*nclu
      nmax = 3*nclu
      do i=1,ndim
      do j=1,ndim
      hessax(i,j)=hessb(i,j)
      enddo
      enddo

      call proj(symb,xx0,mm,mu,nclu,nsurf,zero,gperp3,hessax)
      call dsyev( 'v','u',ndim,hessax,nmax,freq,work,lwork,info )

      write(6,*)"   Index  Force Const(mass-scaled Eh/a0^2) Freq(cm-1)"
      do k=1,ndim
        if (freq(k).gt.0.d0) then
            tmp=dsqrt(freq(k)/mu)
        else
            tmp=-dsqrt(-freq(k)/mu)
        endif
        write(6,150)k,freq(k),tmp*autocmi
      enddo
      write(6,*)

      if (ldebug) write (6,*)"Effective two-state"

      ndim = 3*nclu
      nmax = 3*nclu
      do i=1,ndim
      do j=1,ndim
      hessax(i,j)=(hessb(i,j)+hessa(i,j))/2.d0
      enddo
      enddo

      call proj(symb,xx0,mm,mu,nclu,nsurf,zero,gperp3,hessax)
      call dsyev( 'v','u',ndim,hessax,nmax,freq,work,lwork,info )

      if (ldebug) 
     & write(6,*)"   Index  Force Const(mass-scaled Eh/a0^2) Freq(cm-1)"
      do k=1,ndim
        if (freq(k).gt.0.d0) then
            tmp=dsqrt(freq(k)/mu)
        else
            tmp=-dsqrt(-freq(k)/mu)
        endif
        if (ldebug) write(6,150)k,freq(k),tmp*autocmi
      enddo
      if (ldebug) write(6,*)


      write (6,*)"Effective two-state grad(E1-E2)-projected"

      ndim = 3*nclu
      nmax = 3*nclu
      do i=1,ndim
      do j=1,ndim
      hessax(i,j)=(hessb(i,j)+hessa(i,j))/2.d0
      enddo
      enddo

      call proj(symb,xx0,mm,mu,nclu,nsurf,zero,gperp3,hessax)
      call dsyev( 'v','u',ndim,hessax,nmax,freq,work,lwork,info )

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

c **** SYSTEM SPECIFIC PARAMETERS ****
c numerical grid for nej.dat
      read(5,*)es,emax
      ne=int(emax/es)
      emin=0.  
      emax=emax/autocmi
      emin=emin/autocmi
      read(5,*)js,jmax
      jmin=0.  
      read(5,*)h12,qelec
      h12=h12/autocmi

c MSX properties
      nfreq=(3*nclu-7)
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
      open(33,file="nej_lz.dat")
      open(44,file="nej_wc.dat")
c     initialize density of state arrays t(i) (and at(i))
c     t(i) corresponds to energy bin from E = ESTEP*(i-1) to ESTEP*i
      estep=es/autocmi
      ezero=0.d0
      imax = int(emax/estep)+1
      write(33,*)imax,jmax/js+1,1.d0
      write(44,*)imax,jmax/js+1,1.d0

      do jj=jmin,jmax,js

      do i=1,imax
      at(i) = 0.d0
      enddo
      do k=0,jj
        djk = dble(2*jj+1)         ! degeneracy for j
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
        t_lz(i)=at(i)
        t_ai(i)=at(i)
      enddo
      do i=1,imax
       ee(i)=dble(i)*estep
      enddo
      do i=1,imax
       rholz_lz=0.d0
       rholz_ai=0.d0
c       do j=1,i-1 ! convolute
       do j=1,imax ! convolute
        if (i.eq.j) then
          etmp = estep/2.d0
        elseif (i.lt.j) then
          etmp = -ee(j-i)+estep/2.d0
        else
          etmp = ee(i-j)+estep/2.d0
        endif
c       Landau-Zener transition probability
        if (j.le.i) then
        plz=1.d0-dexp(-2.d0*pi*h12**2/tmp3*dsqrt(0.5d0*mu/etmp))  ! LZ prob
        p2pass_lz=plz+(1.d0-plz)*plz
        else
        plz=0.d0
        p2pass_lz=0.d0
        endif
c        p2pass_lz=1.d0  ! test
c       Weak-coupling transition probability
        tmpf=dsqrt(tmp1*tmp2)
        e0=(tmpf**4/(2.d0*mu*tmp3**2))**(1.d0/3.d0)
        be=(2.d0*h12*tmpf/e0/tmp3)**(3.d0/2.d0)
        ex=etmp*tmp3/(2.d0*h12*tmpf)
        airyarg=-ex*(be**(2.d0/3.d0))
        airypre=pi**2*be**(4.d0/3.d0)
        call airya(airyarg,a1,a2,a3,a4)
        pairy=airypre*a1**2
c        p2pass_ai=pairy+(1.d0-pairy)*pairy
        rholz_lz=rholz_lz+at(j)*p2pass_lz
        rholz_ai=rholz_ai+at(j)*pairy
c       write(6,'(2i5,10f15.8)')i,j,mnstEX,mnstEE,
c     & mnstA0,mnstA1,mnstEC,mnstE0,pmnst
       enddo
       t_lz(i)=rholz_lz
       t_ai(i)=rholz_ai
      enddo

      eee=0.d0
      do while(eee.lt.ezero)
      write(33,133)eee*autocmi,jj,0.d0
      write(44,133)eee*autocmi,jj,0.d0
      eee=eee+estep
      enddo
      write(33,133)eee*autocmi,jj,0.d0
      write(44,133)eee*autocmi,jj,0.d0
      do i=1,imax
      eee=dble(i)*estep+ezero
      if (eee.le.emax) write(33,133)eee*autocmi,jj,t_lz(i)*qelec
      if (eee.le.emax) write(44,133)eee*autocmi,jj,t_ai(i)*qelec
      enddo

      enddo

c create NE.DAT from NEJ.DAT by just summing over J
      close(33)
      close(44)
      rewind(33)
      rewind(44)
      open(33,file="nej_lz.dat")
      read(33,*)
      open(44,file="nej_wc.dat")
      read(44,*)
      open(34,file="ne_lz.dat")
      open(45,file="ne_wc.dat")
      do i=1,imax
         t_lz(i)=0.d0
         t_ai(i)=0.d0
      enddo
      do j=0,jmax,js
      do i=1,imax
      read (33,*)eee,jj,xxx
      t_lz(i)=t_lz(i)+xxx
      read (44,*)eee,jj,xxx
      t_ai(i)=t_ai(i)+xxx
      if (j.eq.jmax) write(34,134)eee,t_lz(i)*js
      if (j.eq.jmax) write(45,134)eee,t_ai(i)*js
      enddo
      enddo

      write(6,*)"Finished writing the nej.dat and ne.dat files."
      write(6,*)

 133  format(f12.2,i12,1pe20.8)
 134  format(f12.2,1pe20.8)

 101  format(i10,f15.5,10e15.5)
 150  format(i10,e15.5,f15.5)
      end
