      subroutine integz(coor,ain,afin,alpha,nsteps,norder,np,
     #                 ng,irun,egtldetr,omegatot,rL,skip,frm,
     #                 pbchoice,hubble,entest,pmass,ng2d,gpscal)

!-----PM time-stepper with z array controlling printouts

      implicit none

!-----scalars

      integer nsteps,norder,np,ng,irun,ii,skip,frm,pbchoice
      integer nsloc,jj,nsum,entest,ng2d

      real tau1,tau2,tau3,tau4,eps,alf,pp,alpha,omegatot,rL
      real ain,afin,pfin,hubble,epsloc,ptime,delp,prange
      real psum,ksum,pold,kold,dce,dke,dpe,phiscal,psum2
      real pold2,kold2,dce2,dke2,dpe2,nn,pmass,gpscal

!-----arrays
   
      real pprnt(irun+2),zprnt(irun)
      real pk(ng)
      real coor(6,np),ct(6,np)

      complex egtldetr(ng,ng,ng)

!-----ct:coor
!-----distribute/aligns

!hpf$ distribute zprnt(block)
!hpf$ distribute pprnt(block)
!hpf$ distribute pk(block)
!hpf$ distribute coor(*,block)
!hpf$ align (*,:) with coor(*,:) :: ct
!hpf$ distribute egtldetr(*,*,block)

!-----read zprnt (arrray contains z values for outputs)

      open(33,file="zprnt",form="formatted")
      do ii=1,irun
      read(33,*)zprnt(ii)
      enddo
      close(33)

      pp=ain**alpha
      pfin=afin**alpha
      prange=pfin-pp
      pprnt(1)=pp
      pprnt(irun+2)=pfin

      forall(ii=1:irun)pprnt(ii+1)=(1.0+zprnt(ii))**(-1.0/alpha)

      nsum=0
      ptime=pp

!-----outer loop

      write(6,*)'before outer loop in integz'

      do ii=2,irun+2
      delp=pprnt(ii)-pprnt(ii-1)
      nsloc=nint((delp/prange)*nsteps)
      epsloc=delp/nsloc
      nsum=nsum+nsloc
      ptime=ptime+epsloc*nsloc

!-----2nd order integrator:
      
      if(norder.eq.2)then

      tau1=epsloc
      tau2=0.5*epsloc

!-----inner loop

      do jj=1,nsloc

!-----coor->IN, ct->OUT
      call map1(alpha,omegatot,tau2,pp,coor,ct,np,ng)
!-----first time update
      pp=pp+tau2

!-----ct->IN, coor->OUT
      call map2(alpha,omegatot,tau1,pp,ct,coor,np,ng,egtldetr,gpscal)
!-----second time update
      pp=pp+tau2

!-----coor->IN, ct->OUT
      call map1(alpha,omegatot,tau2,pp,coor,ct,np,ng)
      coor=ct

!-----Cosmic energy test

      if(entest.eq.0)then

      call cosmic(coor,np,ng,egtldetr,psum,ksum,psum2,gpscal)
     
      phiscal=1.5*omegatot/pp**(1.0/alpha)
 
      psum=0.5*phiscal*psum
      psum2=phiscal*psum2

      write(222,*)pp,ksum,psum,psum2

      psum=psum*pp**(1.0/alpha)

      dke=(ksum-kold)/epsloc
      dpe=(psum-pold)/epsloc
      dce=dke+dpe*pp**(1.0/alpha)

      dke2=(3.0*ksum-4.0*kold+kold2)/(2.0*epsloc)
      dpe2=(3.0*psum-4.0*pold+pold2)/(2.0*epsloc)
      dce2=dke2+dpe2*pp**(1.0/alpha)

      pold2=pold
      pold=psum
      kold2=kold
      kold=ksum

      write(200,*)pp,dke,dpe*pp**(1.0/alpha),dce
      write(201,*)pp,dke2,dpe2*pp**(1.0/alpha),dce2

      endif

      enddo

!-----Print-out of particles and velocities at specified zs

      write(6,*)'time:',' a(t)=',pp**(1.0/alpha),
     #          ' z=',1.0/(pp**(1.0/alpha))-1.0,'ii=',nsloc
      if(frm.eq.1.or.frm.eq.2)then
      call prnt((ii-1)*7,coor,np,rL,ng,skip,
     #          frm,pp,alpha,pbchoice,
     #          omegatot,hubble)
!      call power(500+(ii-1)*7,coor,np,ng,rL,gpscal)
!      call powerphi(600+(ii-1)*7,coor,np,ng,rL,gpscal,egtldetr)

      endif

!-----Print-out 2-dimensional density field at specified zs

      if(frm.eq.3.or.frm.eq.4)then

      nn=(1.0*ng2d/(1.0*ng))

      call dens2d((ii-1)*7,coor,np,ng,hubble,pmass,ng2d,nn,rL,frm)

      endif

      if(ii.eq.(irun+2))then
      write(6,*)'final output:','time:',' a(t)=',pp**(1.0/alpha),
     #          ' z=',1.0/(pp**(1.0/alpha))-1.0,'ii=',nsum
      call prnt(-100,coor,np,rL,ng,skip,
     #          2,pp,alpha,2,omegatot,hubble)
      endif

      endif
      enddo

!      write(6,*)'nsteps=',nsteps, 'nsum=',nsum
!      write(6,*)'pprnt(irun+2)=',pprnt(irun+2), 'ptime=',ptime       
 
      return
      end subroutine
