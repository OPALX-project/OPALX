      subroutine integ(coor,ain,afin,apr,alpha,nsteps,norder,np,
     #                 ng,iprint,egtldetr,omegatot,rL,skip,frm,
     #                 pbchoice,hubble,entest,gpscal)

!-----PM time-stepper

      implicit none

!-----scalars

      integer nsteps,norder,np,ng,iprint,ii,skip,frm,pbchoice,entest

      real tau1,tau2,tau3,tau4,eps,alf,pp,alpha,omegatot,rL
      real ain,afin,apr,ppr,pfin,hubble
      real psum,ksum,pold,kold,dce,dke,dpe,phiscal,psum2
      real pold2,kold2,dce2,dke2,dpe2,gpscal

!-----arrays
   
      real pk(ng)
      real coor(6,np),ct(6,np)

      complex egtldetr(ng,ng,ng)

!-----ct:coor
!-----distribute/aligns

!hpf$ distribute pk(block)
!hpf$ distribute coor(*,block)
!hpf$ align (*,:) with coor(*,:) :: ct
!hpf$ distribute egtldetr(*,*,block)

      pp=ain**alpha
      pfin=afin**alpha
      ppr=apr**alpha

      write(6,*)'time:',' a(t)=',pp**(1.0/alpha),
     #          ' z=',1.0/(pp**(1.0/alpha))-1.0

!-----compute timestep

      eps=(pfin-pp)/(1.0*nsteps)

!-----2nd order integrator:
      
      if(norder.eq.2)then
      
      tau1=eps
      tau2=0.5*eps

      do ii=1,nsteps

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

      dke=(ksum-kold)/eps
      dpe=(psum-pold)/eps
      dce=dke+dpe*pp**(1.0/alpha)

      dke2=(3.0*ksum-4.0*kold+kold2)/(2.0*eps)
      dpe2=(3.0*psum-4.0*pold+pold2)/(2.0*eps)
      dce2=dke2+dpe2*pp**(1.0/alpha)

      pold2=pold
      pold=psum
      kold2=kold
      kold=ksum

      write(200,*)pp,dke,dpe*pp**(1.0/alpha),dce
      write(201,*)pp,dke2,dpe2*pp**(1.0/alpha),dce2

      endif

!-----Print positions and velocities

      if(mod(ii,iprint).eq.0)then

      call prnt(ii,coor,np,rL,ng,skip,frm,pp,alpha,pbchoice,
     #                                omegatot,hubble)
      call power(500+(ii-1)*7,coor,np,ng,rL,gpscal)
      
      endif

      write(6,*)'time:',' a(t)=',pp**(1.0/alpha),
     #          ' z=',1.0/(pp**(1.0/alpha))-1.0,'ii=',ii

      enddo

      endif
 
      return
      end subroutine
