!----This is the main body of subroutine calls (just in case there
!----is any doubt!)

      subroutine solve(ng,np,nsteps,norder,iprint,iseed,jinit,
     #                norm,trans,nxint,nxdint,ns,zin,hubble,
     #                omegadm,rL,alpha,deut,qq,ss8,initp,
     #                skip,frm,pbchoice,zfin,zpr,irun,icbc,pmass,
     #                entest,ng2d,nginit)

      implicit none

!-----scalars

      integer ng,np,nsteps,norder,iprint,iseed,jinit,norm,trans
      integer nxint,nxdint,nbisec,nyint,ll,istp,initp,nx,nq
      integer skip,frm,pbchoice,irun,icbc,entest,ng2d,nginit

      real ns,zin,zfin,zpr,hubble,omegadm,tol,zmin,zmax,dzero
      real rL,alpha,deut,qq,ss8,pp,ain,dplus,ddot,omegab,omegatot
      real afin,apr,pmass,gpscal

!-----arrrays

      real psp(nginit/2+1,nginit/2+1,nginit/2+1)
      real kk(nginit/2+1,nginit/2+1,nginit/2+1)
      real coor(6,np)
      real gradx(nginit,nginit,nginit),grady(nginit,nginit,nginit)
      real gradz(nginit,nginit,nginit)
      complex egtldetr(ng,ng,ng)
      complex egtldetrinit(nginit,nginit,nginit)

!-----distribute/aligns
!-----psp:kk
!-----grady,gradz:gradx

!hpf$ distribute kk(*,*,block)
!hpf$ align (*,*,:) with kk(*,*,:) :: psp
!hpf$ distribute coor(*,block)
!hpf$ distribute gradx(*,*,block)
!hpf$ align (*,*,:) with gradx(*,*,:) :: grady,gradz
!hpf$ distribute egtldetr(*,*,block)
!hpf$ distribute egtldetrinit(*,*,block)

!-----Define nq, gpscal

      nq=nginit/2 + 1
      gpscal=(1.0*ng)/(1.0*nginit)

      write(6,*)'enter solve'
      
!-----Generate periodic Green's function (needed for generation of 
!-----initial conditions and time-stepping, now later!)

      if(jinit.eq.0)then
      call greenf(nginit,egtldetrinit)

!-----Generate initial power spectrum psp and carry on with the
!-----particle initialization, if jinit=0 (default)

      call initspec(nginit,norm,trans,nxint,nxdint,
     #              ns,zin,hubble,omegadm,rL,deut,
     #              qq,ss8,psp,kk,dplus,dzero,ddot)
    
      kk(1,1,1)=0.0

      write(6,*)'done with initspec'

!-----Generate the gradients needed for zeldo with density power 
!-----spectrum psp from initspec. Also produce various power spectra
!-----for diagnostic purposes. 

      call indens(nginit,nq,rL,psp,gradx,grady,gradz,hubble,iseed,
     #            pmass,frm,pbchoice,initp,egtldetrinit)

      write(6,*) 'Initial density field generated'
      
!-----Initial Zeldovich move

      call zeldo(dplus,dzero,ddot,gradx,grady,gradz,nginit,np,rL,zin,
     #           coor,skip,frm,initp,pbchoice,gpscal)
      
      write(6,*) 'Zeldovich step completed'
      
!-----Scale initial conditions to the chosen grid size
 
      coor=gpscal*coor

      endif

!-----Read in externally specified initial conditions, jinit=1

      if(jinit.eq.1)then

      call icread(nginit,np,rL,zin,coor,frm,pbchoice,icbc,gpscal)

!-----Scale initial conditions to the chosen grid size
 
      coor=gpscal*coor

      endif

!-----Read in externally specified intermediate time-output, jinit=2

      if(jinit.eq.2)then

      call cont(ng,np,rL,coor)

      endif

!-----Convert to `symplectic velocities' at a=a_in (x coords already OK)
!-----and define afin and apr

      ain=1.0/(1.0+zin)
      afin=1.0/(1.0+zfin)
      apr=1.0/(1.0+zpr)

      coor(2,:)=coor(2,:)*ain*ain
      coor(4,:)=coor(4,:)*ain*ain
      coor(6,:)=coor(6,:)*ain*ain

!-----Introduce Omega_tot

      omegab=deut/hubble/hubble
      omegatot=omegadm + omegab

!-----calculate initial power spectrum

      if(initp.eq.0.or.initp.eq.1)then
      call power(800,coor,np,ng,rL,gpscal)
      endif

!-----PM time-stepping
!-----free-run from ain to afin with printing after apr set by iprint

      call greenf(ng,egtldetr)

      if(irun.eq.0)then
      call integ(coor,ain,afin,apr,alpha,nsteps,norder,np,ng,iprint,
     #           egtldetr,omegatot,rL,skip,frm,pbchoice,hubble,
     #           entest,gpscal)

!-----run with z array (of size irun) controlling output history

      else
      call integz(coor,ain,afin,alpha,nsteps,norder,np,ng,irun,
     #           egtldetr,omegatot,rL,skip,frm,pbchoice,hubble,
     #           entest,pmass,ng2d,gpscal)
      endif

      return
      end subroutine

