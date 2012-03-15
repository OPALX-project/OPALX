!-----Generate the initial matter power spectrum

      subroutine initspec(ng,norm,trans,nxint,nxdint,ns,zin,
     #                    hubble,omegadm,rL,deut,qq,ss8,psp,
     #                    kk,dplus,dzero,ddot)

      implicit none

!-----scalars

      integer ng,norm,trans,nxint,nxdint,nq,ii,jj,mm

      real ns,zin,hubble,omegadm,omegatot,pnorm,dplus,ddot,dzero
      real rL,rLb,alpha,deut,qq,ss8,tpi,tt,tpiL,scal

!-----arrays

      real psp(ng/2+1,ng/2+1,ng/2+1),kk(ng/2+1,ng/2+1,ng/2+1)
      real pspz(ng/2+1,ng/2+1,ng/2+1)

!-----kk,pspz:psp
!distribute/aligns
!hpf$ distribute psp(*,*,block)
!hpf$ align (*,*,:) with psp(*,*,:) :: kk,pspz

!-----Define constants

      tpi=8.0*atan(1.0)    ! 2 pi
      tt=2.728             ! COBE/FIRAS CMB temperature in K
      rLb=ng*1.0           ! box size as real number
      tpiL=tpi/rL	   ! k0, physical units
      nq=ng/2+1            ! kgrid box size

!-----Define omegatot

      omegatot=omegadm+deut/hubble/hubble

!-----Define kk

      kk=0.0

      do ii=1,nq
      do jj=1,nq
      forall(mm=1:nq) kk(ii,jj,mm)=(tpiL)*sqrt(((ii-1.0)**2
     #                       +(jj-1.0)**2+(mm-1.0)**2))		
      enddo
      enddo


!-----COBE Normalization

      call cobe(omegatot,qq,hubble,nxint,pnorm,ns,norm,tt,dzero)

!-----Pure power law with spectral index ns (note that Peebles COBE 
!-----normalization is valid for ns close to unity but the Bunn/White
!-----fits are more general. When using the sigma8 normalization, it
!-----doesn't matter what ns is.

      if (trans.eq.0) then	
      where(kk.gt.0.0)psp=pnorm*kk**ns
      where(kk.le.0.0)psp=0.0
      psp(1,1,1) = 0.0
      endif

!-----Generate Lambda-CDM power spectra for a flat universe using 
!-----Hu/Sugiyama or Klypin/Holtzman transfer functions.

      if (trans.ne.0) then

      call spectrum(omegatot,deut,ns,qq,tt,nq,hubble,pnorm,
     #              ss8,trans,psp,kk,norm,rL,omegadm)
      endif

!-----Scale by growth factor to get linear power spectrum at initial z

      call initscale(hubble,omegatot,zin,nxdint,scal,dzero,ddot)

      dplus=scal
      write(6,*)'(dplus/dzero)**2=',(dplus/dzero)**2
      pspz=psp*(dplus/dzero)**2

!-----Write 1-D power spectrum for diagnostic purposes

      do jj=1,nq
      do mm=1,nq
      write(21,*)kk(1,jj,mm),psp(1,jj,mm),pspz(1,jj,mm)
      enddo
      enddo

      return
      end subroutine

!-----Auxiliary routines below:

!-----Lambda-CDM power spectra for a flat universe using Hu/Sugiyama
!-----or Klypin/Holtzman transfer functions.

      subroutine spectrum(omegatot,deut,ns,qq,tt,nq,hubble,pnorm,
     #                    ss8,trans,psp,kk,norm,rL,omegadm)

      implicit none

!-----scalars

      integer ii,nk,trans,nq,norm,nmax

      real omegatot,omegab,ns,nn,qq,tt,cc,hubble,aa,bb,pnorm
      real ss8,pi,ss8husug,deut,sig8,newnorm,sigscal
      real alpha,akh1,akh2,s8in,ss8kh,ss8kh2,rL,omegadm

!-----arrays

      real kk(nq,nq,nq),qhusug(nq,nq,nq),qkh(nq,nq,nq)
      real fac(nq,nq,nq),trkh(nq,nq,nq),trhusug(nq,nq,nq)
      real trcmbf(nq,nq,nq)
      real delhusug(nq,nq,nq),delkh(nq,nq,nq),psp(nq,nq,nq)

!-----kk,qhusug,qkh,fac,trkh,trhusug,delhusug,delkh:psp
!-----distribute/aligns

!hpf$ distribute psp(*,*,block)
!hpf$ align (*,*,:) with psp(*,*,:) :: kk,qhusug,qkh,fac,trkh,trcmbf
!hpf$ align (*,*,:) with psp(*,*,:) :: trhusug,delhusug,delkh

!-----Define constants

      pi=4.0*atan(1.0) !  pi
      fac=kk*kk*kk/(2.0*pi*pi)
      omegab=deut/hubble/hubble

!-----COBE-normalized Hu/Sugiyama power spectrum

      if(trans.eq.1)then

      aa=((46.9*omegatot*hubble*hubble)**0.67)
     #   *(1.0 + (32.1*omegatot*hubble*hubble)**(-0.532))
      bb=((12.0*omegatot*hubble*hubble)**0.424)
     #   *(1.0 + (45.0*omegatot*hubble*hubble)**(-0.582))

      qhusug=((kk*(tt/2.7)**2)/(omegatot*hubble))
     #        /sqrt(aa**(-omegab/omegatot)*bb**(-(omegab/omegatot)**3))

!-----Cheap fix to avoid 0/0 problem at kk=0

      where(qhusug.eq.0.0)qhusug=1.0e-8

      trhusug=((log(1.0+2.34*qhusug))/(2.34*qhusug))/(1.0+3.89*qhusug
     #       + (16.1*qhusug)**2 + (5.46*qhusug)**3
     #       + (6.71*qhusug)**4)**(0.25)

      psp=(trhusug**2)*pnorm*kk**ns
      delhusug=fac*psp

      endif

!-----COBE-normalized Klypin/Holtzman power spectrum

      if(trans.eq.2)then

      akh1=((46.9*omegatot*hubble*hubble)**(0.67))
     #     *(1.0+(32.1*omegatot*hubble*hubble)**(-0.532))
      akh2=((12.0*omegatot*hubble*hubble)**(0.424))
     #     *(1.0+(45.0*omegatot*hubble*hubble)**(-0.582))
      alpha=(akh1**(-omegab/omegatot))*(akh2**((-omegab/omegatot)**3))
      qkh=(kk*(tt/2.7)**2)/(omegatot*hubble*(sqrt(alpha))*
     #    (1.0-omegab/omegatot)**(0.6))

!-----Cheap fix to avoid 0/0 problem at kk=0

      where(qkh.eq.0.0)qkh=1.0e-8

      trkh=((log(1.0+2.34*qkh))/(2.34*qkh))*(1.0+13.0*qkh+(10.5*qkh)**2
     #     +(10.4*qkh)**3+(6.51*qkh)**4)**(-0.25)

      psp=(trkh**2)*pnorm*kk**ns
      delkh=fac*psp

      endif

!-----COBE-normalized CMBFAST power spectrum

      if(trans.eq.3)then

      open(7,file="incmb",status='old')
      read(7,*)nmax                      ! number of lines in CMBFAST input file
      close(7)

      call fitlin3d(nmax,nq,rL,trcmbf,omegatot,omegab,omegadm)      
!      call cmbfit(nmax,nq,rL,trcmbf,omegatot,omegab,omegadm)

      psp=(trcmbf**2)*pnorm*kk**ns

      endif


!-----Compute sigma_8 using 1-D normalized distributions

      nk=10001

      call sigma(sig8,pnorm,trans,hubble,omegatot,deut,nk,ns,omegadm)

!-----If desired, normalize to target value of sigma_8, ss8

      write(6,*)'sig8=',sig8,'trans=',trans,'norm=',norm

      if(norm.eq.1)then
      sigscal=ss8*ss8/sig8/sig8
      psp=psp*sigscal
      endif

      newnorm=sigscal*pnorm

      call sigma(sig8,newnorm,trans,hubble,omegatot,deut,nk,ns,omegadm)

      write(6,*)'new sig8=',sig8,'target sig8=',ss8

      return
      end subroutine

!-----COBE Normalization routine

      subroutine cobe(omegatot,qq,hubble,nxint,pnorm,ns,norm,
     #                tt,dzero)

      implicit none

!-----scalars

      integer nxint,ii,trans,norm

      real tt,qq,cc,hubble,pi,dy,omegatot,pnorm,sumsmp
      real abw,bbw,vv,dzero,deltah,ns

!-----arrays

      real integrand(nxint),yy(nxint)

!-----integrand:yy
!-----distribute/aligns

!hpf$ distribute yy(block)
!hpf$ align (:) with yy(:) :: integrand

!-----Define constants

      pi=4.0*atan(1.0)  ! pi
      cc=299792.4583    ! speed of light in km/s

!-----Peebles 1984 normalization 

      dy=1.0/(nxint-1)
      forall(ii=1:nxint) yy(ii)=(ii - 1)*dy
      forall(ii=1:nxint) integrand(ii)=(1.0 + (1.0/omegatot - 1.0)
     #                                 *yy(ii)**(6.0/5.0))**(-1.5)

      call smpint(nxint,integrand,dy,sumsmp)

      dzero=sumsmp/sqrt(omegatot)

      if(norm.eq.2)then
      pnorm=(96.0/5.0)*(cc**4/((100.0*hubble)**4))*(qq*qq/(tt*tt))
     #     *(dzero*dzero/omegatot/omegatot)*pi*pi
      endif

!-----Bunn and White fit for COBE 4-yr data

      if(norm.le.1)then
      abw=-0.95
      bbw=-0.169
      vv=(abw*(ns-1.0)+bbw*(ns-1.0))**2
      deltah=1.94e-5*omegatot**(-0.785-0.05*log(omegatot))
     #       *exp(vv)
      pnorm=deltah**2*2.0*pi**2*cc**4/(hubble*100.0)**4
      endif

      pnorm=pnorm*hubble**4

      return
      end subroutine

!-----Growth factor calculation for Lambda-CDM

      subroutine initscale(hubble,omegatot,zin,nxdint,
     #                    scal,dzero,ddot)

      implicit none

!-----scalars

      integer nxdint,ii

      real dy,sumsmp,zin,scal,dzero,ddot,hubble
      real term1,term2,term3,hzin,zdot,ddz,omegatot

!-----arrrays

      real, dimension (nxdint) :: integrand,yy

!-----integrand:yy
!-----distribute/aligns

!hpf$ distribute yy(block)
!hpf$ align (:) with yy(:) :: integrand

!-----Compute dplus (aka scal)

      dy=1.0/(nxdint-1)
      forall(ii=1:nxdint) yy(ii)=(ii - 1)*dy

      integrand=(1.0 + ((1.0/omegatot-1.0)/(1.0+zin)**3)
     #                 *yy**(6.0/5.0))**(-1.5)

      call smpint(nxdint,integrand,dy,sumsmp)

      scal=(1.0/(1.0+zin))
     #      *sqrt((1.0/omegatot-1.0)/(1.0+zin)**3+1.0)*sumsmp

!-----Compute the time derivative of Dplus, ddot

      term1=-scal/(1.0+zin)

      term2=-(1.5/((1.0/omegatot - 1.0)/(1.0+zin)**3 + 1.0))
     #     *((1.0/omegatot - 1.0)/(1.0+zin)**4)*scal

      integrand=(yy**(6.0/5.0))
     #          *(1.0 + ((1.0/omegatot - 1.0)/(1.0+zin)**3)
     #          *yy**(6.0/5.0))**(-2.5)

      call smpint(nxdint,integrand,dy,sumsmp)

      term3=4.5*((1.0/omegatot - 1.0)/(1.0+zin)**5)
     #     *sqrt((1.0/omegatot - 1.0)/(1.0+zin)**3+1.0)*sumsmp 

      ddz=term1+term2+term3

!-----Time is measured in units of H_0^(-1), thus no multiplication by
!-----H_0 of the RHS for hzin

      hzin=sqrt(1.0+omegatot*((zin+1.0)**3-1.0))

      zdot=-(zin+1.0)*hzin

      ddot=zdot*ddz

      return
      end subroutine

!-----Sigma_8 from 1-D power spectra 

      subroutine sigma(sig8,pnorm,trans,hubble,omegatot,deut,nk,ns,
     #                 omegadm)

      implicit none

!-----scalars

      integer ii,nk,trans,nmax

      real omegatot,omegab,ns,tt,dk,kmax,kmin,hubble,aa,bb,pnorm
      real pi,alpha,akh1,akh2,sig8,rr,deut,sumsmp,omegadm

!-----arrays

      real qhusug(nk),trhusug(nk),qkh(nk),trkh(nk),psptmp(nk)
      real rk(nk),besselj(nk),win(nk),integrand(nk),kk(nk),trcmbf(nk)

!-----qhusug,trhusug,qkh,trkh,psptmp,rk,besselj,win,integrand:kk
!-----distribute/aligns

!hpf$ distribute kk(block)
!hpf$ align(:) with kk(:) :: qhusug,trhusug,qkh,trkh,psptmp,trcmbf
!hpf$ align(:) with kk(:) :: rk,besselj,win,integrand

!-----Define constants

      kmax=10.0         ! max k for sigma_8 integration
      kmin=0.001        ! min k for sigma_8 integration
      pi=4.0*atan(1.0)  ! pi
      tt=2.728          ! COBE/FIRAS CMB temperature in K

      omegab=deut/hubble/hubble

      dk=(kmax - kmin)/(nk-1.0)
      forall(ii=1:nk) kk(ii)=kmin + (ii-1)*dk

!-----Hu/Sugiyama

      if(trans.eq.1)then

      aa=((46.9*omegatot*hubble*hubble)**0.67)
     #   *(1.0 + (32.1*omegatot*hubble*hubble)**(-0.532))
      bb=((12.0*omegatot*hubble*hubble)**0.424)
     #   *(1.0 + (45.0*omegatot*hubble*hubble)**(-0.582))

      qhusug=((kk*(tt/2.7)**2)/(omegatot*hubble))
     #        /sqrt(aa**(-omegab/omegatot)*bb**(-(omegab/omegatot)**3))

      trhusug=((log(1.0+2.34*qhusug))/(2.34*qhusug))/(1.0+3.89*qhusug
     #       + (16.1*qhusug)**2 + (5.46*qhusug)**3
     #       + (6.71*qhusug)**4)**(0.25)

      psptmp=(trhusug**2)*pnorm*kk**ns

      do ii=1,nk
      write(30,*)kk(ii),psptmp(ii)
      enddo

      endif

!----Klypin/Holtzman

      if(trans.eq.2)then

      akh1=((46.9*omegatot*hubble*hubble)**(0.67))
     #     *(1.0+(32.1*omegatot*hubble*hubble)**(-0.532))
      akh2=((12.0*omegatot*hubble*hubble)**(0.424))
     #     *(1.0+(45.0*omegatot*hubble*hubble)**(-0.582))
      alpha=(akh1**(-omegab/omegatot))*(akh2**((-omegab/omegatot)**3))
      qkh=(kk*(tt/2.7)**2)/(omegatot*hubble*(sqrt(alpha))*
     #    (1.0-omegab/omegatot)**(0.6))

      trkh=((log(1.0+2.34*qkh))/(2.34*qkh))*(1.0+13.0*qkh+(10.5*qkh)**2
     #     +(10.4*qkh)**3+(6.51*qkh)**4)**(-0.25)

      psptmp=(trkh**2)*pnorm*kk**ns

      do ii=1,nk
      write(31,*)kk(ii),psptmp(ii)
      enddo

      endif

!-----COBE-normalized CMBFAST power spectrum

      if(trans.eq.3)then

      open(7,file="incmb",status='old')
      read(7,*)nmax                      ! number of lines in CMBFAST input file
      close(7)

!      call cmbfit1d(nmax,nk,kmax,kmin,trcmbf,omegatot,
!     #              omegab,omegadm)

      call fitlin1d(nmax,nk,kmax,kmin,trcmbf,omegatot,
     #              omegab,omegadm)

      psptmp=(trcmbf**2)*pnorm*kk**ns

      do ii=1,nk
      write(31,*)kk(ii),psptmp(ii),trcmbf(ii)
      enddo

      endif


!-----Sigma_8 integration over top hat window function in k-space

      rr=8.0
      rk=rr*kk
      besselj=(sin(rk))/rk/rk - (cos(rk))/rk
      win=9.0*besselj*besselj/rk/rk
      integrand=kk*kk*psptmp*win

      call smpint(nk,integrand,dk,sumsmp)

      sig8=sqrt(sumsmp/(2.0*pi*pi))

      return
      end subroutine

!-----Simpson integrator

      subroutine smpint(nstep,integrand,dx,sumsmp)

      implicit none

!-----scalars

      integer nstep,ii

      real dx,sumsmp

!-----arrays

      real, dimension(nstep) :: integrand,smp

!-----smp:integrand
!-----distribute/aligns

!hpf$ distribute integrand(block)
!hpf$ align(:) with integrand(:) :: smp

      smp=4.0/3.0
      forall(ii=1:nstep:2)smp(ii)=2.0/3.0
      smp(1)=1.0/3.0
      smp(nstep)=1.0/3.0

      sumsmp=dx*sum(smp*integrand)

      return
      end subroutine
