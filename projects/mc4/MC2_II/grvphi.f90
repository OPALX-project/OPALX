!-----This is a mildly truncated version of grvpot which solves
!-----the Poisson equation, finds the scalar potential, and then
!-----calls intscal to return the on-particle potential.

      subroutine grvphi(coor,egtldetr,phi,np,ng,phig,gpscal)

!-----Computes on-particle scalar potential

      implicit none

      integer np,ng,ksign,icpy

      integer indx(np),jndx(np),kndx(np)
      integer indxp1(np),jndxp1(np),kndxp1(np)

      real scal,gpscal

      real coor(6,np),phi(np),ab(np),de(np),gh(np)
      real phig(ng,ng,ng),rho(ng,ng,ng)

      complex erho(ng,ng,ng),erhotr(ng,ng,ng)
      complex egtldetr(ng,ng,ng)

!hpf$ distribute coor(*,block)
!hpf$ align (:) with coor(*,:) :: phi,indx,jndx,kndx,ab,de,gh
!hpf$ align (:) with coor(*,:) :: indxp1,jndxp1,kndxp1
!hpf$ distribute rho(*,*,block)
!hpf$ align (*,*,:) with rho(*,*,:) :: phig,erho
!hpf$ distribute egtldetr(*,*,block)
!hpf$ align (*,*,:) with egtldetr(*,*,:) :: erhotr

!-----Deposit particles on the grid using area weighting (CIC):

      call rhoord(coor,rho,np,ng,ab,de,gh,
     #indx,jndx,kndx,indxp1,jndxp1,kndxp1,gpscal)

      erho=rho

!-----Fourier transform the mass density:

      ksign=1
      scal=1.0
      icpy=0

      call fft3d(ng,ng,ng,ksign,scal,icpy,erho,erhotr)

!-----Convolve the mass density with the Green's function:

      erhotr=erhotr*egtldetr

!-----Inverse FFT

      ksign=-1

      call fft3d (ng,ng,ng,ksign,scal,icpy,erhotr,erho)

!-----Compensate for inverse FFT scaling

      erho=erho/(1.0*ng*ng*ng)

!-----rho is really phi

      rho=real(erho,8)
      phig=(rho)

!-----Interpolate the potential onto the particles:

      call intscal(rho,phi,np,ng,ab,de,gh,
     #indx,jndx,kndx,indxp1,jndxp1,kndxp1)

      write(333,*)sum(phig),sum(phi)

      return
      end subroutine
