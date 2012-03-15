      subroutine powerphi(ifile,coor,np,ng,rL,gpscal,egtldetr)

!-----Reads in particles and positions provided from outside
    
      implicit none

      integer ng,np,mm,ii,ifile,ksign,icpy

      integer indx(np),jndx(np),kndx(np)
      integer indxp1(np),jndxp1(np),kndxp1(np)

      real tpi,tpiL,scal
      real rL,rscale,rLb,gpscal
      real xscal,vscal,vfac,pfac

      real coor(6,np)
      real pkf(ng/2),pkfc(ng/2)
      real rho(ng,ng,ng)
      real ab(np),de(np),gh(np)

      complex erho(ng,ng,ng),erhotr(ng,ng,ng),egtldetr(ng,ng,ng)

!hpf$ distribute pkf(block)
!hpf$ align (:) with pkf(:) :: pkfc 
!hpf$ distribute coor(*,block)
!hpf$ align (:) with coor(*,:) :: indx,jndx,kndx,ab,de,gh
!hpf$ align (:) with coor(*,:) :: indxp1,jndxp1,kndxp1
!hpf$ distribute rho(*,*,block)
!hpf$ align (*,*,:) with rho(*,*,:) :: erho
!hpf$ distribute egtldetr(*,*,block)
!hpf$ align (*,*,:) with egtldetr(*,*,:) :: erhotr

      tpi=8.0*atan(1.0)  ! 2 pi
      rLb=1.0*ng
      write(6,*)'entering power-phi'

!-----Compute CIC density field from particle positions 
!-----Normalized to particle number density.

      call rhoord(coor,rho,np,ng,ab,de,gh,
     #      indx,jndx,kndx,indxp1,jndxp1,kndxp1,gpscal)

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

      call fft3d(ng,ng,ng,ksign,scal,icpy,erhotr,erho)

!-----Compensate for inverse FFT scaling

      erho=erho/(1.0*ng*ng*ng)

!-----rho is really phi

      rho=real(erho,8)

!-----Compute 1-D power spectrum (physical units)

      call pspec1(rho,ng,pkf)

      tpiL=tpi/rL
      rscale=rL/(1.0*ng)

      do mm=1,ng/2
      write(ifile,*) mm*tpiL,pkf(mm)*rscale**3
      enddo


      return
      end subroutine
