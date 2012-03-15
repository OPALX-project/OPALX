      subroutine grvpot(coor,egtldetr,fx,fy,fz,np,ng,gpscal)

!-----Computes gravitational forces on particles

      implicit none

!-----scalars

      integer np,ng,ksign,icpy
      integer left,right,dim,ii

      real scal,gpscal

!-----arrays

      integer indx(np),jndx(np),kndx(np)
      integer indxp1(np),jndxp1(np),kndxp1(np)

      real fx(np),fy(np),fz(np),ab(np),de(np),gh(np)
      real coor(6,np)
      real rho(ng,ng,ng),rhol(ng,ng,ng),rhor(ng,ng,ng)  
      real fxg(ng,ng,ng),fyg(ng,ng,ng),fzg(ng,ng,ng)

      complex erho(ng,ng,ng),erhotr(ng,ng,ng),egtldetr(ng,ng,ng)

!-----indx,jndx,kndx,indxp1,jndxp1,kndxp1:coor
!-----fx,fy,fz,ab,de,gh:coor
!-----fxg,fyg,fzg,erho,rhol,rhor:rho
!-----erhotr:egtldetr

!-----align/distributes

!hpf$ distribute coor(*,block)
!hpf$ align (:) with coor(*,:) :: indx,jndx,kndx,indxp1,jndxp1,kndxp1
!hpf$ align (:) with coor(*,:) :: fx,fy,fz,ab,de,gh
!hpf$ distribute rho(*,*,block)
!hpf$ align (*,*,:) with rho(*,*,:) :: fxg,fyg,fzg,erho,rhol,rhor
!hpf$ distribute egtldetr(*,*,block)
!hpf$ align (*,*,:) with egtldetr(*,*,:) :: erhotr

!-----Deposit particles on the grid using area weighting (CIC):
  
!      write(6,*)'entering grvpot'

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

      call fft3d(ng,ng,ng,ksign,scal,icpy,erhotr,erho)

!-----Compensate for inverse FFT scaling

      erho=erho/(1.0*ng*ng*ng)

!-----rho is really phi

      rho=real(erho,8)

!-----f=-grad(phi)

!      fxg=0.5*(cshift(rho,-1,1)-cshift(rho,1,1))
!      fyg=0.5*(cshift(rho,-1,2)-cshift(rho,1,2))
!      fzg=0.5*(cshift(rho,-1,3)-cshift(rho,1,3))

!-----cshift "done by hand" (cshift not working on qsc)

      right=1
      left=1

      dim=1
      call shift(rho,rhol,rhor,left,right,ng,dim)
      fxg=0.5*(rhor-rhol)

      dim=2
      call shift(rho,rhol,rhor,left,right,ng,dim)
      fyg=0.5*(rhor-rhol)

      dim=3
      call shift(rho,rhol,rhor,left,right,ng,dim)
      fzg=0.5*(rhor-rhol)


!-----Interpolate the forces onto the particles:

      call ntrfas(fxg,fyg,fzg,fx,fy,fz,np,ng,ab,de,gh,
     #indx,jndx,kndx,indxp1,jndxp1,kndxp1)

!      write(6,*)'after ntrfas'

      return
      end subroutine
