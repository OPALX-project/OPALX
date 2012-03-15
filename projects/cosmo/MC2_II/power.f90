      subroutine power(ifile,coor,np,ng,rL,gpscal)

!-----Reads in particles and positions provided from outside
    
      implicit none

      integer ng,np,mm,ii,ifile

      integer indx(np),jndx(np),kndx(np)
      integer indxp1(np),jndxp1(np),kndxp1(np)

      real tpi,tpiL
      real rL,rscale,rLb,gpscal
      real xscal,vscal,vfac,pfac

      real coor(6,np)
      real pkf(ng/2),pkfc(ng/2)
      real rho(ng,ng,ng)
      real ab(np),de(np),gh(np)

!hpf$ distribute pkf(block)
!hpf$ align (:) with pkf(:) :: pkfc 
!hpf$ distribute coor(*,block)
!hpf$ align (:) with coor(*,:) :: indx,jndx,kndx,ab,de,gh
!hpf$ align (:) with coor(*,:) :: indxp1,jndxp1,kndxp1
!hpf$ distribute rho(*,*,block)

      tpi=8.0*atan(1.0)  ! 2 pi
      rLb=1.0*ng
      write(6,*)'entering power'

!-----Compute CIC density field from particle positions 
!-----Normalized to particle number density.

      call rhoord(coor,rho,np,ng,ab,de,gh,
     #      indx,jndx,kndx,indxp1,jndxp1,kndxp1,gpscal)

!-----Compute 1-D power spectrum (physical units)

      rho=(rho-1.0)

      call pspec1(rho,ng,pkf)

      tpiL=tpi/rL
      rscale=rL/(1.0*ng)

      do mm=1,ng/2
      write(ifile,*) mm*tpiL,pkf(mm)*rscale**3
      enddo

!-----Compute 1-D power spectrum with CIC `de-convolution'

      call psp1cic(rho,ng,pkfc)

      do mm=1,ng/2
      write(ifile+1,*)mm*tpiL,pkfc(mm)*rscale**3
      enddo

      return
      end subroutine
