      subroutine rhoord2d(coor,rho,np,ng,nn)

!-----CIC Deposition (smart)
!-----This routine is a 2D-CIC scheme!

      use hpf_library

      implicit none

      integer np,ng

      logical msk(np)

      integer indx(np),jndx(np),indxp1(np)
      integer jndxp1(np)

      real nn

      real coor(2,np),vol(np),ab(np),de(np)
      real rho(ng,ng),tmp(ng,ng)  

!hpf$ distribute coor(*,block)
!hpf$ align (:) with coor(*,:) :: msk,vol
!hpf$ align (:) with coor(*,:) :: ab,de,indx,jndx
!hpf$ align (:) with coor(*,:) :: indxp1,jndxp1
!hpf$ distribute rho(*,block)
!hpf$ align (*,:) with rho(*,:) :: tmp

      msk=.true.

      indx=int(coor(1,:)) + 1
      jndx=int(coor(2,:)) + 1

      indxp1=indx+1
      where( indxp1.eq.ng+nn)indxp1=1
      jndxp1=jndx+1
      where( jndxp1.eq.ng+nn)jndxp1=1

      ab=indx-coor(1,:)
      de=jndx-coor(2,:)

      rho=0.0

!-----1 (i,j):

      vol=ab*de
      tmp=0.0
      tmp=sum_scatter(vol,tmp,indx,jndx,msk)
      rho=rho+tmp

!-----2 (i,j+1):

      vol=ab*(1.0-de)
      tmp=0.0
      tmp=sum_scatter(vol,tmp,indx,jndxp1,msk)
      rho=rho+tmp

!-----3 (i+1,j+1):

      vol=(1.0-ab)*(1.0-de)
      tmp=0.0
      tmp=sum_scatter(vol,tmp,indxp1,jndxp1,msk)
      rho=rho+tmp

!-----4 (i+1,j):

      vol=(1.0-ab)*de
      tmp=0.0
      tmp=sum_scatter(vol,tmp,indxp1,jndx,msk)
      rho=rho+tmp

      return
      end subroutine










