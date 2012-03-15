      subroutine rhoord(coor,rho,np,ng,ab,de,gh,
     #      indx,jndx,kndx,indxp1,jndxp1,kndxp1,gpscal)

!-----CIC Deposition (Basic, slow routine)

      use hpf_library

      implicit none

!-----scalars

      integer np,ng

      real gpscal

!-----arrrays

      logical msk(np)

      integer indx(np),jndx(np),kndx(np)
      integer indxp1(np),jndxp1(np),kndxp1(np)

      real vol(np),ab(np),de(np),gh(np)
      real coor(6,np)
      real rho(ng,ng,ng),tmp(ng,ng,ng)  

!-----msk,indx,jndx,kndx,indxp1,jndxp1,kndxp1,vol,ab,de,gh:coor
!-----tmp:rho
!-----distribute/aligns

!hpf$ distribute coor(*,block)
!hpf$ align (:) with coor(*,:) :: msk
!hpf$ align (:) with coor(*,:) :: indx,jndx,kndx,indxp1,jndxp1,kndxp1
!hpf$ align (:) with coor(*,:) :: vol,ab,de,gh
!hpf$ distribute rho(*,*,block)
!hpf$ align (*,*,:) with rho(*,*,:) :: tmp

      msk=.true.

      indx=int(coor(1,:)) + 1
      jndx=int(coor(3,:)) + 1
      kndx=int(coor(5,:)) + 1

      indxp1=indx+1
      where( indxp1.eq.ng+1)indxp1=1
      jndxp1=jndx+1
      where( jndxp1.eq.ng+1)jndxp1=1
      kndxp1=kndx+1
      where( kndxp1.eq.ng+1)kndxp1=1

      ab=indx-coor(1,:)
      de=jndx-coor(3,:)
      gh=kndx-coor(5,:)

      rho=0.0

!-----1 (i,j,k):

      vol=ab*de*gh
      tmp=0.0
      tmp=sum_scatter(vol,tmp,indx,jndx,kndx,msk)
      rho=rho+tmp

!-----2 (i,j+1,k):

      vol=ab*(1.0-de)*gh
      tmp=0.0
      tmp=sum_scatter(vol,tmp,indx,jndxp1,kndx,msk)
      rho=rho+tmp

!-----3 (i,j+1,k+1):

      vol=ab*(1.0-de)*(1.0-gh)
      tmp=0.0
      tmp=sum_scatter(vol,tmp,indx,jndxp1,kndxp1,msk)
      rho=rho+tmp

!-----4 (i,j,k+1):

      vol=ab*de*(1.0-gh)
      tmp=0.0
      tmp=sum_scatter(vol,tmp,indx,jndx,kndxp1,msk)
      rho=rho+tmp

!-----5 (i+1,j,k+1):

      vol=(1.0-ab)*de*(1.0-gh)
      tmp=0.0
      tmp=sum_scatter(vol,tmp,indxp1,jndx,kndxp1,msk)
      rho=rho+tmp

!-----6 (i+1,j+1,k+1):

      vol=(1.0-ab)*(1.0-de)*(1.0-gh)
      tmp=0.0
      tmp=sum_scatter(vol,tmp,indxp1,jndxp1,kndxp1,msk)
      rho=rho+tmp

!-----7 (i+1,j+1,k):

      vol=(1.0-ab)*(1.0-de)*gh
      tmp=0.0
      tmp=sum_scatter(vol,tmp,indxp1,jndxp1,kndx,msk)
      rho=rho+tmp

!-----8 (i+1,j,k):

      vol=(1.0-ab)*de*gh
      tmp=0.0
      tmp=sum_scatter(vol,tmp,indxp1,jndx,kndx,msk)
      rho=rho+tmp

      rho=gpscal**3*rho
      
      return
      end subroutine










