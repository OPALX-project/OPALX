      subroutine ntrfas(fxg,fyg,fzg,fx,fy,fz,np,ng,ab,de,gh,
     #indx,jndx,kndx,indxp1,jndxp1,kndxp1)

!-----CIC Interpolation onto particles (SLOW)

      implicit none

!-----scalars

      integer np,ng,n

!-----arrays

      integer indx(np),jndx(np),kndx(np)
      integer indxp1(np),jndxp1(np),kndxp1(np)

      real fx(np),fy(np),fz(np),ab(np),de(np),gh(np)
      real fxg(ng,ng,ng),fyg(ng,ng,ng),fzg(ng,ng,ng)

!-----indx,jndx,kndx,indxp1,jndxp1,kndxp1,fy,fz,ab,de,gh:fx
!-----fyg,fzg:fxg
!-----distribute/aligns

!hpf$ distribute fx(block)
!hpf$ align (:) with fx(:) :: indx,jndx,kndx,indxp1,jndxp1,kndxp1
!hpf$ align (:) with fx(:) :: fy,fz,ab,de,gh

!hpf$ distribute fxg(*,*,block)
!hpf$ align (*,*,:) with fxg(*,*,:) :: fyg,fzg

      forall(n=1:np)fx(n)=
     # fxg(indx(n),jndx(n),kndx(n))*ab(n)*de(n)*gh(n)
     #+fxg(indx(n),jndxp1(n),kndx(n))*ab(n)*(1.-de(n))*gh(n)
     #+fxg(indx(n),jndxp1(n),kndxp1(n))*ab(n)*(1.-de(n))*(1.-gh(n))
     #+fxg(indx(n),jndx(n),kndxp1(n))*ab(n)*de(n)*(1.-gh(n))
     #+fxg(indxp1(n),jndx(n),kndxp1(n))*(1.-ab(n))*de(n)*(1.-gh(n))
     #+fxg(indxp1(n),
     #jndxp1(n),kndxp1(n))*(1.-ab(n))*(1.-de(n))*(1.-gh(n))
     #+fxg(indxp1(n),jndxp1(n),kndx(n))*(1.-ab(n))*(1.-de(n))*gh(n)
     #+fxg(indxp1(n),jndx(n),kndx(n))*(1.-ab(n))*de(n)*gh(n)
 
      forall(n=1:np)fy(n)=
     # fyg(indx(n),jndx(n),kndx(n))*ab(n)*de(n)*gh(n)
     #+fyg(indx(n),jndxp1(n),kndx(n))*ab(n)*(1.-de(n))*gh(n)
     #+fyg(indx(n),jndxp1(n),kndxp1(n))*ab(n)*(1.-de(n))*(1.-gh(n))
     #+fyg(indx(n),jndx(n),kndxp1(n))*ab(n)*de(n)*(1.-gh(n))
     #+fyg(indxp1(n),jndx(n),kndxp1(n))*(1.-ab(n))*de(n)*(1.-gh(n))
     #+fyg(indxp1(n),
     #jndxp1(n),kndxp1(n))*(1.-ab(n))*(1.-de(n))*(1.-gh(n))
     #+fyg(indxp1(n),jndxp1(n),kndx(n))*(1.-ab(n))*(1.-de(n))*gh(n)
     #+fyg(indxp1(n),jndx(n),kndx(n))*(1.-ab(n))*de(n)*gh(n)
 
      forall(n=1:np)fz(n)=
     # fzg(indx(n),jndx(n),kndx(n))*ab(n)*de(n)*gh(n)
     #+fzg(indx(n),jndxp1(n),kndx(n))*ab(n)*(1.-de(n))*gh(n)
     #+fzg(indx(n),jndxp1(n),kndxp1(n))*ab(n)*(1.-de(n))*(1.-gh(n))
     #+fzg(indx(n),jndx(n),kndxp1(n))*ab(n)*de(n)*(1.-gh(n))
     #+fzg(indxp1(n),jndx(n),kndxp1(n))*(1.-ab(n))*de(n)*(1.-gh(n))
     #+fzg(indxp1(n),
     #jndxp1(n),kndxp1(n))*(1.-ab(n))*(1.-de(n))*(1.-gh(n))
     #+fzg(indxp1(n),jndxp1(n),kndx(n))*(1.-ab(n))*(1.-de(n))*gh(n)
     #+fzg(indxp1(n),jndx(n),kndx(n))*(1.-ab(n))*de(n)*gh(n)

      return
      end subroutine
