      subroutine intscal(phig,phi,np,ng,ab,de,gh,
     #indx,jndx,kndx,indxp1,jndxp1,kndxp1)

!-----CIC Interpolation onto particles (SLOW)
!-----Scalar field from grid (phig) to scalar field on particles (phi)

      implicit none

      integer np,ng,nn

      integer indx(np),jndx(np),kndx(np)
      integer indxp1(np),jndxp1(np),kndxp1(np)

      real phi(np),ab(np),de(np),gh(np)
      real phig(ng,ng,ng)

!hpf$ distribute phi(block)
!hpf$ align (:) with phi(:) :: ab,de,gh,indx,jndx,kndx
!hpf$ align (:) with phi(:) :: indxp1,jndxp1,kndxp1
!hpf$ distribute phig(*,*,block)

      forall(nn=1:np)phi(nn)=
     # phig(indx(nn),jndx(nn),kndx(nn))*ab(nn)*de(nn)*gh(nn)
     #+phig(indx(nn),jndxp1(nn),kndx(nn))*ab(nn)*(1.-de(nn))*gh(nn)
     #+phig(indx(nn),jndxp1(nn),kndxp1(nn))
     #     *ab(nn)*(1.-de(nn))*(1.-gh(nn))
     #+phig(indx(nn),jndx(nn),kndxp1(nn))*ab(nn)*de(nn)*(1.-gh(nn))
     #+phig(indxp1(nn),jndx(nn),kndxp1(nn))
     #     *(1.-ab(nn))*de(nn)*(1.-gh(nn))
     #+phig(indxp1(nn),
     #jndxp1(nn),kndxp1(nn))*(1.-ab(nn))*(1.-de(nn))*(1.-gh(nn))
     #+phig(indxp1(nn),jndxp1(nn),kndx(nn))
     #     *(1.-ab(nn))*(1.-de(nn))*gh(nn)
     #+phig(indxp1(nn),jndx(nn),kndx(nn))*(1.-ab(nn))*de(nn)*gh(nn)

      return
      end subroutine
