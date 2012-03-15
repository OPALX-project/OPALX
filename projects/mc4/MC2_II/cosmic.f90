!-----This routine computes the on-particle potential and kinetic
!-----energies and checks for conservation of the total cosmic energy.

      subroutine cosmic(coor,np,ng,egtldetr,psum,ksum,psum2,gpscal)

      implicit none

      integer np,ng,ii

      real psum,ksum,psum2,gpscal

      real coor(6,np),phi(np),phig(ng,ng,ng)

      complex egtldetr(ng,ng,ng)

!hpf$ distribute coor(*,block)
!hpf$ align (:) with coor(*,:) :: phi
!hpf$ distribute egtldetr(*,*,block)
!hpf$ align (*,*,:) with egtldetr(*,*,:) :: phig

!-----Get on-particle potential from grvphi
	
      call grvphi(coor,egtldetr,phi,np,ng,phig,gpscal)

      psum=sum(phi)

      psum2=sum(phig)

!-----Compute on-particle kinetic energy

      ksum=0.5*sum(coor(2,:)**2+coor(4,:)**2+coor(6,:)**2)

      return
      end subroutine
