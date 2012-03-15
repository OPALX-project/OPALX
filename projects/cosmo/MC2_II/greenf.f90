      subroutine greenf(ng,egtldetr)

!-----Computes periodic Green's function for del^2 phi=delta for
!-----the centered-difference Laplacian, i.e., second-order accuracy
!-----for a smooth density field.

      implicit none

!-----scalars

      integer ng,ii,jj,mm

      real tpin

!-----arrays

      complex egtlde(ng,ng,ng),egtldetr(ng,ng,ng)

!-----distribute/align

!hpf$ distribute egtlde(*,*,block)
!hpf$ distribute egtldetr(*,*,block)

      tpin=4.*asin(1.0)/(1.0*ng)

      do ii=1,ng
      forall(jj=1:ng,mm=1:ng,ii+jj+mm.ne.3)egtlde(ii,jj,mm)=
     #  0.5/(cos(tpin*(ii-1))+cos(tpin*(jj-1))+cos(tpin*(mm-1))-3.0)
      enddo

      egtlde(1,1,1)=0.0

!-----Transposing the Green's function

      do ii=1,ng
!hpf$ independent, new(mm)
      do jj=1,ng
!hpf$ independent
      do mm=1,ng
      egtldetr(ii,jj,mm)=egtlde(mm,ii,jj)
      enddo
      enddo
      enddo

      return
      end
