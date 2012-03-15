      subroutine map2(alpha,omegatot,tau,pp,coor,ct,np,ng,
     #                egtldetr,gpscal)

!-----Kick step

      implicit none

!-----scalars

      integer np,ng,ii

      real tau,fscal,pp,omegatot,phiscal,alpha,adot,gpscal

!-----arrays

      real fx(np),fy(np),fz(np)
      real coor(6,np),ct(6,np)

      complex egtldetr(ng,ng,ng)

!-----fx,fy,fz,ct:coor
!-----distribute/aligns

!hpf$ distribute coor(*,block)
!hpf$ align (:) with coor(*,:) :: fx,fy,fz
!hpf$ align (*,:) with coor(*,:) :: ct

!hpf$ distribute egtldetr(*,*,block)

      ct=0.0

!-----adot in units of H_0^(-1)

      adot=sqrt((omegatot+(1.0-omegatot)*pp**(3.0/alpha))
     #         /pp**(1.0/alpha))

      call grvpot(coor,egtldetr,fx,fy,fz,np,ng,gpscal)

!-----Compute scaling functions

      phiscal=1.5*omegatot/pp**(1.0/alpha)
      fscal=phiscal/(alpha*adot*pp**(1.0-1.0/alpha))

!      write(112,*)1./pp-1.,fscal

!-----Compute self-fields and kick


      ct(1,:)=coor(1,:)
      ct(2,:)=coor(2,:)+fx(:)*fscal*tau

      ct(3,:)=coor(3,:)
      ct(4,:)=coor(4,:)+fy(:)*fscal*tau

      ct(5,:)=coor(5,:)
      ct(6,:)=coor(6,:)+fz(:)*fscal*tau

      return
      end subroutine



