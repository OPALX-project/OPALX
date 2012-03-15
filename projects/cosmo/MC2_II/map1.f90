      subroutine map1(alpha,omegatot,tau,pp,coor,ct,np,ng)

!-----Stream step
!-----Historical note: map1 rewritten from CMF to HPF, Aug'97

      implicit none

!-----scalars

      integer np,ng,ii

      real tau,fn,pp,omegatot,adot,alpha,pre

!-----arrays

      real vec(np)
      real coor(6,np),ct(6,np)

!-----vec,ct:coor
!-----distribute/aligns

!hpf$ distribute coor(*,block)
!hpf$ align (:) with coor(*,:) :: vec
!hpf$ align (*,:) with coor(*,:) :: ct

      ct=0.0

!-----adot in units of H_0^(-1)

      adot=sqrt((omegatot+(1.0-omegatot)*pp**(3.0/alpha))
     #         /pp**(1.0/alpha))

!-----prefactor for stream

      pre=1.0/(alpha*adot*pp**(1.0+1.0/alpha))

!-----implement stream

      ct(1,:)=coor(1,:)+pre*tau*coor(2,:)

      ct(2,:)=coor(2,:)

      ct(3,:)=coor(3,:)+pre*tau*coor(4,:)

      ct(4,:)=coor(4,:)

      ct(5,:)=coor(5,:)+pre*tau*coor(6,:)

      ct(6,:)=coor(6,:)


!-----Enforce periodic BCs:

      fn=1.0*ng

      vec=ct(1,:)
      where(vec.lt.0.)vec=vec+fn
      where(vec.ge.fn)vec=vec-fn
      ct(1,:)=vec

      vec=ct(3,:)
      where(vec.lt.0.)vec=vec+fn
      where(vec.ge.fn)vec=vec-fn
      ct(3,:)=vec

      vec=ct(5,:)
      where(vec.lt.0.)vec=vec+fn
      where(vec.ge.fn)vec=vec-fn
      ct(5,:)=vec

      return
      end subroutine
