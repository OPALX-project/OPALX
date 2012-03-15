      subroutine bin1d(phi,ng,pkf)

!-----This routine computes bin a 3-D positive field into a
!-----1-D object. 

      use hpf_library
     
      implicit none

!-----scalars

      integer ii,jj,mm,ng

      real psp

!-----arrays

      logical msk(ng,ng,ng)

      integer ndx(ng,ng,ng)                  

      real pkf(ng/2),pk(ng*2+1),nk(ng*2+1)
      real phi(ng,ng,ng),tmpp(ng,ng,ng)      

!-----nk:pk
!-----msk,ndx,tmpp:phi
!-----distribute/align

!hpf$ distribute pkf(block)
!hpf$ distribute pk(block)
!hpf$ align (:) with pk(:) :: nk
!hpf$ distribute phi(*,*,block)
!hpf$ align (*,*,:) with phi(*,*,:) :: msk,ndx,tmpp

!-----Bin to get 1-D object 
 
      msk=.true.
      pk=0.0
      nk=0.0

      forall(ii=1:ng/2,jj=1:ng/2,mm=1:ng/2)ndx(ii,jj,mm)=     
     #   1+nint(sqrt(real((ii-1)**2+(jj-1)**2+(mm-1)**2)))
      forall(ii=ng/2+1:ng,jj=1:ng/2,mm=1:ng/2)ndx(ii,jj,mm)=  
     #   1+nint(sqrt(real((ii-ng-1)**2+(jj-1)**2+(mm-1)**2)))
      forall(ii=1:ng/2,jj=ng/2+1:ng,mm=1:ng/2)ndx(ii,jj,mm)=  
     #   1+nint(sqrt(real((ii-1)**2+(jj-ng-1)**2+(mm-1)**2)))
      forall(ii=1:ng/2,jj=1:ng/2,mm=ng/2+1:ng)ndx(ii,jj,mm)= 
     #   1+nint(sqrt(real((ii-1)**2+(jj-1)**2+(mm-ng-1)**2)))
      forall(ii=1:ng/2,jj=ng/2+1:ng,mm=ng/2+1:ng)ndx(ii,jj,mm)=  
     #   1+nint(sqrt(real((ii-1)**2+(jj-ng-1)**2+(mm-ng-1)**2)))
      forall(ii=ng/2+1:ng,jj=ng/2+1:ng,mm=1:ng/2)ndx(ii,jj,mm)= 
     #   1+nint(sqrt(real((ii-ng-1)**2+(jj-ng-1)**2+(mm-1)**2)))
      forall(ii=ng/2+1:ng,jj=1:ng/2,mm=ng/2+1:ng)ndx(ii,jj,mm)=
     #   1+nint(sqrt(real((ii-ng-1)**2+(jj-1)**2+(mm-ng-1)**2)))
      forall(ii=ng/2+1:ng,jj=ng/2+1:ng,mm=ng/2+1:ng)ndx(ii,jj,mm)= 
     #   1+nint(sqrt(real((ii-ng-1)**2+(jj-ng-1)**2+(mm-ng-1)**2)))

      pk=sum_scatter(phi,pk,ndx,msk)
      tmpp=1.0
      nk=sum_scatter(tmpp,nk,ndx,msk)

      do ii=1,ng/2
      pkf(ii)=pk(ii+1)/nk(ii+1)
      enddo

      pkf=pkf/(1.0*ng)**3

      return
      end subroutine








