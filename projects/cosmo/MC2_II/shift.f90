      subroutine shift(rho,rhol,rhor,left,right,ng,dim)

!-----This routine does exactly the same operation c-shift would normally
!-----do, just with foralls; necessary routine for running on QSC in parallel
!-----with the COMPAQ-compiler

!-----A three-dimensional array (usually rho or phi) comes from the
!-----outside into this routine, left and right gives the direction 
!-----(left=1, right=1, nothing happens otherwise) in which the array
!-----will be shifted

      implicit none

      integer left,right,ng,ii,dim
      
      real rho(ng,ng,ng),rhol(ng,ng,ng),rhor(ng,ng,ng),lambda(ng,ng,ng)

!hpf$ distribute rho(*,*,block)
!hpf$ align (*,*,:) with rho(*,*,:) :: lambda,rhol,rhor

!-----Initialize arrays

      rhol=0.0
      rhor=0.0

!-----Left shift

      if(left.eq.1)then

      if(dim.eq.1)forall(ii=1:ng-1)lambda(ii,:,:)=rho(ii+1,:,:)      
      if(dim.eq.2)forall(ii=1:ng-1)lambda(:,ii,:)=rho(:,ii+1,:)      
      if(dim.eq.3)forall(ii=1:ng-1)lambda(:,:,ii)=rho(:,:,ii+1)      

      if(dim.eq.1)lambda(ng,:,:)=rho(1,:,:)
      if(dim.eq.2)lambda(:,ng,:)=rho(:,1,:)
      if(dim.eq.3)lambda(:,:,ng)=rho(:,:,1)

      rhol=lambda

      endif

!-----Right shift

      if(right.eq.1)then

      if(dim.eq.1)forall(ii=2:ng)lambda(ii,:,:)=rho(ii-1,:,:)      
      if(dim.eq.2)forall(ii=2:ng)lambda(:,ii,:)=rho(:,ii-1,:)      
      if(dim.eq.3)forall(ii=2:ng)lambda(:,:,ii)=rho(:,:,ii-1)      

      if(dim.eq.1)lambda(1,:,:)=rho(ng,:,:)
      if(dim.eq.2)lambda(:,1,:)=rho(:,ng,:)
      if(dim.eq.3)lambda(:,:,1)=rho(:,:,ng)

      rhor=lambda

      endif
      
      return
      end subroutine	
