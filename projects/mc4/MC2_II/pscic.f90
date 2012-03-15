      subroutine pscic(phi,ng,pkf)

!-----This routine computes the 1-D power spectrum of a 3-D input
!-----field (given in real space), Fourier transformed and then 
!-----convolved with the 3-D CIC filter. 

      use hpf_library
     
      implicit none

!-----scalars
 
      integer ii,jj,mm,ksign,icpy,ng

      real scal,hh,tpi

!-----arrays

      logical msk(ng,ng,ng)

      integer ndx(ng,ng,ng)

      real pkf(ng/2)
      real mult(ng),kk(ng),kknew(ng)
      real pk(ng*2+1),nk(ng*2+1)
      real psp(ng,ng,ng),phi(ng,ng,ng),rho(ng,ng,ng)

      complex erho(ng,ng,ng),erhotr(ng,ng,ng)

!-----mult:kk
!-----nk:pk
!-----msk,ndx,phi,erho:rho
!-----psp:erhotr
!-----distribute/aligns

!hpf$ distribute pkf(block)
!hpf$ distribute kk(block)
!hpf$ align (:) with kk(:) :: mult,kknew
!hpf$ distribute pk(block)
!hpf$ align (:) with pk(:) :: nk
!hpf$ distribute rho(*,*,block)
!hpf$ align (*,*,:) with rho(*,*,:) :: msk,ndx,phi,erho
!hpf$ distribute erhotr(*,*,block)
!hpf$ align (*,*,:) with erhotr(*,*,:) :: psp

      msk=.true.
      erho=phi     
      tpi=8.0*atan(1.0)

!-----Forward fft. because icpy=0, output goes to erhotr.

      ksign=1
      scal=1.0
      icpy=0

      call fft3d(ng,ng,ng,ksign,scal,icpy,erho,erhotr)

!-----3-D CIC filter. Box units.

      forall(ii=1:ng)kk(ii)=(ii-1)*tpi/(1.0*ng)
     
      do ii=1,ng
      if(ii.ge.ng/2+1)kk(ii)=(ii-ng-1)*tpi/(1.0*ng)
      enddo

      do ii=1,ng/2

      jj=ii+ng/2
      mm=jj-ng

      kknew(ii)=tpi*(ii-1)/(1.0*ng)
      kknew(jj)=tpi*(mm-1)/(1.0*ng)      

      enddo

!      do ii=1,ng
!      write(77,*)ii*tpi,kk(ii),kknew(ii)
!      enddo

      mult(1)=1.0

      forall(ii=2:ng)mult(ii)=
     #     (sin(kk(ii)/2.0)/(kk(ii)/2.0))**2.0
                          
      forall(ii=1:ng,jj=1:ng,mm=1:ng)erhotr(ii,jj,mm)=
     #     mult(ii)*mult(jj)*mult(mm)*erhotr(ii,jj,mm)


!-----3-D power spectrum.

      psp=real(erhotr*conjg(erhotr),8)

!-----Bin to get 1-D spectrum

      pk=0.0
      nk=0.0

      forall(ii=1:ng/2,jj=1:ng/2,mm=1:ng/2)ndx(ii,jj,mm)=
     #  1+nint(sqrt(real((ii-1)**2+(jj-1)**2+(mm-1)**2)))
      forall(ii=ng/2+1:ng,jj=1:ng/2,mm=1:ng/2)ndx(ii,jj,mm)=
     #  1+nint(sqrt(real((ii-ng-1)**2+(jj-1)**2+(mm-1)**2)))
      forall(ii=1:ng/2,jj=ng/2+1:ng,mm=1:ng/2)ndx(ii,jj,mm)=
     #  1+nint(sqrt(real((ii-1)**2+(jj-ng-1)**2+(mm-1)**2)))
      forall(ii=1:ng/2,jj=1:ng/2,mm=ng/2+1:ng)ndx(ii,jj,mm)=
     #  1+nint(sqrt(real((ii-1)**2+(jj-1)**2+(mm-ng-1)**2)))
      forall(ii=1:ng/2,jj=ng/2+1:ng,mm=ng/2+1:ng)ndx(ii,jj,mm)=
     #  1+nint(sqrt(real((ii-1)**2+(jj-ng-1)**2+(mm-ng-1)**2)))
      forall(ii=ng/2+1:ng,jj=ng/2+1:ng,mm=1:ng/2)ndx(ii,jj,mm)=
     #  1+nint(sqrt(real((ii-ng-1)**2+(jj-ng-1)**2+(mm-1)**2)))
      forall(ii=ng/2+1:ng,jj=1:ng/2,mm=ng/2+1:ng)ndx(ii,jj,mm)=
     #  1+nint(sqrt(real((ii-ng-1)**2+(jj-1)**2+(mm-ng-1)**2)))
      forall(ii=ng/2+1:ng,jj=ng/2+1:ng,mm=ng/2+1:ng)ndx(ii,jj,mm)=
     #  1+nint(sqrt(real((ii-ng-1)**2+(jj-ng-1)**2+(mm-ng-1)**2)))

      pk=sum_scatter(psp,pk,ndx,msk)
      rho=1
      nk=sum_scatter(rho,nk,ndx,msk)

      do ii=1,ng/2
      pkf(ii)=pk(ii+1)/nk(ii+1)
      end do

      pkf=pkf/(1.0*ng)**3

      return
      end subroutine
