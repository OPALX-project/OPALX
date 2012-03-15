      subroutine psp1cic(phi,ng,pkf)

!-----This routine computes the 1-D power spectrum of a 3-D input
!-----field (given in real space). 

      use hpf_library
     
      implicit none

      integer ii,jj,mm,ksign,icpy,ndx,ng

      real scal,psp,pk,pkf,nk,phi,tmpp,tpi,kk,mult

      complex erho,erhotr

      logical msk(ng,ng,ng)

      dimension psp(ng,ng,ng),tmpp(ng,ng,ng),ndx(ng,ng,ng)
      dimension erho(ng,ng,ng),erhotr(ng,ng,ng),phi(ng,ng,ng)
      dimension pk(ng*2+1),nk(ng*2+1),pkf(ng/2),kk(ng),mult(ng)

!hpf$ distribute pk(block)
!hpf$ distribute psp(*,*,block)
!hpf$ distribute nk(block)
!hpf$ distribute pkf(block)
!hpf$ distribute phi(*,*,block)
!hpf$ align (*,*,:) with phi(*,*,:) :: ndx,erho,tmpp
!hpf$ distribute erhotr(*,*,block)
!hpf$ distribute kk(*)
!hpf$ align (:) with kk(:) :: mult

      pkf=0.0

      erho=phi     
      tpi=8.0*atan(1.0)
      msk=.true.

!-----Forward fft. Because icpy=0, output goes to erhotr.

      ksign=1
      scal=1.0
      icpy=0

      call fft3d(ng,ng,ng,ksign,scal,icpy,erho,erhotr)

!-----3-D anti-CIC filter for deconvolution

      forall(ii=1:ng)kk(ii)=(ii-1)*tpi/(1.0*ng)

      do ii=1,ng
      if(ii.ge.ng/2+1)kk(ii)=(ii-ng-1)*tpi/(1.0*ng)
      enddo

      mult(1)=1.0

      forall(ii=2:ng)mult(ii)=
     #      1.0/(sin(kk(ii)/2.0)/(kk(ii)/2.0))**(2.0)

      forall(ii=1:ng,jj=1:ng,mm=1:ng)erhotr(ii,jj,mm)=
     #      mult(ii)*mult(jj)*mult(mm)*erhotr(ii,jj,mm)

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
      tmpp=1
      nk=sum_scatter(tmpp,nk,ndx,msk)

      do ii=1,ng/2
      pkf(ii)=pk(ii+1)/nk(ii+1)
      end do

      pkf=pkf/(1.0*ng**3)

      return
      end subroutine
