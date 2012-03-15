!-----Generates the gradients needed for zeldo for a single realization
!-----of an input matter power spectrum from initspec. For this
!-----realization, it also produces the potential power spectrum and the 
!-----density field power spectrum.

      subroutine indens(ng,nq,rL,psp,gradxv,gradyv,gradzv,hubble,iseed,
     #                  pmass,frm,pbchoice,initp,egtldetr)

      implicit none

!-----scalars

      integer ii,jj,mm,ng,nq,iseed,np
      integer left,right,dim      
      integer seed_size,clock,ksign,icpy
      integer frm,pbchoice,initp

      real rL,tpiL,rLb
      real rscale,tpi,ns,w,scal,hubble,pmass

      complex imi

!-----arrays

      integer, dimension(:), allocatable :: seed

      real pk(ng/2)
      real psp(nq,nq,nq),aran(nq,nq,nq),bran(nq,nq,nq)
      real amp(nq,nq,nq),ramp1(nq,nq,nq),ramp2(nq,nq,nq)
      real dens(ng,ng,ng),phi(ng,ng,ng),phil(ng,ng,ng),phir(ng,ng,ng)
      real gradx(ng,ng,ng),grady(ng,ng,ng),gradz(ng,ng,ng)
      real gradxv(ng,ng,ng),gradyv(ng,ng,ng),gradzv(ng,ng,ng)

      complex ee(ng)
      complex aarg(nq,nq,nq)
      complex egtldetr(ng,ng,ng),dc(ng,ng,ng),dctr(ng,ng,ng)
      complex kdc(ng,ng,ng)

!-----aran,bran,amp,ramp1,ramp2,aarg:psp
!-----phi,phil,phir,gradx,grady,gradz:dens
!-----align/distributes

!hpf$ distribute ee(block)
!hpf$ distribute pk(block)
!hpf$ distribute psp(*,*,block)
!hpf$ align (*,*,:) with psp(*,*,:) :: aran,bran,amp,ramp1,ramp2
!hpf$ align (*,*,:) with psp(*,*,:) :: aarg
!hpf$ distribute dens(*,*,block)
!hpf$ align (*,*,:) with dens(*,*,:) :: phi,phil,phir,gradx,grady,gradz
!hpf$ align (*,*,:) with dens(*,*,:) :: gradxv,gradyv,gradzv
!hpf$ distribute dc(*,*,block)
!hpf$ align (*,*,:) with dc(*,*,:) :: kdc
!hpf$ distribute dctr(*,*,block)
!hpf$ distribute egtldetr(*,*,block)

!-----Define constants

      tpi=8.0*atan(1.0) ! 2 pi
      rLb=ng*1.0        ! Box length converted to real length
      rscale=rL/rLb     ! scaling from box to physical units
      tpiL=tpi/rL       ! k0, physical units

!-----Density power spectrum amplitude

      amp=0.0
      amp=sqrt(psp)
      amp(1,1,1)=0.0

!-----Begin generation of realization

      aran=0.0
      bran=0.0

      call random_seed(SIZE=seed_size)
      allocate (seed(seed_size))
      CALL random_seed(GET=seed(1:seed_size))
      CALL system_clock(COUNT=clock)
      do ii=1,seed_size	
      	 seed(ii)=iseed+(ii-1)*1000
      enddo
      CALL random_seed(PUT=seed(1:seed_size))

!-----First set of Rayleigh/uniform deviates

      CALL random_number(aran)
      CALL random_number(bran)

      ramp1=amp*sqrt(-log(aran))
      aarg=ramp1*cmplx(cos(tpi*bran),sin(tpi*bran))

!-----Enforce reality condition at corners

      aarg(nq,1,1)=ramp1(nq,1,1)
      aarg(1,nq,1)=ramp1(1,nq,1)
      aarg(1,1,nq)=ramp1(1,1,nq)
      aarg(nq,nq,1)=ramp1(nq,nq,1)
      aarg(1,nq,nq)=ramp1(1,nq,nq)
      aarg(nq,1,nq)=ramp1(nq,1,nq)
      aarg(nq,nq,nq)=ramp1(nq,nq,nq)

      do ii=1,nq
      do jj=1,nq
         forall(mm=1:nq)dc(ii,jj,mm)=aarg(ii,jj,mm)
	 forall(mm=1:nq)dc(mod(ng+1-ii,ng)+1,mod(ng+1-jj,ng)+1,
     #   mod(ng+1-mm,ng)+1)=conjg(dc(ii,jj,mm))
      enddo
      enddo

!-----Second set of Rayleigh/uniform deviates

      CALL random_number(aran)
      CALL random_number(bran)

      ramp1=amp*sqrt(-log(aran))
      aarg=ramp1*cmplx(cos(tpi*bran),sin(tpi*bran))

!-----Enforce reality condition at corners

      aarg(nq,1,1)=ramp1(nq,1,1)
      aarg(1,nq,1)=ramp1(1,nq,1)
      aarg(1,1,nq)=ramp1(1,1,nq)
      aarg(nq,nq,1)=ramp1(nq,nq,1)
      aarg(1,nq,nq)=ramp1(1,nq,nq)
      aarg(nq,1,nq)=ramp1(nq,1,nq)
      aarg(nq,nq,nq)=ramp1(nq,nq,nq)

      do ii=1,nq
      do jj=1,nq
	 forall(mm=1:nq)dc(ii,jj,mod(ng+1-mm,ng)+1)=aarg(ii,jj,mm)
	 forall(mm=1:nq)dc(mod(ng+1-ii,ng)+1,mod(ng+1-jj,ng)+1,mm)= 
     #   conjg(dc(ii,jj,mod(ng+1-mm,ng)+1))
      enddo
      enddo

!-----Third set of Rayleigh/uniform deviates

      CALL random_number(aran)
      CALL random_number(bran)

      ramp1=amp*sqrt(-log(aran))
      aarg=ramp1*cmplx(cos(tpi*bran),sin(tpi*bran))

!-----Enforce reality condition at corners

      aarg(nq,1,1)=ramp1(nq,1,1)
      aarg(1,nq,1)=ramp1(1,nq,1)
      aarg(1,1,nq)=ramp1(1,1,nq)
      aarg(nq,nq,1)=ramp1(nq,nq,1)
      aarg(1,nq,nq)=ramp1(1,nq,nq)
      aarg(nq,1,nq)=ramp1(nq,1,nq)
      aarg(nq,nq,nq)=ramp1(nq,nq,nq)

      do ii=1,nq
      do jj=1,nq
	 forall(mm=1:nq)dc(ii,mod(ng+1-jj,ng)+1,mm)=aarg(ii,jj,mm)
	 forall(mm=1:nq)dc(mod(ng+1-ii,ng)+1,jj,mod(ng+1-mm,ng)+1)
     #               =conjg(dc(ii,mod(ng+1-jj,ng)+1,mm))
      enddo
      enddo

!-----Fourth set of Rayleigh/uniform deviates  
     
      CALL random_number(aran)
      CALL random_number(bran)

      ramp1=amp*sqrt(-log(aran))
      aarg=ramp1*cmplx(cos(tpi*bran),sin(tpi*bran))

!-----Enforce reality condition at corners

      aarg(nq,1,1)=ramp1(nq,1,1)
      aarg(1,nq,1)=ramp1(1,nq,1)
      aarg(1,1,nq)=ramp1(1,1,nq)
      aarg(nq,nq,1)=ramp1(nq,nq,1)
      aarg(1,nq,nq)=ramp1(1,nq,nq)
      aarg(nq,1,nq)=ramp1(nq,1,nq)
      aarg(nq,nq,nq)=ramp1(nq,nq,nq)

      do ii=1,nq
      do jj=1,nq
         forall(mm=1:nq)dc(ii,mod(ng+1-jj,ng)+1,
     #        mod(ng+1-mm,ng)+1)=aarg(ii,jj,mm)
         forall(mm=1:nq)dc(mod(ng+1-ii,ng)+1,jj,mm)=
     #      conjg(dc(ii,mod(ng+1-jj,ng)+1,mod(ng+1-mm,ng)+1))
      enddo
      enddo
      
!-----Inverse FFT to get dens(x) (icpy=1 means that the result 
!-----of fft3d, dctr, gets put back in dc)

      ksign=-1
      scal=1.0
      icpy=1

      call fft3d(ng,ng,ng,ksign,scal,icpy,dc,dctr)

!-----Define dens and scale by volume^3/2 factor to get true density

      dens=real(dc,8)/rL**1.5

      write(6,*) 'maxdens=', maxval(dens)

!-----NOTE: Conventions are such that there is no compensation above for 
!-----inverse FFT scaling!

!-----Why rscale**3? This is because the density field needs to
!-----be multiplied by rL**1.5, hence the power spectrum by rL**3, to 
!-----return to the original definition of P(k). At the same time, the
!-----average over the total number of points requires a division by ng**3
!-----just as for computing the Gaussian checks above. Consequently, we
!-----have a multiplicative factor of rL**3/ng**3, i.e., rscale**3.

      do ii=1,ng
      write(39,*) ii,dens(ng/2,ng/2,ii)
      enddo

!-----Write the initial density field

      if(initp.eq.2)then

!-----In box units

      if(pbchoice.eq.1)then

!-----In ASCII

      if(frm.eq.1)then

      open(11,form='formatted',status='unknown')
      
      do ii=1,ng
      do jj=1,ng
      do mm=1,ng
      write(11,*)dens(ii,jj,mm)
      enddo
      enddo
      enddo

      close(11)

      endif

!-----In Binary

      if(frm.eq.2)then

      open(11,form='unformatted',access='sequential',
     #     recordtype='stream') 

      write(11)dens(:,:,:)
      
      close(11)

      endif

      endif

!-----In physical units

      if(pbchoice.eq.2)then

!-----In ASCII

      if(frm.eq.1)then

      open(11,form='formatted',status='unknown')
      
      do ii=1,ng
      do jj=1,ng
      do mm=1,ng
      write(11,*)dens(ii,jj,mm)*pmass*ng**3/(rL**3)
      enddo
      enddo
      enddo

      close(11)

      endif

!-----In Binary

      if(frm.eq.2)then

      open(11,form='unformatted',access='sequential',
     #     recordtype='stream') 

      write(11)dens(:,:,:)*pmass*ng**3/(rL**3)
      
      close(11)

      endif

      endif

      endif

!-----End of writing initial density field

      call pspec1(dens,ng,pk)
      do mm=1,ng/2
      write(41,*) mm*tpiL,(pk(mm))*rscale**3
      enddo

      

!-----Output CIC-Filtered density power spectrum (physical units)

      call pscic(dens,ng,pk)

      do mm=1,ng/2
      write(42,*)mm*tpiL,(pk(mm))*rscale**3
      enddo

!-----Solve Poisson equation for `bare' phi
!-----Step I: Fourier transform the mass density
!-----instead of icpy=0, icpy=1 means that the result of fft3d, 
!-----dctr, gets put back in dc. This involves an extra transpose.

!-----Reuse dc and dctr as temporary arrays

      dc=dens

      ksign=1
      scal=1.0
      icpy=0

      call fft3d(ng,ng,ng,ksign,scal,icpy,dc,dctr)

!-----Step II: Convolve the mass density with the Green's function

      dctr=dctr*egtldetr

!-----Step III: Carry out inverse FFT

      ksign=-1

      call fft3d (ng,ng,ng,ksign,scal,icpy,dctr,dc)

!-----Compensate for inverse FFT scaling

      dc=dc/(1.0*ng*ng*ng)

!-----Step IV: Final result

      phi=real(dc,8)

!-----Get power spectrum of bare potential by calling pspec1

      call pspec1(phi,ng,pk)

      do mm=1,ng/2
      write(43,*) mm*tpiL,(pk(mm))*rscale**3,
     #            pk(mm)*(mm*tpiL)**4*rscale**7
      enddo

!-----Compute centered gradients in real space from phi

!      gradx=0.5*(cshift(phi,1,1)-cshift(phi,-1,1))
!      grady=0.5*(cshift(phi,1,2)-cshift(phi,-1,2))
!      gradz=0.5*(cshift(phi,1,3)-cshift(phi,-1,3))

!-----cshift "done by hand" (cshift not working on qsc)

      right=1
      left=1

      dim=1
      call shift(phi,phil,phir,left,right,ng,dim)
      gradx=0.5*(phil-phir)

      dim=2
      call shift(phi,phil,phir,left,right,ng,dim)
      grady=0.5*(phil-phir)

      dim=3
      call shift(phi,phil,phir,left,right,ng,dim)
      gradz=0.5*(phil-phir)

!-----Compute the gradient of the potential via FFT, first Fourier 
!-----transform phi(x)

      dc=phi

      ksign=1
      scal=1.0
      icpy=1

      call fft3d(ng,ng,ng,ksign,scal,icpy,dc,dctr)

!-----icpy=1 means that the result of fft3d, dctr, gets put back in dc.
!-----This involves an extra transpose. (Historical remark: Set icpy=1 
!-----for now to make the code behave like the original cmf code.)

!-----This piece of code produces the correct k-vector in the usual 
!-----screwed-up FFT format (can use mod instead of do's but doesn't 
!-----really matter)

      imi=(0.0,1.0)
      
      do ii=1,ng/2

      jj=ii + ng/2
      mm=jj - ng
      w=tpi*(ii-1)/(1.0*ng)
      ee(ii)=-imi*w
      w=tpi*(mm-1)/(1.0*ng)
      ee(jj) = -imi*w

      enddo

      

!-----End of little piece of code.

!-----Compute x gradient from k-space potential field.

      forall(ii=1:ng)kdc(ii,:,:)=ee(ii)*dc(ii,:,:)

!-----Go back to real space.

      icpy=1
      ksign=-1
      scal=1.0

      call fft3d (ng,ng,ng,ksign,scal,icpy,kdc,dctr)

!-----Compensation for inverse FFT scaling

      kdc=kdc/(rLb**3)

      gradxv=real(kdc)

!-----Now compute y gradient from original k-space potential field:

      forall(ii=1:ng)kdc(:,ii,:)=ee(ii)*dc(:,ii,:)

!-----Go back to real space.

      icpy=1
      ksign=-1
      scal=1.0

      CALL fft3d (ng,ng,ng,ksign,scal,icpy,kdc,dctr)

!-----Compensation for inverse FFT scaling

      kdc=kdc/(rLb**3)

      gradyv=real(kdc)

!-----Finally, compute z gradient from original k-space potential 
!-----field:

      forall(ii=1:ng)kdc(:,:,ii)=ee(ii)*dc(:,:,ii)

!-----Go back to real space.

      icpy=1
      ksign=-1
      scal=1.0

      CALL fft3d (ng,ng,ng,ksign,scal,icpy,kdc,dctr)

!-----Compensation for inverse FFT scaling

      kdc=kdc/(rLb**3)

      gradzv=real(kdc)
      
      do ii=1,ng
      write(60,*)ii,gradz(ii,1,1),gradzv(ii,1,1)
      enddo

      return
      end subroutine








