!-----HPF 3D FFT routine

      subroutine fft3d(n1,n2,n3,ksign,scale,icpy,x,x3)

      implicit none

!-----scalars

      integer, intent(in) :: n1, n2, n3, ksign, icpy
      real, intent(in) :: scale

!-----arrays

      complex, dimension(n1,n2,n3), intent(inout) :: x
      complex, dimension(n3,n2,n1), intent(out) :: x3

!-----distribute/aligns

!hpf$ distribute x(*,*,block)
!hpf$ distribute x3(*,*,block)

!-----Local variables
!-----scalars

      integer :: i,j,k

!-----arrays

      complex, dimension(n2,n1,n3) :: x2

!-----distribute/aligns

!hpf$ distribute x2(*,*,block)

!-----interface mfft_local1

      interface
      extrinsic (hpf_local) subroutine mfft_local1 (ksign,x)
      implicit none
      integer, intent(in) :: ksign
      complex, dimension(:,:,:), intent(inout) :: x
!hpf$ distribute(*,*,block) :: x
      end subroutine mfft_local1
      end interface

!-----interface mfft_local2

      interface
      extrinsic (hpf_local) subroutine mfft_local2 (ksign,x,y)
      implicit none
      integer, intent(in) :: ksign
      complex, dimension(:,:,:), intent(in) :: x
      complex, dimension(:,:,:), intent(out) :: y
!hpf$ distribute(*,*,block) :: y
!hpf$ distribute(*,*,block) :: x
      end subroutine mfft_local2
      end interface

!-----FFTs along all dimensions except the last:

      call mfft_local2 (ksign,x,x2)

!-----The 3D and higher, the previous call transposes the subarray 
!-----of all but the last index; Now move the last index from the 
!-----right-most position to the left-most position in preparation 
!-----for final, in-place transformation: In 2D this is simply a 
!-----transpose of the entire array.

!hpf$ independent, new(j)
      do i=1,n1
!hpf$ independent, new(k)
      do j=1,n2
!hpf$ independent
      do k=1,n3
      x3(k,j,i)=x2(j,i,k)
      enddo
      enddo
      enddo

!-----Final transform along the left-most direction

      call mfft_local1 (ksign,x3)

      if(scale.ne.1.0)x3=x3*scale

!-----if icpy.eq.1, store the final result back in x
!-----if icpy.ne.1, store with zeros (ensure users know it's filled 
!-----with junk)

      if(icpy.eq.1)then
!hpf$ independent, new(j)
      do i=1,n1
!hpf$ independent, new(k)
      do j=1,n2
!hpf$ independent
      do k=1,n3
      x(i,j,k)=x3(k,j,i)
      enddo
      enddo
      enddo

      else
      x=(0.,0.)

      endif
      return
      end subroutine fft3d

!-----HPF_Local subroutine to perform multiple FFTs on the inner 
!-----dimension

      extrinsic (hpf_local) subroutine mfft_local1 (ksign,x)
          
      implicit none
      
!-----scalars
    
      integer, intent(in) :: ksign

!-----arrays

      complex, dimension(:,:,:), intent(inout) :: x

!-----distribute/aligns

!hpf$ distribute(*,*,block) :: x

      integer :: j,k

      real :: uscal

      complex, dimension(size(x,1)) :: table,tempi,tempo
      complex, dimension(2*size(x,1)) :: work

!-----interface ccfftnr

      interface
      extrinsic (hpf_local) subroutine ccfftnr(cdata,ksign)
      implicit none
      complex, intent(inout) :: cdata(:)
!hpf$ distribute cdata(*)
      integer, intent(in) :: ksign
      end subroutine ccfftnr
      end interface

!-----Initialize:

      uscal=1.0

!-----Perform multiple FFTs without scaling:

!hpf$ independent, new(j)
      do k = 1,size(x,3)
!hpf$ independent
      do j = 1,size(x,2)
             call ccfftnr(x(:,j,k),ksign)
      end do
      end do
      return
      end subroutine mfft_local1

!-----HPF_Local subroutine to perform multiple FFTs on all but 
!-----last dimension.

      extrinsic (hpf_local) subroutine mfft_local2(ksign,x,y)

      implicit none

!-----scalars

      integer, intent(in) :: ksign

!-----arrays
          
      complex, dimension(:,:,:), intent(in) :: x
      complex, dimension(:,:,:), intent(out) :: y

!-----distribute/aligns

!hpf$ distribute(*,*,block) :: x
!hpf$ distribute(*,*,block) :: y

      integer :: i, j, k
      
      real :: uscal

      complex, dimension(max(size(y,1),size(y,2))) :: table
      complex, dimension(size(y,1)) :: temp1
      complex, dimension(size(y,2)) :: temp2
      complex, dimension(2*max(size(y,1),size(y,2))) :: work

!hpf$ distribute (*) :: temp1
!hpf$ distribute (*) :: temp2
!hpf$ distribute (*) :: work

!-----interface

      interface
      extrinsic (hpf_local) subroutine ccfftnr(cdata,ksign)
      implicit none
      complex, intent(inout) :: cdata(:)
!hpf$ distribute (*) :: cdata
      integer, intent(in) :: ksign
      end subroutine ccfftnr
      end interface

!-----Initialize:

      uscal=1.0

!-----Perform multiple FFTs without scaling.

!hpf$ independent, new(i)
      do k = 1, size (x,3)
!hpf$ independent
      do i = 1, size (x,2)
      temp2 = x(:,i,k)
      call ccfftnr(temp2,ksign)
      y(i,:,k) = temp2
      end do
      end do

!-----Perform multiple FFTs without scaling.
        
!hpf$ independent, new(j)  
      do k=1,size(y,3)
!hpf$ independent
      do j=1,size(y,2)
      temp1 = y(:,j,k)
      call ccfftnr(temp1,ksign)
      y(:,j,k) = temp1
      end do
      end do
      return

      end subroutine mfft_local2

!-----Numerical Recipes FFT (sigh --)

      extrinsic (hpf_local) subroutine ccfftnr(cdata,ksign)
      
!-----Note Rynian wierdness
  
      implicit real(a-h,o-z)

      complex cdata(:)
!hpf$ distribute (*) :: cdata
      integer :: nn
      dimension data(2*size(cdata))

      nn = size(cdata)

!hpf$ independent
      do i=1,nn
      data(2*i-1)=real(cdata(i))
      data(2*i) =aimag(cdata(i))
      enddo

!-----bit reversal:

      n=2*nn
      j=1

      do i=1,n,2

      if(j.gt.i)then
      tempr=data(j)
      tempi=data(j+1)
      data(j)=data(i)
      data(j+1)=data(i+1)
      data(i)=tempr
      data(i+1)=tempi
      endif

      m=n/2
  1   if((m.ge.2).and.(j.gt.m))then
      j=j-m
      m=m/2
      goto 1
      endif
      j=j+m
      enddo

!-----Danielson-Lanczos:

      twopi=4.0*asin(1.0)
      mmax=2
  2   if(n.gt.mmax)then
      istep=2*mmax
      theta=twopi/(ksign*mmax)
      wpr=-2.*sin(0.5*theta)**2
      wpi=sin(theta)
      wr=1.0
      wi=0.0

      do m=1,mmax,2
      do i=m,n,istep
      j=i+mmax
      tempr=wr*data(j)-wi*data(j+1)
      tempi=wr*data(j+1)+wi*data(j)
      data(j)=data(i)-tempr
      data(j+1)=data(i+1)-tempi
      data(i)=data(i)+tempr
      data(i+1)=data(i+1)+tempi
      enddo
      wtemp=wr
      wr=wr*wpr-wi*wpi+wr
      wi=wi*wpr+wtemp*wpi+wi
      enddo

      mmax=istep
      goto 2
      endif

      do i=1,nn
      cdata(i)=data(2*i-1)+(0.,1.)*data(2*i)
      enddo
      
      return
      end


