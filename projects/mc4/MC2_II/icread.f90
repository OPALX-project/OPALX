      subroutine icread(ng,np,rL,zin,coor,frm,pbchoice,icbc,gpscal)

!-----Reads in initial conditions provided from outside, icbc controls
!-----whether before or after enforcing periodic BCs

      implicit none

      integer ng,np,mm,ii,jj,nopt,skip,frm,pbchoice,icbc

      integer indx(np),jndx(np),kndx(np)
      integer indxp1(np),jndxp1(np),kndxp1(np)

      real fn,tpi,tpiL,zin
      real rL,rscale,rLb,gpscal
      real xscal,vscal,alpha,omegatot,pp

      real coor(6,np),vec(np),xpot(ng,ng,ng)
      real pkf(ng/2)
      real rho(ng,ng,ng)
      real ab(np),de(np),gh(np)

!hpf$ distribute pkf(block)
!hpf$ distribute coor(*,block)
!hpf$ align (:) with coor(*,:) :: vec
!hpf$ align (:) with coor(*,:) :: indx,jndx,kndx,ab,de,gh
!hpf$ align (:) with coor(*,:) :: indxp1,jndxp1,kndxp1
!hpf$ distribute rho(*,*,block)
!hpf$ align (*,*,:) with rho(*,*,:) :: xpot

      tpi=8.0*atan(1.0)  ! 2 pi
      rLb=1.0*ng

!-----Read in Initial Conditions

!-----Initial conditions in box units

      if(pbchoice.eq.1)then

!-----Initial conditions given in binary

      if(frm.eq.2)then

!-----Initial conditions in box units
 
      open(1,file='fort.66',form='unformatted',access='sequential',
     #     recordtype='stream')
      open(2,file='fort.67',form='unformatted',access='sequential',
     #     recordtype='stream')
      open(3,file='fort.68',form='unformatted',access='sequential',
     #     recordtype='stream')
      open(4,file='fort.69',form='unformatted',access='sequential',
     #     recordtype='stream')
      open(5,file='fort.70',form='unformatted',access='sequential',
     #     recordtype='stream')
      open(6,file='fort.71',form='unformatted',access='sequential',
     #     recordtype='stream')

      read(1)coor(1,:)
      read(2)coor(2,:)
      read(3)coor(3,:)
      read(4)coor(4,:)
      read(5)coor(5,:)
      read(6)coor(6,:)

      close(1)
      close(2)
      close(3)
      close(4)
      close(5)
      close(6)

      endif     

!-----Initial conditions given in ASCI

      if(frm.eq.1)then

!-----Initial conditions in box units

      open(1,file="fort.66",status="old")

      do ii=1,np
      read(1,103)coor(1,ii),coor(2,ii),coor(3,ii),
     #           coor(4,ii),coor(5,ii),coor(6,ii)
      enddo

      close(1)

      endif     

      endif

!-----Initial conditions in physical units

      if(pbchoice.eq.2)then

!-----in binary 

      if(frm.eq.2)then

      open(1,file='fort.66',form='unformatted',access='sequential',
     #     recordtype='stream')
      open(2,file='fort.67',form='unformatted',access='sequential',
     #     recordtype='stream')
      open(3,file='fort.68',form='unformatted',access='sequential',
     #     recordtype='stream')
      open(4,file='fort.69',form='unformatted',access='sequential',
     #     recordtype='stream')
      open(5,file='fort.70',form='unformatted',access='sequential',
     #     recordtype='stream')
      open(6,file='fort.71',form='unformatted',access='sequential',
     #     recordtype='stream')

      read(1)coor(1,:)
      read(2)coor(2,:)
      read(3)coor(3,:)
      read(4)coor(4,:)
      read(5)coor(5,:)
      read(6)coor(6,:)

      close(1)
      close(2)
      close(3)
      close(4)
      close(5)
      close(6)

      endif     

!-----in ASCII

      if(frm.eq.1)then

      open(1,file="fort.66",status="old")

      do ii=1,np
      read(1,103)coor(1,ii),coor(2,ii),coor(3,ii),
     #           coor(4,ii),coor(5,ii),coor(6,ii)
      enddo

      close(1)

      endif

!-----Converting to box-units for further use in the code:
     
      xscal=rL/(1.0*ng)
      vscal=100.0*xscal

      coor(1,:)=coor(1,:)/xscal
      coor(2,:)=coor(2,:)/vscal
      coor(3,:)=coor(3,:)/xscal
      coor(4,:)=coor(4,:)/vscal
      coor(5,:)=coor(5,:)/xscal
      coor(6,:)=coor(6,:)/vscal

      endif

!-----Print first 100 particle positions and velocities (for test purposes)
 
      do ii=1,100
      write(11,103)coor(1,ii),coor(2,ii),coor(3,ii),
     #             coor(4,ii),coor(5,ii),coor(6,ii)
      enddo

!-----Enforce periodic boundary conditions if needed

      if(icbc.eq.0)then

      fn=1.0*ng
      vec=coor(1,:)
      where(vec.lt.0.)vec=vec+fn
      where(vec.ge.fn)vec=vec-fn
      coor(1,:)=vec

      fn=1.0*ng
      vec=coor(3,:)
      where(vec.lt.0.)vec=vec+fn
      where(vec.ge.fn)vec=vec-fn
      coor(3,:)=vec

      fn=1.0*ng
      vec=coor(5,:)
      where(vec.lt.0.)vec=vec+fn
      where(vec.ge.fn)vec=vec-fn
      coor(5,:)=vec

      endif

!-----Compute CIC density field from particle positions
!-----Normalized to particle number density.

      call rhoord(coor,rho,np,ng,ab,de,gh,
     #      indx,jndx,kndx,indxp1,jndxp1,kndxp1,1.0)

!      do ii=1,ng
!      write(32,*) ii,rho(ng/2,ng/2,ii)
!      enddo

!-----Compute 1-D power spectrum (physical units)

      call pspec1(rho,ng,pkf)

      tpiL=tpi/rL
      rscale=rL/(1.0*ng)

      do mm=1,ng/2
      write(33,*) mm*tpiL,pkf(mm)*rscale**3
      enddo

103   format(7(e12.5,1x))

      return
      end subroutine
