      subroutine prnt(istp,coor,np,rL,ng,skip,frm,pp,alpha,
     #                pbchoice,omegatot,hubble)
           
!-----This subroutine prints out coordinates and nothing else!

      implicit none

!-----scalars

      integer istp,np,ifile,k,mm,ng,skip,frm,pbchoice

      real rL,pp,alpha,omegatot,hubble,aa,xscal,vscal

!-----arrays

      real coor(6,np),tmp(np)

!-----tmp:coor
!-----distribute/aligns

!hpf$ distribute coor(*,block)
!hpf$ align (:) with coor(*,:) :: tmp

      aa=pp**(1.0/alpha)

      ifile=istp+111

!-----Print out coordinates and velocities of the particles (box units)

      if(pbchoice.eq.1)then

!-----Print out coordinates and velocities of the particles in ASCI

      if(frm.eq.1)then

      open(ifile,status='unknown')

      do  k=1,np,skip
      write(ifile,104)coor(1,k),coor(2,k)/pp**(2.0/alpha),
     #                coor(3,k),coor(4,k)/pp**(2.0/alpha),
     #                coor(5,k),coor(6,k)/pp**(2.0/alpha)
      enddo

      close(ifile)

      endif

!-----Print out in binary
     
      if(frm.eq.2)then

      open(ifile,form='unformatted',access='sequential',
     #     recordtype='stream') 
      open(ifile+1,form='unformatted',access='sequential',
     #     recordtype='stream') 
      open(ifile+2,form='unformatted',access='sequential',
     #     recordtype='stream') 
      open(ifile+3,form='unformatted',access='sequential',
     #     recordtype='stream') 
      open(ifile+4,form='unformatted',access='sequential',
     #     recordtype='stream') 
      open(ifile+5,form='unformatted',access='sequential',
     #     recordtype='stream') 

      write(ifile)coor(1,:)
      write(ifile+1)coor(2,:)
      write(ifile+2)coor(3,:)
      write(ifile+3)coor(4,:)
      write(ifile+4)coor(5,:)
      write(ifile+5)coor(6,:)

      close(ifile)
      close(ifile+1)
      close(ifile+2)
      close(ifile+3)
      close(ifile+4)
      close(ifile+5)

      endif

      endif

!-----Print out coordinates and velocities of the particles (physical units)
!-----Position is in h^(-1) Mpc, comoving peculiar velocity in km/s.

      if(pbchoice.eq.2)then

      xscal=rL/(1.0*ng)
      vscal=xscal*100.0/pp**2

!-----Print out coordinates and velocities of the particles in ASCI

      if(frm.eq.1)then

      open(ifile,status='unknown')
     
      do  k=1,np,skip
      write(ifile,104)xscal*coor(1,k),vscal*coor(2,k),xscal*coor(3,k),
     #                vscal*coor(4,k),xscal*coor(5,k),vscal*coor(6,k)
      enddo

      
!      do  k=1,np,skip
!      write(ifile+2,103)k,xscal*coor(1,k),xscal*coor(3,k),
!     #                  xscal*coor(5,k),vscal*coor(2,k),
!     #                  vscal*coor(4,k),vscal*coor(6,k),
!     #                  2.3217082E+11
!      write(ifile,103)k,xscal*coor(1,k),xscal*coor(3,k),
!     #                  xscal*coor(5,k),2.3217082E+11
!      enddo
      

      close(ifile)

      endif

!-----Print out coordinates and velocities of the particles in binary

      if(frm.eq.2)then

      tmp(:)=xscal*coor(1,:)

      open(ifile,form='unformatted',access='sequential',
     #     recordtype='stream')
      write(ifile)tmp(:)
      close(ifile)

      tmp(:)=vscal*coor(2,:)

      open(ifile+1,form='unformatted',access='sequential',
     #     recordtype='stream')
      write(ifile+1)tmp(:)
      close(ifile+1)

      tmp(:)=xscal*coor(3,:)

      open(ifile+2,form='unformatted',access='sequential',
     #     recordtype='stream')
      write(ifile+2)tmp(:)
      close(ifile+2)

      tmp(:)=vscal*coor(4,:)

      open(ifile+3,form='unformatted',access='sequential',
     #     recordtype='stream')
      write(ifile+3)tmp(:)
      close(ifile+3)

      tmp(:)=xscal*coor(5,:)

      open(ifile+4,form='unformatted',access='sequential',
     #     recordtype='stream')
      write(ifile+4)tmp(:)
      close(ifile+4)

      tmp(:)=vscal*coor(6,:)

      open(ifile+5,form='unformatted',access='sequential',
     #     recordtype='stream')
      write(ifile+5)tmp(:)
      close(ifile+5)

      close(ifile)
      close(ifile+1)
      close(ifile+2)
      close(ifile+3)
      close(ifile+4)
      close(ifile+5)

      endif

      endif     

  103 format(1X,1(i5),7(e12.4))
  104 format(1X,6(e12.4))

      return
      end subroutine
