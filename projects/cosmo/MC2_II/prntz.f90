      subroutine prntz(coor,np,rL,ng,skip,frm,pbchoice)
          
!-----This subroutine prints out coordinates and nothing else
!-----during the Zeldovich step

      implicit none

!-----scalars

      integer np,ii,mm,ng,skip,frm,pbchoice

      real rL,xscal,vscal

!-----arrays

      real coor(6,np),tmp(np)

!-----tmp:coor
!-----distribute/aligns

!hpf$ distribute coor(*,block)
!hpf$ align (:) with coor(*,:) :: tmp 

      write(6,*)'entering prntz'

!-----Print out coordinates and velocities of the particles (box units)

      if(pbchoice.eq.1)then
      if(frm.eq.1)then

      open(66,status='unknown')

      do  ii=1,np,skip
      write(66,103)coor(1,ii),coor(2,ii),coor(3,ii),coor(4,ii),
     #             coor(5,ii),coor(6,ii)
      enddo

      close(66)

      endif

      if(frm.eq.2)then

      open(66,form='unformatted',access='sequential',
     #     recordtype='stream')
      open(67,form='unformatted',access='sequential',
     #     recordtype='stream')
      open(68,form='unformatted',access='sequential',
     #     recordtype='stream')
      open(69,form='unformatted',access='sequential',
     #     recordtype='stream')
      open(70,form='unformatted',access='sequential',
     #     recordtype='stream')
      open(71,form='unformatted',access='sequential',
     #     recordtype='stream')

      write(66)coor(1,:)
      write(67)coor(2,:)
      write(68)coor(3,:)
      write(69)coor(4,:)
      write(70)coor(5,:)
      write(71)coor(6,:)

      close(66)
      close(67)
      close(68)
      close(69)
      close(70)
      close(71)

      endif
      endif

!-----Print out coordinates and velocities of the particles (physical units)
!-----Position is in h^(-1) Mpc, comoving peculiar velocity in km/s.

      if(pbchoice.eq.2)then

      xscal=rL/(1.0*ng)
      vscal=100.0*xscal

      if(frm.eq.1)then
      open(66,status='unknown')

      do  ii=1,np,skip
      write(66,103)xscal*coor(1,ii),vscal*coor(2,ii),xscal*coor(3,ii),
     #             vscal*coor(4,ii),xscal*coor(5,ii),vscal*coor(6,ii)
      enddo

      close(66)

      endif

      if(frm.eq.2)then

      tmp(:)=xscal*coor(1,:)

      open(66,form='unformatted',access='sequential',
     #     recordtype='stream')
      write(66)tmp(:)
      close(66)

      tmp(:)=vscal*coor(2,:)

      open(67,form='unformatted',access='sequential',
     #     recordtype='stream')
      write(67)tmp(:)
      close(67)

      tmp(:)=xscal*coor(3,:)

      open(68,form='unformatted',access='sequential',
     #     recordtype='stream')
      write(68)tmp(:)
      close(68)

      tmp(:)=vscal*coor(4,:)

      open(69,form='unformatted',access='sequential',
     #     recordtype='stream')
      write(69)tmp(:)
      close(69)

      tmp(:)=xscal*coor(5,:)

      open(70,form='unformatted',access='sequential',
     #     recordtype='stream')
      write(70)tmp(:)
      close(70)

      tmp(:)=vscal*coor(6,:)

      open(71,form='unformatted',access='sequential',
     #     recordtype='stream')
      write(71)tmp(:)
      close(71)

      endif
      endif

103   format(7(e12.5,1x))

      return
      end subroutine
