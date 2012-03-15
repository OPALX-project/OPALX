      subroutine cont(ng,np,rL,coor)

!-----Reads in initial conditions from a previous code run

      implicit none

!-----scalars

      integer ng,np

      real rL,xscal,vscal

!-----arrays

      real coor(6,np)

!hpf$ distribute coor(*,block)

!-----Read in Initial Conditions
!-----Only in physical units and binary format!

      write(6,*)'entering read-in'

      open(1,file='in.191',form='unformatted',access='sequential',
     #     recordtype='stream')
      read(1)coor(1,:)
      close(1)

      open(2,file='in.192',form='unformatted',access='sequential',
     #     recordtype='stream')
      read(2)coor(2,:)
      close(2)

      open(3,file='in.193',form='unformatted',access='sequential',
     #     recordtype='stream')
      read(3)coor(3,:)
      close(3)

      open(4,file='in.194',form='unformatted',access='sequential',
     #     recordtype='stream')
      read(4)coor(4,:)
      close(4)

      open(5,file='in.195',form='unformatted',access='sequential',
     #     recordtype='stream')
      read(5)coor(5,:)
      close(5)

      open(6,file='in.196',form='unformatted',access='sequential',
     #     recordtype='stream')
      read(6)coor(6,:)
      close(6)

!-----Converting to box-units for further use in the code:

      xscal=rL/(1.0*ng)
      vscal=100.0*xscal

      coor(1,:)=coor(1,:)/xscal
      coor(2,:)=coor(2,:)/vscal
      coor(3,:)=coor(3,:)/xscal
      coor(4,:)=coor(4,:)/vscal
      coor(5,:)=coor(5,:)/xscal
      coor(6,:)=coor(6,:)/vscal

      write(6,*)'read-in finished'

103   format(7(e12.5,1x))

      return
      end subroutine












