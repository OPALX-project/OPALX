      subroutine dens2d(ifile,coor,np,ng,hubble,pmass,ng2d,nn,rL,frm)

      implicit none

      integer ii,jj,mm,np,ng,ng2d,ifile,frm

      real hubble,pmass,nn,rL

      real coor3d(3,np),coor2d(2,np),rho2d(ng2d,ng2d)
      real coordsc(3,np),coor(6,np)

!hpf$ distribute coor3d(*,block)
!hpf$ distribute coor2d(*,block)
!hpf$ distribute coor(*,block)
!hpf$ align (*,:) with coor3d(*,:) :: coordsc
!hpf$ distribute rho2d(*,block)

      coor3d(1,:)=coor(1,:)
      coor3d(2,:)=coor(3,:)
      coor3d(3,:)=coor(5,:)

!-----Generate scaled coordinate array

      coordsc=nn*coor3d

!-----Projection into xy-plane

      coor2d(1,:)=coordsc(1,:)
      coor2d(2,:)=coordsc(2,:)

      call rhoord2d(coor2d,rho2d,np,ng2d,nn)

      if(frm.eq.3)then      

      open(ifile+100,form='formatted')

      do ii=1,ng2d
      do jj=1,ng2d

      write(ifile+100,103)ii,jj,rho2d(ii,jj)*pmass*(ng2d/rL)**2

      enddo

      write(ifile+100,*)

      enddo

      close(ifile+100)

      do ii=1,np
      write(1,*)coor2d(1,ii),coor2d(2,ii)
      enddo

      endif

      if(frm.eq.4)then

      open(ifile+100,form='unformatted',access='sequential',
     #     recordtype='stream')

      write(ifile+100)rho2d(:,:)*pmass*(ng2d/rL)**2

      close(ifile+100)

      endif

!-----Projection into xz-plane

      coor2d(1,:)=coordsc(1,:)
      coor2d(2,:)=coordsc(3,:)

      call rhoord2d(coor2d,rho2d,np,ng2d,nn)

      if(frm.eq.3)then

      open(ifile+101,form='formatted')
     
      do ii=1,ng2d
      do jj=1,ng2d

      write(ifile+101,103)ii,jj,rho2d(ii,jj)*pmass*(ng2d/rL)**2

      enddo

      write(ifile+101,*)

      enddo

      close(ifile+101)

      do ii=1,np
      write(2,*)coor2d(1,ii),coor2d(2,ii)
      enddo

      endif

      if(frm.eq.4)then

      open(ifile+101,form='unformatted',access='sequential',
     #     recordtype='stream')

      write(ifile+101)rho2d(:,:)*pmass*(ng2d/rL)**2

      close(ifile+101)

      endif

!-----Projection into yz-plane

      coor2d(1,:)=coordsc(2,:)
      coor2d(2,:)=coordsc(3,:)

      call rhoord2d(coor2d,rho2d,np,ng2d,nn)

      if(frm.eq.3)then

      open(ifile+102,form='formatted')
     
      do ii=1,ng2d
      do jj=1,ng2d

      write(ifile+102,103)ii,jj,rho2d(ii,jj)*pmass*(ng2d/rL)**2

      enddo

      write(ifile+102,*)

      enddo

      close(ifile+102)

      do ii=1,np
      write(3,*)coor2d(1,ii),coor2d(2,ii)
      enddo

      endif

      if(frm.eq.4)then

      open(ifile+102,form='unformatted',access='sequential',
     #     recordtype='stream')

      write(ifile+102)rho2d(:,:)*pmass*(ng2d/rL)**2

      close(ifile+102)

      endif

103   format(1x,2(i5),e12.5)

      return
      end subroutine     
