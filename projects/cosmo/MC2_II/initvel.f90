      subroutine initvel(coor,np,ng,egtldetr,dplus,ddot)

!-----Generation of `self-consistent' Zeldovich velocities from the 
!-----numerically obtained gravitational potential after the initial 
!-----Zeldovich displacement. Not obvious why one needs this; some 
!-----testing is required (velocity dispersion).

      implicit none

!-----scalars

      integer np,ng,ii,jj,mm,ksign,icpy

      real ct,dplus,ddot,scal,fn
      real meanvx,meanvy,meanvz,sigmavx,sigmavy,sigmavz
      real skewvx,skewvy,skewvz,kurtvx,kurtvy,kurtvz

!-----arrays

      integer indx(np),jndx(np),kndx(np),indxp1(np)
      integer jndxp1(np),kndxp1(np)

      real coor(6,np),fx(np),fy(np),fz(np)
      real ab(np),de(np),gh(np),vec(np)
      real fxg(ng,ng,ng),fyg(ng,ng,ng),fzg(ng,ng,ng),rho(ng,ng,ng)

      complex egtldetr(ng,ng,ng),erhotr(ng,ng,ng),erho(ng,ng,ng)

!-----indx,jndx,kndx,indxp1,jndxp1,kndxp1,fx,fy,fz,ab,de,gh,vec:coor
!-----fxg,fyg,fzg,erho:rho
!-----erhotr:egtldetr
!-----distribute/aligns

!hpf$ distribute coor(*,block)
!hpf$ align (:) with coor(*,:) :: indx,jndx,kndx,indxp1,jndxp1,kndxp1
!hpf$ align (:) with coor(*,:) :: fx,fy,fz,ab,de,gh,vec
!hpf$ distribute rho(*,*,block)
!hpf$ align (*,*,:) with rho(*,*,:) :: fxg,fyg,fzg,erho
!hpf$ distribute egtldetr(*,*,block)
!hpf$ align (*,*,:) with egtldetr(*,*,:) :: erhotr

!-----Step I: Obtain new grid density field from the particle
!-----distribution. [Note: This is not obviously the initial density 
!-----field which we used to obtain the particle positions in the
!-----first place.] Deposit particles on the grid using CIC:

      call rhoord(coor,rho,np,ng,ab,de,gh,
     #            indx,jndx,kndx,indxp1,jndxp1,kndxp1)

!-----Step II: Now solve the Poisson equation. First, Fourier 
!-----transform the mass density:

      erho=rho

      ksign=1
      scal=1.0
      icpy=0

      call fft3d(ng,ng,ng,ksign,scal,icpy,erho,erhotr)

!-----Then convolve the mass density with the Green's function:

      erhotr=erhotr*egtldetr

      ksign=-1

      call fft3d (ng, ng, ng, ksign, scal, icpy, erhotr, erho)

!-----Compensate for inverse FFT scaling

      erho=erho/(1.0*ng*ng*ng)

!-----Final solution (rho is really phi)

      rho=real(erho,8)

!-----A useful diagnostic to compare against the initial on-grid
!-----bare potential field

      do mm=1,ng
      write(92,*) mm,rho(ng/2,ng/2,mm)
      enddo

!-----Compute the resulting gradient field (on-grid)

!------DISABLED TEMPORARILY!
!      fxg=0.5*(cshift(rho,1,1)-cshift(rho,-1,1))
!      fyg=0.5*(cshift(rho,1,2)-cshift(rho,-1,2))
!      fzg=0.5*(cshift(rho,1,3)-cshift(rho,-1,3))

!-----Start particles once again on-grid 

      forall(ii=1:ng,jj=1:ng,mm=1:ng) 
     #        coor(1,ng*ng*(mm-1)+ng*(jj-1)+ii)=ii-1
      forall(ii=1:ng,jj=1:ng,mm=1:ng) 
     #        coor(3,ng*ng*(mm-1)+ng*(jj-1)+ii)=jj-1
      forall(ii=1:ng,jj=1:ng,mm=1:ng) 
     #        coor(5,ng*ng*(mm-1)+ng*(jj-1)+ii)=mm-1

!-----Re-index force on particles on grid point (i.e., convert fxg, 
!-----fyg, and fzg to vectors in particle index space)

      forall(ii=1:ng,jj=1:ng,mm=1:ng)fx(ng*ng*(mm-1)+ng*(jj-1)+ii)
     #       =fxg(ii,jj,mm)
      forall(ii=1:ng,jj=1:ng,mm=1:ng)fy(ng*ng*(mm-1)+ng*(jj-1)+ii)
     #       =fyg(ii,jj,mm)
      forall(ii=1:ng,jj=1:ng,mm=1:ng)fz(ng*ng*(mm-1)+ng*(jj-1)+ii)
     #       =fzg(ii,jj,mm)

!-----Zeldovich move using new gradient field

      forall(ii=1:np)coor(1,ii)=coor(1,ii)-dplus*fx(ii)
      forall(ii=1:np)coor(3,ii)=coor(3,ii)-dplus*fy(ii)
      forall(ii=1:np)coor(5,ii)=coor(5,ii)-dplus*fz(ii)

!-----Enforce periodic boundary conditions

      fn=1.0*ng
      vec=coor(1,:)
      where(vec.lt.0.)vec=vec+fn
      where(vec.ge.fn)vec=vec-fn
      coor(1,:)=vec

      vec=coor(3,:)
      where(vec.lt.0.)vec=vec+fn
      where(vec.ge.fn)vec=vec-fn
      coor(3,:)=vec

      vec=coor(5,:)
      where(vec.lt.0.)vec=vec+fn
      where(vec.ge.fn)vec=vec-fn
      coor(5,:)=vec
  
!-----Compute new velocities

      coor(2,:)=-fx(:)*ddot
      coor(4,:)=-fy(:)*ddot
      coor(6,:)=-fz(:)*ddot

!-----Compute the first four moments of the velocity distribution

      meanvx=sum(coor(2,:))
      meanvy=sum(coor(4,:))
      meanvz=sum(coor(6,:))

      sigmavx=sqrt(sum((coor(2,:)-meanvx)**2))
      sigmavy=sqrt(sum((coor(4,:)-meanvy)**2))
      sigmavz=sqrt(sum((coor(6,:)-meanvz)**2))
      skewvx=sum((coor(2,:)-meanvx)**3)
      skewvy=sum((coor(4,:)-meanvy)**3)
      skewvz=sum((coor(6,:)-meanvz)**3)
      kurtvx=np*sum((coor(2,:)-meanvx)**4)/sigmavx**4-3.0
      kurtvy=np*sum((coor(4,:)-meanvy)**4)/sigmavy**4-3.0
      kurtvz=np*sum((coor(6,:)-meanvz)**4)/sigmavz**4-3.0

      print*,'meanvx,meanvy,meanvy   =',meanvx/(1.0*np),meanvy/(1.0*np),
     #                                meanvz/(1.0*np)
      print*,'sigmavx,sigmavy,sigmavz=',sigmavx/(1.0*np),
     #                                sigmavy/(1.0*np),sigmavz/(1.0*np)
      print*,'skewvx,skewvy,skewvz   =',skewvx/(1.0*np),skewvy/(1.0*np),
     #                                skewvz/(1.0*np)
      print*,'kurtvx,kurtvy,kurtvz   =',kurtvx,kurtvy,kurtvz

!-----Compute CIC density field from particle positions
!-----Note: Not normalized to physical density!

      call rhoord(coor,rho,np,ng,ab,de,gh,
     #      indx,jndx,kndx,indxp1,jndxp1,kndxp1)


      do mm=1,np
      write(17,*) mm,coor(1,mm),coor(3,mm),coor(5,mm)
      enddo

      do mm=1,np
      write(18,*) mm,coor(2,mm),coor(4,mm),coor(6,mm)
      enddo

      do ii=1,ng
      write(70,*) ii,fxg(ng/2,ng/2,ii)
      write(71,*) ii,fyg(ng/2,ng/2,ii)
      write(72,*) ii,fzg(ng/2,ng/2,ii)
      enddo

      return
      end subroutine







