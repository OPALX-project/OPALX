      subroutine zeldo(dplus,dzero,ddot,gradx,grady,gradz,ng,
     #                 np,rL,zin,coor,skip,frm,initp,pbchoice,
     #                 gpscal)

!-----Implements Zeldovich displacements of particles.
    
      implicit none

!-----scalars

      integer ng,np,mm,ii,jj,nopt,skip,frm,initp,pbchoice

      real ain,fn,tpi,tpiL,zin
      real ppin,dplus,ddot,dzero
      real rL,rscale,gpscal

!-----arrays

      integer indx(np),jndx(np),kndx(np)
      integer indxp1(np),jndxp1(np),kndxp1(np)

      real pkf(ng/2)
      real ab(np),de(np),gh(np)
      real vec(np),fx(np),fy(np),fz(np)
      real coor(6,np)
      real gradx(ng,ng,ng),grady(ng,ng,ng),gradz(ng,ng,ng)
      real fxg(ng,ng,ng),fyg(ng,ng,ng),fzg(ng,ng,ng),rho(ng,ng,ng)

!-----indx,jndx,kndx,indxp1,jndxp1,kndxp1,ab,de,gh,vec,fx,fy,fz:coor
!-----gradx,grady,gradz,fxg,fyg,fzg:rho
!-----distribute/aligns

!hpf$ distribute pkf(block)
!hpf$ distribute coor(*,block)
!hpf$ align (:) with coor(*,:) :: indx,jndx,kndx,indxp1,jndxp1,kndxp1
!hpf$ align (:) with coor(*,:) :: ab,de,gh,vec,fx,fy,fz
!hpf$ distribute rho(*,*,block)
!hpf$ align (*,*,:) with rho(*,*,:) :: gradx,grady,gradz
!hpf$ align (*,*,:) with rho(*,*,:) :: fxg,fyg,fzg

!-----Define ain

      ain=1.0/(1.0+zin)

!-----Put particles on grid

      coor=0.0

      do ii=1,ng
!hpf$ independent, new(mm)
      do jj=1,ng
!hpf$ independent
      do mm=1,ng

      coor(1,ng*ng*(mm-1)+ng*(jj-1)+ii)=1.0*(ii-1)
      coor(3,ng*ng*(mm-1)+ng*(jj-1)+ii)=1.0*(jj-1)
      coor(5,ng*ng*(mm-1)+ng*(jj-1)+ii)=1.0*(mm-1)

      enddo
      enddo
      enddo

!-----Re-index force on particles on grid point (i.e., convert gradx, 
!-----grady, and gradz to vectors in particle index space)

      do ii=1,ng
!hpf$ independent, new(mm)
      do jj=1,ng
!hpf$ independent
      do mm=1,ng

      fx(ng*ng*(mm-1)+ng*(jj-1)+ii)=gradx(ii,jj,mm)
      fy(ng*ng*(mm-1)+ng*(jj-1)+ii)=grady(ii,jj,mm)
      fz(ng*ng*(mm-1)+ng*(jj-1)+ii)=gradz(ii,jj,mm)

      enddo
      enddo
      enddo

!-----Zeldovich displacement = (dplus/dzero)*fi. 

      forall(ii=1:np)coor(1,ii)=coor(1,ii)-(dplus/dzero)*fx(ii)
      forall(ii=1:np)coor(3,ii)=coor(3,ii)-(dplus/dzero)*fy(ii)
      forall(ii=1:np)coor(5,ii)=coor(5,ii)-(dplus/dzero)*fz(ii)

!-----Zeldovich (comoving, peculiar) velocities. Time measured in units
!-----of H_0^(-1).

      forall(ii=1:np)coor(2,ii)=0.0-(ddot/dzero)*fx(ii) 
      forall(ii=1:np)coor(4,ii)=0.0-(ddot/dzero)*fy(ii) 
      forall(ii=1:np)coor(6,ii)=0.0-(ddot/dzero)*fz(ii)        

!-----Print particle positions and velocities before periodic BCs
!-----initp=0

      if(initp.eq.0)then
      call prntz(coor,np,rL,ng,skip,frm,pbchoice)
      endif

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

!-----Print particle positions and velocities after periodic BCs
!-----initp=1
      
      if(initp.eq.1)then
      call prntz(coor,np,rL,ng,skip,frm,pbchoice)
      endif

!-----Compute CIC density field from particle positions
!-----Normalized to particle number density!
      
      call rhoord(coor,rho,np,ng,ab,de,gh,
     #      indx,jndx,kndx,indxp1,jndxp1,kndxp1,1.0)
      
!-----Compute 1-D power spectrum (physical units), remember to normalize
!-----the density 

      rho=(dzero/dplus)*(rho-1.0)

      do ii=1,ng
      write(32,*) ii,rho(ng/2,ng/2,ii)
      enddo

      call pspec1(rho,ng,pkf)

      tpi=8.0*atan(1.0)
      tpiL=tpi/rL
      rscale=rL/(1.0*ng)

      do mm=1,ng/2
      write(33,*) mm*tpiL,pkf(mm)*rscale**3
      enddo

      call psp1cic(rho,ng,pkf)

      do mm=1,ng/2
      write(34,*) mm*tpiL,pkf(mm)*rscale**3
      enddo

      return
      end subroutine










