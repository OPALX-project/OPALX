      program universe

!-----Even more cleaned up version (October 22, 2002)
!-----Cleaned up version (begins June 30, 2002)

!-----Historical Note: 3D large scale structure code 
!-----(8/28/95) Habib, Ryne, Sathya

!-----Machines: CM-5, SGI Origin 2000, T3E, SP Power 3, Compaq Alpha
!-----Beowulf (Intel Pentium/Xeon)

!-----This routine reads in the input information from indat and passes
!-----it on to main. Basically it does nothing (except for the particle
!-----hardwire) - its job is to function as a dynamical memory allocation
!-----tool.

      implicit none

      integer ng,np,nsteps,norder,iprint,iseed,jinit,norm,trans
      integer nxint,nxdint,initp,skip,frm,pbchoice,nginit
      integer hdfform,rprec,iprec,byteo,irun,icbc,entest,ng2d

      real ns,zin,zfin,zpr,hubble,omegadm,omegab
      real rL,rLmin,alpha,deut,qq,ss8,pmass,omegatot

      open(7,file="indat",status='old')
      read(7,*)skip     ! stride in particle data outputs
      read(7,*)jinit    ! ICs (0=>internal, 1=>read, 2=>restart)
      read(7,*)icbc     ! 0=>enforces periodic BCs, otherwise not (jinit=1)
      read(7,*)norm     ! 0=>COBE 4-yr B/W, 1=>sigma8, 2=> COBE Peebles
      read(7,*)trans    ! 0=>T(k)=1, 1=>T(k)=H/S, 2=>T(k)=K/H, 3=>CMBFAST
      read(7,*)initp    ! 0=>part-ICs before BCs,1=>after,2=>dens,else no print
      read(7,*)norder	! order of time integration (2 or 4)
      read(7,*)frm      ! output format, 1=> ascii, 2=> binary
      read(7,*)pbchoice ! units, 1=> box units, 2=> physical units
      read(7,*)entest   ! Layzer-Irvine test 0=> yes, else not
      read(7,*)ng       ! number of grid points in one dimension
      read(7,*)ng2d     ! ng for 2-dim density projection, power of two	
      read(7,*)ns       ! slope of the PS
      read(7,*)iseed	! seed for random # generation
      read(7,*)alpha	! power of scale factor in timestepping
      read(7,*)zin	! initial scale factor 
      read(7,*)zfin     ! final redshift
      read(7,*)zpr      ! redshift after which print-out starts
      read(7,*)np	! number of particles per dimension
      read(7,*)nsteps	! number of time-steps
      read(7,*)iprint   ! output interval
      read(7,*)irun     ! 0=>free-run, else print z determined by zprnt(irun)
      read(7,*)hubble	! hubble constant/100 km/s/Mpc
      read(7,*)omegadm	! dark matter density
      read(7,*)deut	! omegab=deut/h/h
      read(7,*)rL	! physical box size [h^(-1)Mpc)]
      read(7,*)ss8      ! sigma8
      read(7,*)qq       ! 18.e-6=COBE/DMR (4-yr for n=1) Q_rms (K)
      read(7,*)nxint    ! number of evaluation points (ODD) for D_0 integral
      read(7,*)nxdint   ! evaluation points (ODD) for D integral
      read(7,*)hdfform  ! version number of format of hdf 
      read(7,*)rprec    ! precision of the real variables
      read(7,*)iprec    ! precision of the integers
      read(7,*)byteo    ! byteorder
      close(7)

!-----initial conditions grid size, typically locked to np

      nginit=np

!-----Change np to true number of particles

      np=np**3

!-----Particle Mass
 
      omegatot=omegadm+deut/hubble/hubble
      pmass=2.77536627e11*(rL**3)*(omegatot)/hubble/np
 
      write(6,*)'Particle Mass=',pmass

!-----Write input-file for HDF

      rLmin=0.0

      open(14,file='input',form='formatted')

      write(14,*)ng
      write(14,*)hubble
      write(14,*)pmass
      write(14,*)((1.0/(1.0+zin))**alpha)**(1.0/alpha)
      write(14,*)omegadm
      write(14,*)deut
      write(14,*)hdfform
      write(14,*)byteo
      write(14,*)rprec
      write(14,*)iprec
      write(14,*)rL
      write(14,*)rLmin
      write(14,*)ss8
      write(14,*)ns

      close(14)

!-----Brute force (STREAMLINE LATER WITH MODULES)

      call solve(ng,np,nsteps,norder,iprint,iseed,jinit,norm,
     #          trans,nxint,nxdint,ns,zin,hubble,omegadm,rL,alpha,
     #          deut,qq,ss8,initp,skip,frm,pbchoice,zfin,
     #          zpr,irun,icbc,pmass,entest,ng2d,nginit)

      stop
      end




