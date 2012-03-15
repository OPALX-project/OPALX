#!/usr/bin/python
# Filename: runOBLAWithGun.py
# 
# script that launches gun simulation and further simulatin of OBLA phase I setup
# -OPAL creates the OPAL inputfile
import sys,re,os,string
from math import sqrt

##########################################################################
# input parameters with default values                                   #
##########################################################################

N=-1              # a running number; if given use it to label directory!

DT='1.0e-12'      # time step for beam line
DTGUN='1.0e-13'   # time step for gun
NTSTEP='8000'     # number of time steps for beam line
NTSTGUN='2000'    # number of time steps for gun
NBUNCH='39'       # number of energy bins (gun)

GAP='12'          # diode gap in mm  !!!!!!!! not an argument yet !!!!!
VGAP='100'        # diode voltage in V

NPTCL='50000'    # total number of particles

NX='32'           # mesh size in x
NY='32'           # mesh size in y
NZ='32'           # mesh size in z

NSLICES='390'     # number of emission slices
TEMISS='39.e-12'  # emission time, will be computed from SIGT, unless SIGT > 1E-9

SIGT='6.5E-12'    # longitudinal Gauss sigma - if SIGT>1e-9, uniform distribution according to TEMISS
SIGX='0.00030'
SIGY='0.00030'
SIGZ='0.'         # will be computed from SIGT, unless SIGT > 1E-9, then from TEMISS

QBUNCH='6.0E-12'  # bunch charge in C

EKIN='1.0'        # initial kinetic energy in eV
FREQ='1498.956e6' # rf frequency

RESTARTF='1'      # restart flag (should always be 1)
RSTSTEP='1800'    # restart step (typically 1000)

CORES='1'         # number of cores (cpu's) to use

GREENSF='INTEGRATED'
GREENSF='STANDARD'

SP1B='0.'         # SP magnets field strength (in mT)
SP2B='0.'
SP3B='0.'
SP4B='0.'
SP5B='0.'

SP1S='-0.5385'         # SP magnets positions FROM ANODE (in m)
SP2S='-0.4040'         # need to add gap width!
SP3S='-0.2740'
SP4S='-0.1640'
SP5S='-0.0540'

SP1S=str(float(SP1S)+0.001*float(GAP))
SP2S=str(float(SP2S)+0.001*float(GAP))
SP3S=str(float(SP3S)+0.001*float(GAP))
SP4S=str(float(SP4S)+0.001*float(GAP))
SP5S=str(float(SP5S)+0.001*float(GAP))

SP1X='0.'         # SP magnets misalignment (in um)
SP1Y='0.'
SP2X='0.'
SP2Y='0.'
SP3X='0.'
SP3Y='0.'
SP4X='0.'
SP4Y='0.'
SP5X='0.'
SP5Y='0.'

H5FREQ='10'        # h5 dump frequency
H5PSFQ='10'        # h5 phase space dump frequency
H5KEEP='1'         # keep h5 file flag

dirname=''

doOpal=False       # default mode IMPACT-T (we will change that soon!)

##########################################################################
# collect input from command line                                        #
##########################################################################

valre=re.compile('^[a-zA-Z]\w*=-?\d+\.?\d*([eE]-?\d*)?$')
varre=re.compile('^[a-zA-Z]\w*')
numre=re.compile('=(-?\d+\.?\d*([eE]-?\d*)?)$')
for arg in sys.argv[1:]:
    if arg == "-OPAL":
        doOpal=True
    if valre.search(arg):
        var=varre.match(arg)
        num=numre.search(arg)
        if var.group() == "N":
            N=num.group(1)
        if var.group() == "DT":
            DT=num.group(1)
            dirname+='_DT='+DT
        if var.group() == "DTGUN":
            DTGUN=num.group(1)
            dirname+='_DTGUN='+DTGUN
        if var.group() == "NTSTEP":
            NTSTEP=num.group(1)
            dirname+='_NTSTEP='+NTSTEP
        if var.group() == "CORES":
            CORES=num.group(1)
            dirname+='_CORES='+CORES
        if var.group() == "GREENSF":
            GREENSF=num.group(1)
            dirname+='_GREENSF='+GREENSF
        if var.group() == "NTSTGUN":
            NTSTGUN=num.group(1)
            dirname+='_NTSTGUN='+NTSTGUN
        if var.group() == "NBUNCH":
            NBUNCH=num.group(1)
            dirname+='_NBUNCH='+NBUNCH
        if var.group() == "VGAP":
            VGAP=num.group(1)
            dirname+='_VGAP='+VGAP
        if var.group() == "NPTCL":
            NPTCL=num.group(1)
            dirname+='_NPTCL='+NPTCL
        if var.group() == "NX":
            NX=num.group(1)
            dirname+='_NX='+NX
        if var.group() == "NY":
            NY=num.group(1)
            dirname+='_NY='+NY
        if var.group() == "NZ":
            NZ=num.group(1)
            dirname+='_NZ='+NZ
        if var.group() == "NSLICES":
            NSLICES=num.group(1)
            dirname+='_NSLICES='+NSLICES
        if var.group() == "TEMISS":
            TEMISS=num.group(1)
            dirname+='_TEMISS='+TEMISS
        if var.group() == "SIGT":
            SIGT=num.group(1)
            dirname+='_SIGT='+SIGT
        if var.group() == "SIGX":
            SIGX=num.group(1)
            dirname+='_SIGX='+SIGX
        if var.group() == "SIGY":
            SIGY=num.group(1)
            dirname+='_SIGY='+SIGY
        if var.group() == "SIGZ":
            SIGZ=num.group(1)
            dirname+='_SIGZ='+SIGZ
        if var.group() == "QBUNCH":
            QBUNCH=num.group(1)
            dirname+='_QBUNCH='+QBUNCH
        if var.group() == "EKIN":
            EKIN=num.group(1)
            dirname+='_EKIN='+EKIN
        if var.group() == "FREQ":
            FREQ=num.group(1)
            dirname+='_FREQ='+FREQ
        if var.group() == "SP1B":
            SP1B=num.group(1)
            dirname+='_SP1B='+SP1B
        if var.group() == "SP2B":
            SP2B=num.group(1)
            dirname+='_SP2B='+SP2B
        if var.group() == "SP3B":
            SP3B=num.group(1)
            dirname+='_SP3B='+SP3B
        if var.group() == "SP4B":
            SP4B=num.group(1)
            dirname+='_SP4B='+SP4B
        if var.group() == "SP5B":
            SP5B=num.group(1)
            dirname+='_SP5B='+SP5B
        if var.group() == "SP1X":
            SP1X=num.group(1)
            dirname+='_SP1X='+SP1X
        if var.group() == "SP1Y":
            SP1Y=num.group(1)
            dirname+='_SP1Y='+SP1Y
        if var.group() == "SP2X":
            SP2X=num.group(1)
            dirname+='_SP2X='+SP2X
        if var.group() == "SP2Y":
            SP2Y=num.group(1)
            dirname+='_SP2Y='+SP2Y
        if var.group() == "SP3X":
            SP3X=num.group(1)
            dirname+='_SP3X='+SP3X
        if var.group() == "SP3Y":
            SP3Y=num.group(1)
            dirname+='_SP3Y='+SP3Y
        if var.group() == "SP4X":
            SP4X=num.group(1)
            dirname+='_SP4X='+SP4X
        if var.group() == "SP4Y":
            SP4Y=num.group(1)
            dirname+='_SP4Y='+SP4Y
        if var.group() == "SP5X":
            SP5X=num.group(1)
            dirname+='_SP5X='+SP5X
        if var.group() == "SP5Y":
            SP5Y=num.group(1)
            dirname+='_SP5Y='+SP5Y
        if var.group() == "H5FREQ":
            H5FREQ=num.group(1)
        if var.group() == "H5PSFQ":
            H5PSFQ=num.group(1)
        if var.group() == "H5KEEP":
            H5KEEP=num.group(1)
    else:
        print 'bad argument: '+arg+' -- will be ignored!'

# multiply B fields and misalignements by 1.e-6
mu=1.e-6
SP1B=str(float(SP1B)*mu)
SP2B=str(float(SP2B)*mu)
SP3B=str(float(SP3B)*mu)
SP4B=str(float(SP4B)*mu)
SP5B=str(float(SP5B)*mu)

SP1X=str(float(SP1X)*mu)
SP1Y=str(float(SP1Y)*mu)
SP2X=str(float(SP2X)*mu)
SP2Y=str(float(SP2Y)*mu)
SP3X=str(float(SP3X)*mu)
SP3Y=str(float(SP3Y)*mu)
SP4X=str(float(SP4X)*mu)
SP4Y=str(float(SP4Y)*mu)
SP5X=str(float(SP5X)*mu)
SP5Y=str(float(SP5Y)*mu)


me = 510999.  # electron mass in eV
c=2.99e8      # c in m/s


if doOpal == True:
    print 'In OPAL mode'

    if N>=0:
        obladirname='OBLA-OPAL'+str(N)
    else:
        obladirname='OBLA-OPAL'+dirname

# create directories
    os.mkdir(obladirname)
    os.chdir(obladirname)
    os.mkdir('Gun')
    os.chdir('Gun')
    infile='OBLAGun-OPAL.in'
    runscript='run-gun'

    nb = int(NBUNCH)
    qb=float(QBUNCH)
    f = float(FREQ)
    T = float(TEMISS)
    Ek = float(EKIN)

    # compute I for input files:
    I = str(qb*f)
    IBUNCH = str(qb*f)

    gamma = 1.+Ek/me
    beta = sqrt(1.-(1./(gamma*gamma)))
    bg = beta*gamma

    v0 = beta*c

    sigma_t = float(SIGT)

    CENTPZ=str(bg)

    print 'CENTPZ = '+CENTPZ 
    print 'I = '+I 

    CENTZlist=range(1,nb+1)

    DIST='23'

    if sigma_t < 1.e-9:
      # Gaussian distribution:
      # SIGZ follows from SIGT
      # CENTZ is at -3sigma for all bunches
      # TEMISS is 6 sigma 
      sigma_z = sigma_t*v0
      cent_z = -3.*sigma_z
      SIGZ=str(sigma_z)
      for b in range(1,nb+1):
        CENTZlist[b-1]=str(-3.*sigma_z)  
      TEMISS=str(6.*sigma_t)
      print 'GAUSSIAN DISTRIBUTION:'
      print 'quantities derived from SIGT = '+SIGT+':'
      print 'SIGZ = '+SIGZ 
      print 'CENTZ = '+CENTZlist[0]
      print 'TEMISS = '+TEMISS 
    else:
      # Uniform distribution:
      # SIGZ follows from TEMISS given by user
      # CENTZ is a list computed from velocity
      #
      T = float(TEMISS)
      sigma_z = T*v0/2.
      SIGZ=str(sigma_z)
      DIST='17'
      for b in range(1,nb+1):
        s0 = v0*T/float(nb)
        sb = s0*(0.5 + b -1)
        CENTZlist[b-1]=str(-sb)  
      print 'UNIFORM DISTRIBUTION:'
      print 'quantities derived from TEMISS = '+TEMISS+':'
      print 'SIGZ = '+SIGZ 

    # particles per bunch:
    NPBUNCH=str(int(NPTCL))

    # length of diode field map (gap+6 mm)
    LDIODE=str(0.001*(float(GAP)+6.))

    # voltage scale factor (depends on gap width
    # and gap voltage):
    v=float(VGAP)
    VSCALE='0'
    
    if GAP=='4':
        VSCALE=str(-0.0002*v)    
    # For 12mm and 500kV the peal field is -46.84 MV/m 
    if GAP=='12':
        VSCALE=str(-46.84)     

    # position for merging bins (10 mm after anode):
    #SMERGE=str(0.01*(float(GAP)+10.))
    SMERGE=str(1000.0)
    os.system('cp /scratch0/shared/adelmann/Gun/Impact-t/input/OBLAGun-OPAL.in '+infile)
    os.system('cp /scratch0/shared/adelmann/Gun/Impact-t/input/run-gun '+runscript)

    os.system('sed -i \'s/_DTGUN_/'+DTGUN+'/\' '+infile)
    os.system('sed -i \'s/_NTSTGUN_/'+NTSTGUN+'/\' '+infile)
    os.system('sed -i \'s/_NBUNCH_/'+NBUNCH+'/\' '+infile)

    os.system('sed -i \'s/_NPBUNCH_/'+NPBUNCH+'/\' '+infile)

    os.system('sed -i \'s/_NX_/'+NX+'/\' '+infile)
    os.system('sed -i \'s/_NY_/'+NY+'/\' '+infile)
    os.system('sed -i \'s/_NZ_/'+NZ+'/\' '+infile)

    os.system('sed -i \'s/_DIST_/'+DIST+'/\' '+infile)
    os.system('sed -i \'s/_NSLICES_/'+NSLICES+'/\' '+infile)
    os.system('sed -i \'s/_TEMISS_/'+TEMISS+'/\' '+infile)
    os.system('sed -i \'s/_SIGX_/'+SIGX+'/\' '+infile)
    os.system('sed -i \'s/_SIGY_/'+SIGY+'/\' '+infile)
    os.system('sed -i \'s/_SIGT_/'+SIGT+'/\' '+infile)

    os.system('sed -i \'s/_IBUNCH_/'+IBUNCH+'/\' '+infile)
    os.system('sed -i \'s/_QBUNCH_/'+QBUNCH+'/\' '+infile)
    os.system('sed -i \'s/_EKIN_/'+EKIN+'/\' '+infile)
    os.system('sed -i \'s/_FREQ_/'+FREQ+'/\' '+infile)

    os.system('sed -i \'s/_SMERGE_/'+SMERGE+'/\' '+infile)
    os.system('sed -i \'s/_LDIODE_/'+LDIODE+'/\' '+infile)
    os.system('sed -i \'s/_VSCALE_/'+VSCALE+'/\' '+infile)

    os.system('sed -i \'s/_SP1B_/'+SP1B+'/\' '+infile)
    os.system('sed -i \'s/_SP1S_/'+SP1S+'/\' '+infile)
    os.system('sed -i \'s/_SP2B_/'+SP2B+'/\' '+infile)
    os.system('sed -i \'s/_SP2S_/'+SP2S+'/\' '+infile)
    os.system('sed -i \'s/_SP3B_/'+SP3B+'/\' '+infile)
    os.system('sed -i \'s/_SP3S_/'+SP3S+'/\' '+infile)

    os.system('sed -i \'s/_SP1X_/'+SP1X+'/\' '+infile)
    os.system('sed -i \'s/_SP1Y_/'+SP1Y+'/\' '+infile)

    os.system('sed -i \'s/_GREENSF_/'+GREENSF+'/\' '+infile)

    print 'wrote '+obladirname+'/Gun/'+infile

    os.system('ln -s /home5/schietinger/ImpactT/OBLA/magnet/1T1.T7 1T1.T7')
    os.system('ln -s /home5/schietinger/ImpactT/OBLA/magnet/1T2.T7 1T2.T7')
    os.system('ln -s /home5/schietinger/ImpactT/OBLA/rf/GAP'+GAP+'mm.T7 1T3.T7')

    print 'linked magnet and rf files - now starting the gun simulation'

    ### run the simulation, first argument is the number of processors to use!
    command = './run-gun '+CORES
    print command
    os.system(command)

    print 'done with OPAL'

else:
    print 'In IMPACT-T mode'

    if N>=0:
        obladirname='OBLA-'+str(N)
    else:
        obladirname='OBLA'+dirname

# create directories
    os.mkdir(obladirname)
    os.chdir(obladirname)
    os.mkdir('Gun')


    os.chdir('Gun')
    infile0='ImpactT'

    nb = int(NBUNCH)
    qb=float(QBUNCH)
    f = float(FREQ)
    T = float(TEMISS)
    Ek = float(EKIN)

    # compute I for input files:
    I = str(qb*f)
    IBUNCH = str(qb*f/float(nb))

    gamma = 1.+Ek/me
    beta = sqrt(1.-(1./(gamma*gamma)))
    bg = beta*gamma

    v0 = beta*c

    sigma_t = float(SIGT)

    CENTPZ=str(bg)

    print 'CENTPZ = '+CENTPZ 
    print 'I = '+I 

    CENTZlist=range(1,nb+1)

    DIST='23'

    if sigma_t < 1.e-9:
      # Gaussian distribution:
      # SIGZ follows from SIGT
      # CENTZ is at -3sigma for all bunches
      # TEMISS is 6 sigma 
      sigma_z = sigma_t*v0
      cent_z = -3.*sigma_z
      SIGZ=str(sigma_z)
      for b in range(1,nb+1):
        CENTZlist[b-1]=str(-3.*sigma_z)  
      TEMISS=str(6.*sigma_t)
      print 'GAUSSIAN DISTRIBUTION:'
      print 'quantities derived from SIGT = '+SIGT+':'
      print 'SIGZ = '+SIGZ 
      print 'CENTZ = '+CENTZlist[0]
      print 'TEMISS = '+TEMISS 
    else:
      # Uniform distribution:
      # SIGZ follows from TEMISS given by user
      # CENTZ is a list computed from velocity
      #
      T = float(TEMISS)
      sigma_z = T*v0/2.
      SIGZ=str(sigma_z)
      DIST='17'
      for b in range(1,nb+1):
        s0 = v0*T/float(nb)
        sb = s0*(0.5 + b -1)
        CENTZlist[b-1]=str(-sb)  
      print 'UNIFORM DISTRIBUTION:'
      print 'quantities derived from TEMISS = '+TEMISS+':'
      print 'SIGZ = '+SIGZ 

    print 'CENTZ list:'
    print CENTZlist

    # particles per bunch:
    NPBUNCH=str(int(NPTCL)/nb)

    # length of diode field map (gap+6 mm)
    LDIODE=str(0.001*(float(GAP)+6.))

    # voltage scale factor (depends on gap width
    # and gap voltage):
    v=float(VGAP)
    VSCALE='0'
    if GAP=='4':
      VSCALE=str(-0.0002*v)    
    if GAP=='12':
      VSCALE=str(-0.00000019976*v)

    # position for merging bins (10 mm after anode):
    #SMERGE=str(0.01*(float(GAP)+10.))
    SMERGE=str(1000.)

    for i in range(1,nb+1):
      infile = infile0  
      if i!=1:
        infile+=str(i)
      infile+='.in'

      CENTZ=CENTZlist[i-1]

      os.system('cp /home5/schietinger/ImpactT/OBLA/input/OBLAGun.in '+infile)

      os.system('sed -i \'s/_DTGUN_/'+DTGUN+'/\' '+infile)
      os.system('sed -i \'s/_NTSTGUN_/'+NTSTGUN+'/\' '+infile)
      os.system('sed -i \'s/_NBUNCH_/'+NBUNCH+'/\' '+infile)

      os.system('sed -i \'s/_NPBUNCH_/'+NPBUNCH+'/\' '+infile)

      os.system('sed -i \'s/_NX_/'+NX+'/\' '+infile)
      os.system('sed -i \'s/_NY_/'+NY+'/\' '+infile)
      os.system('sed -i \'s/_NZ_/'+NZ+'/\' '+infile)

      os.system('sed -i \'s/_DIST_/'+DIST+'/\' '+infile)
      os.system('sed -i \'s/_NSLICES_/'+NSLICES+'/\' '+infile)
      os.system('sed -i \'s/_TEMISS_/'+TEMISS+'/\' '+infile)
      os.system('sed -i \'s/_SIGX_/'+SIGX+'/\' '+infile)
      os.system('sed -i \'s/_SIGY_/'+SIGY+'/\' '+infile)
      os.system('sed -i \'s/_SIGZ_/'+SIGZ+'/\' '+infile)
      os.system('sed -i \'s/_CENTZ_/'+CENTZ+'/\' '+infile)
      os.system('sed -i \'s/_CENTPZ_/'+CENTPZ+'/\' '+infile)

      os.system('sed -i \'s/_IBUNCH_/'+IBUNCH+'/\' '+infile)
      os.system('sed -i \'s/_EKIN_/'+EKIN+'/\' '+infile)
      os.system('sed -i \'s/_FREQ_/'+FREQ+'/\' '+infile)

      os.system('sed -i \'s/_SMERGE_/'+SMERGE+'/\' '+infile)
      os.system('sed -i \'s/_LDIODE_/'+LDIODE+'/\' '+infile)
      os.system('sed -i \'s/_VSCALE_/'+VSCALE+'/\' '+infile)

      os.system('sed -i \'s/_SP1B_/'+SP1B+'/\' '+infile)
      os.system('sed -i \'s/_SP1S_/'+SP1S+'/\' '+infile)

      os.system('sed -i \'s/_SP1X_/'+SP1X+'/\' '+infile)
      os.system('sed -i \'s/_SP1Y_/'+SP1Y+'/\' '+infile)

      print 'wrote '+obladirname+'/Gun/'+infile

    else:
      print 'finished writing files'  

    os.system('ln -s /home5/schietinger/ImpactT/OBLA/magnet/1T1.T7 1T1.T7')
    os.system('ln -s /home5/schietinger/ImpactT/OBLA/magnet/1T2.T7 1T2.T7')
    #os.system('ln -s /home5/schietinger/ImpactT/OBLA/rf/rfdata1 rfdata1')
    os.system('ln -s /home5/schietinger/ImpactT/OBLA/rf/GAP'+GAP+'mm.T7 1T3.T7')

    print 'linked magnet and rf files - now starting the gun simulation'

    ### run the simulation!
    command = '/home5/schietinger/bin/impactt ImpactT.in '+H5FREQ+' '+H5PSFQ+' | tee sim.out'
    print command
    os.system(command)

    ### remove h5 file unless it is requested
    if (H5KEEP=='0') or (H5KEEP=='NO') or (H5KEEP=='no'):
        print 'removing h5 file'
        os.system('rm -f *.h5')

    print 'done simulating gun'

    ##########################################################################
    # simulate the rest                                                      #
    ##########################################################################

    os.chdir('..')

    infile='OBLAPostGun.in'

    os.system('cp /home5/schietinger/ImpactT/OBLA/input/OBLAPostGun.in '+infile)

    os.system('sed -i \'s/_DT_/'+DT+'/\' '+infile)
    os.system('sed -i \'s/_NTSTEP_/'+NTSTEP+'/\' '+infile)

    os.system('sed -i \'s/_NPTCL_/'+NPTCL+'/\' '+infile)

    os.system('sed -i \'s/_NX_/'+NX+'/\' '+infile)
    os.system('sed -i \'s/_NY_/'+NY+'/\' '+infile)
    os.system('sed -i \'s/_NZ_/'+NZ+'/\' '+infile)

    os.system('sed -i \'s/_RESTARTF_/'+RESTARTF+'/\' '+infile)
    os.system('sed -i \'s/_RSTSTEP_/'+RSTSTEP+'/\' '+infile)

    os.system('sed -i \'s/_I_/'+I+'/\' '+infile)
    os.system('sed -i \'s/_EKIN_/'+EKIN+'/\' '+infile)
    os.system('sed -i \'s/_FREQ_/'+FREQ+'/\' '+infile)

    os.system('sed -i \'s/_SP1B_/'+SP1B+'/\' '+infile)
    os.system('sed -i \'s/_SP2B_/'+SP2B+'/\' '+infile)
    os.system('sed -i \'s/_SP3B_/'+SP3B+'/\' '+infile)
    os.system('sed -i \'s/_SP4B_/'+SP4B+'/\' '+infile)
    os.system('sed -i \'s/_SP5B_/'+SP5B+'/\' '+infile)

    os.system('sed -i \'s/_SP1S_/'+SP1S+'/\' '+infile)
    os.system('sed -i \'s/_SP2S_/'+SP2S+'/\' '+infile)
    os.system('sed -i \'s/_SP3S_/'+SP3S+'/\' '+infile)
    os.system('sed -i \'s/_SP4S_/'+SP4S+'/\' '+infile)
    os.system('sed -i \'s/_SP5S_/'+SP5S+'/\' '+infile)

    os.system('sed -i \'s/_SP1X_/'+SP1X+'/\' '+infile)
    os.system('sed -i \'s/_SP1Y_/'+SP1Y+'/\' '+infile)
    os.system('sed -i \'s/_SP2X_/'+SP2X+'/\' '+infile)
    os.system('sed -i \'s/_SP2Y_/'+SP2Y+'/\' '+infile)
    os.system('sed -i \'s/_SP3X_/'+SP3X+'/\' '+infile)
    os.system('sed -i \'s/_SP3Y_/'+SP3Y+'/\' '+infile)
    os.system('sed -i \'s/_SP4X_/'+SP4X+'/\' '+infile)
    os.system('sed -i \'s/_SP4Y_/'+SP4Y+'/\' '+infile)
    os.system('sed -i \'s/_SP5X_/'+SP5X+'/\' '+infile)
    os.system('sed -i \'s/_SP5Y_/'+SP5Y+'/\' '+infile)

    print 'wrote '+obladirname+'/'+infile

    print 'linking restart file from gun'

    # check this first...
    linkcmd='ln -s Gun/ImpactT000001.h5 OBLAPostGun000001.h5'
    print linkcmd
    os.system(linkcmd)

    # old days: copy restart file
    #sp1f=float(SP1B)
    #sp1str = '%(x)03d' % {"x": (int(round(sp1f*1.e6)))}
    #restartfile='restart-'+sp1str+'.dat'
    #os.system('ln -s /home5/schietinger/ImpactT/OBLA/restart/'+restartfile+' restart.dat')

    #print 'copied '+restartfile

    os.system('ln -s /home5/schietinger/ImpactT/OBLA/magnet/1T1.T7 1T1.T7')
    os.system('ln -s /home5/schietinger/ImpactT/OBLA/magnet/1T2.T7 1T2.T7')

    print 'linked magnet files - now starting the simulation'

    ### run the simulation!
    command='/home5/schietinger/bin/impactt '+infile+' '+H5FREQ+' '+H5PSFQ+' | tee sim.out'
    print command
    os.system(command)

    ### remove h5 file unless it is requested
    if (H5KEEP=='0') or (H5KEEP=='NO') or (H5KEEP=='no'):
        print 'removing h5 file'
        os.system('rm -f *.h5')

    print 'done with IMPACT-T!'
