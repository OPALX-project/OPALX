#!/usr/bin/python
# Filename: runOBLA2WithGun.py
# 
# script that launches gun simulation and further simulatin of OBLA phase II setup
# -OPAL creates the OPAL inputfile
import sys,re,os,string
from math import sqrt

##########################################################################
# input parameters with default values                                   #
##########################################################################

N=-1              # a running number; if given use it to label directory!

DT='1.0e-12'      # time step for beam line
DTGUN='1.0e-13'   # time step for gun
NTSTEP='12000'    # number of time steps for beam line
NTSTGUN='2000'    # number of time steps for gun
NBUNCH='39'       # number of energy bins (gun)

GAP='5'           # diode gap in mm (only integer values between 4 and 10)
VGAP='300'        # diode voltage in kV

NPTCL='50000'     # total number of particles

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
RSTSTEP='1800'    # restart step (typically 1800)

CORES='1'         # number of cores (cpu's) to use

GREENSF='INTEGRATED'
GREENSF='STANDARD'

PSLB='20.'         # solenoid magnets field strength (in mT)
SL10B='85.'        # PSL   = pulsed solenoid (FINEG-MSL10)
SL20B='50.'        # SL10B = FINLB01-MSL10 etc.
                   # default values are Kevin's values

PSLS='-0.004'      # solenoid magnets positions FROM ANODE (in m)
SL10S='0.33'       # need to add gap width!
SL11S='0.43'       # (SL11 is second field for SL10)
SL20S='0.83'
SL21S='0.93'       # (SL21 is second field for SL20)

PSLX='0.'         # solenoid magnets misalignment (in um)
PSLY='0.'
SL10X='0.'
SL10Y='0.'
SL11X='0.'
SL11Y='0.'
SL20X='0.'
SL20Y='0.'
SL21X='0.'
SL21Y='0.'

RACS='0.099'
PHASE='70.'

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
        if var.group() == "GAP":
            GAP=num.group(1)
            dirname+='_GAP='+GAP
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
        if var.group() == "PSLB":
            PSLB=num.group(1)
            dirname+='_PSLB='+PSLB
        if var.group() == "SL10B":
            SL10B=num.group(1)
            dirname+='_SL10B='+SL10B
        if var.group() == "SL20B":
            SL20B=num.group(1)
            dirname+='_SL20B='+SL20B
        if var.group() == "PSLX":
            PSLX=num.group(1)
            dirname+='_PSLX='+PSLX
        if var.group() == "PSLY":
            PSLY=num.group(1)
            dirname+='_PSLY='+PSLY
        if var.group() == "SL10X":
            SL10X=num.group(1)
            dirname+='_SL10X='+SL10X
        if var.group() == "SL10Y":
            SL10Y=num.group(1)
            dirname+='_SL10Y='+SL10Y
        if var.group() == "SL20X":
            SL20X=num.group(1)
            dirname+='_SL20X='+SL20X
        if var.group() == "SL20Y":
            SL20Y=num.group(1)
            dirname+='_SL20Y='+SL20Y
        if var.group() == "PHASE":
            PHASE=num.group(1)
            dirname+='_PHASE='+PHASE
        if var.group() == "RACS":
            RACS=num.group(1)
            dirname+='_RACS='+PHASE
        if var.group() == "H5FREQ":
            H5FREQ=num.group(1)
        if var.group() == "H5PSFQ":
            H5PSFQ=num.group(1)
        if var.group() == "H5KEEP":
            H5KEEP=num.group(1)
    else:
        print 'bad argument: '+arg+' -- will be ignored!'

# add gap to magnet positions
PSLS=str(float(PSLS)+0.001*float(GAP))
SL10S=str(float(SL10S)+0.001*float(GAP))
SL11S=str(float(SL11S)+0.001*float(GAP))
SL20S=str(float(SL20S)+0.001*float(GAP))
SL21S=str(float(SL21S)+0.001*float(GAP))

# multiply B fields and misalignements by 1.e-6
mu=1.e-6
m=1.e-3
PSLB=str(float(PSLB)*m) 
SL10B=str(float(SL10B)*m)
SL20B=str(float(SL20B)*m)
SL11B=str(-float(SL10B))
SL21B=str(-float(SL20B))

PSLX=str(float(PSLX)*mu)
PSLY=str(float(PSLY)*mu)
SL10X=str(float(SL10X)*mu)
SL10Y=str(float(SL10Y)*mu)
SL20X=str(float(SL20X)*mu)
SL20Y=str(float(SL20Y)*mu)
SL11X=str(float(SL10X))
SL11Y=str(float(SL10Y))
SL21X=str(float(SL20X))
SL21Y=str(float(SL20Y))

# add gap to cavity position
RACS=str(float(RACS)+0.001*float(GAP))

me = 510999.  # electron mass in eV
c=2.99e8      # c in m/s


if doOpal == True:
    print 'In OPAL mode'

    if N>=0:
        obladirname='OBLA2-OPAL'+str(N)
    else:
        obladirname='OBLA2-OPAL'+dirname

# create directories
    os.mkdir(obladirname)
    os.chdir(obladirname)
    os.mkdir('Gun')
    os.chdir('Gun')
    infile='OBLA2Gun-OPAL.in'
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
    # adjust voltage to yield correct energy (empirically...)
    # (factor is 1000 larger than in IMPACTT - why?)
    VSCALE=str(-1.995e-4*v)
    
    # position for merging bins (10 mm after anode):
    #SMERGE=str(0.01*(float(GAP)+10.))
    #SMERGE=str(1000.0)
    # OPAL uses DEBIN (set to 1 keV in OBLA2Gun-OPAL.in)
    
    os.system('cp $OPAL_ROOT/opal-Tests/RegressionTests/Gun/Impact-t/input/OBLA2Gun-OPAL.in '+infile)
    os.system('cp $OPAL_ROOT/opal-Tests/RegressionTests/Gun/Impact-t/input/run-gun '+runscript)

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

    #os.system('sed -i \'s/_SMERGE_/'+SMERGE+'/\' '+infile)
    os.system('sed -i \'s/_LDIODE_/'+LDIODE+'/\' '+infile)
    os.system('sed -i \'s/_VSCALE_/'+VSCALE+'/\' '+infile)

    os.system('sed -i \'s/_PSLS_/'+PSLS+'/\' '+infile)
    os.system('sed -i \'s/_SL10S_/'+SL10S+'/\' '+infile)
    os.system('sed -i \'s/_SL20S_/'+SL20S+'/\' '+infile)
    os.system('sed -i \'s/_SL11S_/'+SL11S+'/\' '+infile)
    os.system('sed -i \'s/_SL21S_/'+SL21S+'/\' '+infile)
    os.system('sed -i \'s/_PSLB_/'+PSLB+'/\' '+infile)
    os.system('sed -i \'s/_SL10B_/'+SL10B+'/\' '+infile)
    os.system('sed -i \'s/_SL20B_/'+SL20B+'/\' '+infile)
    os.system('sed -i \'s/_SL11B_/'+SL11B+'/\' '+infile)
    os.system('sed -i \'s/_SL21B_/'+SL21B+'/\' '+infile)

    os.system('sed -i \'s/_PSLX_/'+PSLX+'/\' '+infile)
    os.system('sed -i \'s/_PSLY_/'+PSLY+'/\' '+infile)
    os.system('sed -i \'s/_SL10X_/'+SL10X+'/\' '+infile)
    os.system('sed -i \'s/_SL10Y_/'+SL10Y+'/\' '+infile)
    os.system('sed -i \'s/_SL20X_/'+SL20X+'/\' '+infile)
    os.system('sed -i \'s/_SL20Y_/'+SL20Y+'/\' '+infile)
    os.system('sed -i \'s/_SL11X_/'+SL11X+'/\' '+infile)
    os.system('sed -i \'s/_SL11Y_/'+SL11Y+'/\' '+infile)
    os.system('sed -i \'s/_SL21X_/'+SL21X+'/\' '+infile)
    os.system('sed -i \'s/_SL21Y_/'+SL21Y+'/\' '+infile)

    os.system('sed -i \'s/_PHASE_/'+PHASE+'/\' '+infile)
    os.system('sed -i \'s/_RACS_/'+RACS+'/\' '+infile)

    os.system('sed -i \'s/_GREENSF_/'+GREENSF+'/\' '+infile)

    print 'wrote '+obladirname+'/Gun/'+infile

    os.system('ln -s /home5/schietinger/ImpactT/OBLA/rf/OPAL/DIODE_'+GAP+'_MM.T7 1T1.T7')
    os.system('ln -s /home5/schietinger/ImpactT/OBLA/magnet/OPAL/FINEG-MSL10.T7 1T2.T7')
    os.system('ln -s /home5/schietinger/ImpactT/OBLA/magnet/OPAL/FINLB01-RACF.T7 1T3.T7')
    os.system('ln -s /home5/schietinger/ImpactT/OBLA/magnet/OPAL/FINLB01-RACH.T7 1T4.T7')
    os.system('ln -s /home5/schietinger/ImpactT/OBLA/magnet/OPAL/FINLB01-MSL10.T7 1T5.T7')

    print 'linked magnet and rf files - now starting the gun simulation'

    ### run the simulation, first argument is the number of processors to use!
    command = './run-gun '+CORES
    print command
    os.system(command)

    print 'done with OPAL'

else:
    print 'In IMPACT-T mode'

    if N>=0:
        obladirname='OBLA2-'+str(N)
    else:
        obladirname='OBLA2'+dirname

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
    g=float(GAP)
    # adjust voltage to yield correct energy (empirically...)
    VSCALE=str(-1.9996e-7*v)
    
    # position for merging bins (10 mm after anode):
    SMERGE=str(0.001*(g+10.))

    for i in range(1,nb+1):
        infile = infile0  
        if i!=1:
            infile+=str(i)
        infile+='.in'

        CENTZ=CENTZlist[i-1]
        
        os.system('cp $OPAL_ROOT/opal-Tests/RegressionTests/Gun/Impact-t/input/OBLA2Gun.in '+infile)
        
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

        os.system('sed -i \'s/_PSLS_/'+PSLS+'/\' '+infile)
        os.system('sed -i \'s/_SL10S_/'+SL10S+'/\' '+infile)
        os.system('sed -i \'s/_SL20S_/'+SL20S+'/\' '+infile)
        os.system('sed -i \'s/_SL11S_/'+SL11S+'/\' '+infile)
        os.system('sed -i \'s/_SL21S_/'+SL21S+'/\' '+infile)
        os.system('sed -i \'s/_PSLB_/'+PSLB+'/\' '+infile)
        os.system('sed -i \'s/_SL10B_/'+SL10B+'/\' '+infile)
        os.system('sed -i \'s/_SL20B_/'+SL20B+'/\' '+infile)
        os.system('sed -i \'s/_SL11B_/'+SL11B+'/\' '+infile)
        os.system('sed -i \'s/_SL21B_/'+SL21B+'/\' '+infile)

        os.system('sed -i \'s/_PSLX_/'+PSLX+'/\' '+infile)
        os.system('sed -i \'s/_PSLY_/'+PSLY+'/\' '+infile)
        os.system('sed -i \'s/_SL10X_/'+SL10X+'/\' '+infile)
        os.system('sed -i \'s/_SL10Y_/'+SL10Y+'/\' '+infile)
        os.system('sed -i \'s/_SL20X_/'+SL20X+'/\' '+infile)
        os.system('sed -i \'s/_SL20Y_/'+SL20Y+'/\' '+infile)
        os.system('sed -i \'s/_SL11X_/'+SL11X+'/\' '+infile)
        os.system('sed -i \'s/_SL11Y_/'+SL11Y+'/\' '+infile)
        os.system('sed -i \'s/_SL21X_/'+SL21X+'/\' '+infile)
        os.system('sed -i \'s/_SL21Y_/'+SL21Y+'/\' '+infile)

        os.system('sed -i \'s/_PHASE_/'+PHASE+'/\' '+infile)
        os.system('sed -i \'s/_RACS_/'+RACS+'/\' '+infile)
    
        print 'wrote '+obladirname+'/Gun/'+infile
        
    else:
        print 'finished writing files'  

    os.system('ln -s /home5/schietinger/ImpactT/OBLA/rf/DIODE_'+GAP+'_MM.T7 1T1.T7')
    os.system('ln -s /home5/schietinger/ImpactT/OBLA/magnet/FINEG-MSL10.T7 1T2.T7')
    os.system('ln -s /home5/schietinger/ImpactT/OBLA/magnet/FINLB01-RACF.T7 1T3.T7')
    os.system('ln -s /home5/schietinger/ImpactT/OBLA/magnet/FINLB01-RACH.T7 1T4.T7')
    os.system('ln -s /home5/schietinger/ImpactT/OBLA/magnet/FINLB01-MSL10.T7 1T5.T7')

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

    infile='OBLA2PostGun.in'

    os.system('cp $OPAL_ROOT/opal-Tests/RegressionTests/Gun/Impact-t/input/OBLA2PostGun.in '+infile)

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

    os.system('sed -i \'s/_PSLS_/'+PSLS+'/\' '+infile)
    os.system('sed -i \'s/_SL10S_/'+SL10S+'/\' '+infile)
    os.system('sed -i \'s/_SL20S_/'+SL20S+'/\' '+infile)
    os.system('sed -i \'s/_SL11S_/'+SL11S+'/\' '+infile)
    os.system('sed -i \'s/_SL21S_/'+SL21S+'/\' '+infile)
    os.system('sed -i \'s/_PSLB_/'+PSLB+'/\' '+infile)
    os.system('sed -i \'s/_SL10B_/'+SL10B+'/\' '+infile)
    os.system('sed -i \'s/_SL20B_/'+SL20B+'/\' '+infile)
    os.system('sed -i \'s/_SL11B_/'+SL11B+'/\' '+infile)
    os.system('sed -i \'s/_SL21B_/'+SL21B+'/\' '+infile)

    os.system('sed -i \'s/_PSLX_/'+PSLX+'/\' '+infile)
    os.system('sed -i \'s/_PSLY_/'+PSLY+'/\' '+infile)
    os.system('sed -i \'s/_SL10X_/'+SL10X+'/\' '+infile)
    os.system('sed -i \'s/_SL10Y_/'+SL10Y+'/\' '+infile)
    os.system('sed -i \'s/_SL20X_/'+SL20X+'/\' '+infile)
    os.system('sed -i \'s/_SL20Y_/'+SL20Y+'/\' '+infile)
    os.system('sed -i \'s/_SL11X_/'+SL11X+'/\' '+infile)
    os.system('sed -i \'s/_SL11Y_/'+SL11Y+'/\' '+infile)
    os.system('sed -i \'s/_SL21X_/'+SL21X+'/\' '+infile)
    os.system('sed -i \'s/_SL21Y_/'+SL21Y+'/\' '+infile)

    os.system('sed -i \'s/_PHASE_/'+PHASE+'/\' '+infile)
    os.system('sed -i \'s/_RACS_/'+RACS+'/\' '+infile)

    print 'wrote '+obladirname+'/'+infile

    print 'linking restart file from gun'

    # check this first...
    linkcmd='ln -s Gun/ImpactT000001.h5 OBLA2PostGun000001.h5'
    print linkcmd
    os.system(linkcmd)

    # old days: copy restart file
    #sp1f=float(SP1B)
    #sp1str = '%(x)03d' % {"x": (int(round(sp1f*1.e6)))}
    #restartfile='restart-'+sp1str+'.dat'
    #os.system('ln -s /home5/schietinger/ImpactT/OBLA/restart/'+restartfile+' restart.dat')

    #print 'copied '+restartfile

    os.system('ln -s /home5/schietinger/ImpactT/OBLA/magnet/FINEG-MSL10.T7 1T2.T7')
    os.system('ln -s /home5/schietinger/ImpactT/OBLA/magnet/FINLB01-RACF.T7 1T3.T7')
    os.system('ln -s /home5/schietinger/ImpactT/OBLA/magnet/FINLB01-RACH.T7 1T4.T7')
    os.system('ln -s /home5/schietinger/ImpactT/OBLA/magnet/FINLB01-MSL10.T7 1T5.T7')

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
