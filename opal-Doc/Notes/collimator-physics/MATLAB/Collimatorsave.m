clc
clear

%constant-------------------
c=3e8;
h=6.626e-34;%J s
epsilon0=8.854e-12;%F m-1
me=0.510998918e6;%eV/c2
Avo=6.022e23;%mol-1
re=2.82e-15;%m

%Particle proton-----------------
pmass=0.938272e9;%eV/c^2
zin=1;
    
%Material Cu------------------------
%Energy range 0 to 72MeV
%Cu63 75% Cu65 25%
%*Notice: the parameters below is only valid for Cu*
Z=29;
A=63.546;
rho=8.96;%g/cm^3

X0=12.86/rho/100;%m
I=10*Z;%eV
n=rho/A*Avo;%Atoms/cm^3


%particle number and steps-------------
mparticle=100;
step=1500;
deltat=1e-12;%s
%when a particle is lost, loss(i)=1. rec is used for record particle at a
%fixed position
loss=zeros(1,mparticle);
rec=zeros(1,mparticle);


%Initialize=========================================
%xm records the x value of particles after each step

xm=zeros(1,mparticle);
pxm=zeros(1,mparticle);
ym=zeros(1,mparticle);
pym=zeros(1,mparticle);
zm=zeros(1,mparticle);
pzm=zeros(1,mparticle);
Em=zeros(1,mparticle);


for j=1:step
for i=1:mparticle
      if loss(i)==1
          continue
      end
    if j==1
    %Initial Beam
    x=normrnd(0,0.005);%m
    px=0;
    y=normrnd(0,0.005);%m
    py=0;
    z=0;
    pz=sqrt(1-px^2-py^2);%rad   
     E=72;
    p=sqrt((pmass+E*1e6)^2-pmass^2);%eV/c total
    gamma=(E*1e6+pmass)/pmass;
    beta=sqrt(1-(1/gamma^2));
     
    ptot=sqrt(px^2+py^2+pz^2); %rad
    
    deltas=deltat*beta*c;%m
    deltasrho=deltas*100*rho;%g/cm^2
  
    else
    x=xm(i);
    px=pxm(i);
    y=ym(i);
    py=pym(i);
    z=zm(i);
    pz=pzm(i);
    E=Em(i);

    p=sqrt((pmass+E*1e6)^2-(pmass)^2); %eV/c  
    gamma=(E*1e6+pmass)/pmass;
    beta=sqrt(1-(1/gamma^2)); 
    
    ptot=sqrt(px^2+py^2+pz^2);
    
    deltas=deltat*beta*c;%m
    deltasrho=deltas*100*rho;%g/cm^2

    end



%ellipse collimator
if(z>0.01&&z<0.1&&(x^2+y^2)<0.003^2&&sqrt(0.003^2-(x^2+y^2))<deltas)
    deltas=deltas/10;
       deltasrho=deltas*100*rho;%g/cm^2
end
if (z>0.01&&z<0.1&&(x^2+y^2)>0.003^2&&sqrt((x^2+y^2)-0.003^2)>1e-5)
    if (sqrt((x^2+y^2)-0.003^2)<deltas)
         deltas=deltas/10;
       deltasrho=deltas*100*rho;%g/cm^2
    end



    
%energy loss=================================================== 
K=0.307075;%MeV cm^2 mol-1
Tmax=2*me*beta^2*gamma^2/(1+2*gamma*me/pmass+(me/pmass)^2);%eV
dEdx=-K*zin^2*Z/(A*beta^2)*(1/2*log(2*me*beta^2*gamma^2*Tmax/I^2)-beta^2);%MeV g-1 cm^2
%Compared with the NIST table
%No need to consider density effect since it is not ultrarelativistic particles
delta_Eave=deltasrho*dEdx;
sigma_E=sqrt(4*pi*Avo*re^2*(me/1e6)^2*rho*Z/A*deltas*1e6); %MeV
delta_E=normrnd(delta_Eave,sigma_E);
E=E+delta_E;%MeV

if E<0.1
    loss(i)=1;
    continue
end
p=sqrt((pmass+E*1e6)^2-(pmass)^2); %eV/c
    gamma=(E*1e6+pmass)/pmass;
    beta=sqrt(1-(1/gamma^2));
ptot=1;
pz=sqrt(ptot^2-px^2-py^2);%rad


%multiple Coulomb scattering===================================
theta0=13.6e6/(beta*p)*zin*sqrt(deltas/X0)*(1+0.038*log(deltas/X0)); %rad *Notice log*
%x direction------------------------------
z1=normrnd(0,1);
z2=normrnd(0,1);
thetaplanex=z2*theta0;
while abs(thetaplanex)>3.5*theta0
    z1=normrnd(0,1);
    z2=normrnd(0,1);
    thetaplanex=z2*theta0;
end
xplane=z1*deltas*theta0/sqrt(12)+z2*deltas*theta0/2;
 [rotang,pxpar,pzpar] = rot(px,pz,thetaplanex);

 
 x=x+deltas*px/ptot+xplane*cos(rotang);
 z=z-xplane*sin(rotang);
 px=pxpar*cos(rotang)+pzpar*sin(rotang);
 pz=-pxpar*sin(rotang)+pzpar*cos(rotang);

P2=rand(1);
if P2<0.0047;
    pop=1;
    P3=rand(1);

     thetarux=2.5*sqrt(1/P3)*sqrt(2)*theta0;  
     P4=rand(1);    
    if P4>0.5
    thetarux=thetarux;
    else
    thetarux=-thetarux;
    end
    [rotang,pxpar,pzpar] = rot(px,pz,thetarux);

 
    
 px=pxpar*cos(rotang)+pzpar*sin(rotang);
 pz=-pxpar*sin(rotang)+pzpar*cos(rotang);
 end
%y direction------------------------------
z1=normrnd(0,1);
z2=normrnd(0,1);
thetaplaney=z2*theta0;
while abs(thetaplaney)>3.5*theta0
    z1=normrnd(0,1);
    z2=normrnd(0,1);
    thetaplaney=z2*theta0;
end
yplane=z1*deltas*theta0/sqrt(12)+z2*deltas*theta0/2;
 [rotang,pypar,pzpar] = rot(py,pz,thetaplaney);

 
 y=y+deltas*py/ptot+yplane*cos(rotang);
 z=z+deltas*pz/ptot-yplane*sin(rotang);
 py=pypar*cos(rotang)+pzpar*sin(rotang);
 pz=-pypar*sin(rotang)+pzpar*cos(rotang);

P2=rand(1);
if P2<0.0047;
    P3=rand(1);

     thetaruy=2.5*sqrt(1/P3)*sqrt(2)*theta0;  
     P4=rand(1);    
    if P4>0.5
    thetaruy=thetaruy;
    else
    thetaruy=-thetaruy;
    end
[rotang,pypar,pzpar] = rot(py,pz,thetaruy);

 
 py=pypar*cos(rotang)+pzpar*sin(rotang);
 pz=-pypar*sin(rotang)+pzpar*cos(rotang);

end


%==========================didn't hit on the collimator
else x=x+px/ptot*deltas;
     y=y+py/ptot*deltas;
     z=z+pz/ptot*deltas;
end
%============================For next step

xm(i)=x;
pxm(i)=px;
ym(i)=y;
pym(i)=py;
zm(i)=z;
pzm(i)=pz;
Em(i)=E;

end
plot(zm,ym,'r.','MarkerSize',1)
hold on
end