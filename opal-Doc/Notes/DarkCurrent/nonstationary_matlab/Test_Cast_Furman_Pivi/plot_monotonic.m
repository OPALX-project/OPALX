e=1.6E-019;
m=9.10938188E-031;
JtoeV=6.2415E+018;%eV/J
d=0.005;%m
f=2.0e8;%Hz
w=2*pi*f;
Vo=120;%V
vw=e*Vo/m/w/d;
lambda=w*d/vw;
% fid = fopen('Gds_Gss.txt', 'w', 'n');
phi_s=122/180*pi;
phi=mod(phi_s,2*pi);
vt=7.268929821E+5;
global uss_coeff;
uss_coeff=load( 'myuss_coeff.txt');
global uds_coeff;
uds_coeff=load( 'myuds_coeff.txt');
global u_min;
u_min=load( 'u_min.txt');
h=pi/180;
n=ceil(2*pi/h);
t0=0;%start phase;

myphi_s=zeros(n,1);
myphi_s(1,1)=t0;
for i=1:n
    
    myphi_s(i,1)=myphi_s(1,1)+(i-1)*h;
end
for k=1:20:n
    phi_s=myphi_s(k,1);
    xx=zeros(950,1);
    yy=zeros(950,1);
    zz=zeros(950,1);
    for i=1:950,
        xx(i,1)=0.1*i;
        yy(i,1)=myuss(xx(i,1),lambda,phi_s);
%         yy(i,1)=mygds(xx(i,1),0.01*h,vw,vt,lambda,phi_s);
%         zz(i,1)=mygss(xx(i,1),0.01*h,vw,vt,lambda,phi_s);
%         yy(i,1)=(xx(i,1)+cos(phi_s))*(2*pi+acos(xx(i,1)+cos(phi_s))-phi_s)+sin(phi_s)-sin(acos(xx(i,1)+cos(phi_s)));
%         zz(i,1)=xx(i,1)+cos(phi_s);
        zz(i,1)=myuds(xx(i,1),lambda,phi_s);
%         zz(i,1)=(0-sin(phi_s)+sin(xx(i,1)+phi_s))/xx(i,1)-cos(phi_s);
%         yy(i,1)=(lambda-sin(phi_s)+sin(xx(i,1)+phi_s))/xx(i,1)-cos(phi_s);
    end
     scatter(xx,yy);
     hold on;
     scatter(xx,zz);
     hold on;
%     plot(xx,yy);
%     hold on;
%     plot(xx,zz);
%     hold on;
% 
     pause
end


                

