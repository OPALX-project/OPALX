%function Y=myvoltra(n,t0,h,vw,vt,lambda)
clear;
e=1.6E-019;
m=9.10938188E-031;
JtoeV=6.2415E+018;%eV/J
d=0.001;%m
f=1.64e9;%Hz
w=2*pi*f;
Vo=120;%V
vw=e*Vo/m/w/d;
lambda=w*d/vw;
%fid = fopen('Gds_Gss.txt', 'w', 'n');
phi_s=122/180*pi;
phi=mod(phi_s,2*pi);
vt=7.268929821E+5;
m=2;
h=pi/180;
% h=w*5*1E-12;
n=ceil(5*1E-9*w/h);% total time step equals to time period you want, divided by the time interval. here 5ns/dt
t0=0;%start phase;
myt=zeros(n,1);
%myphi_s=zeros(n,1);
myf=zeros(m,n);
myt(1,1)=t0;
A=zeros(m);
B=zeros(m,1);
BB=zeros(m,1);
myI=zeros(m,n);
global uss_coeff;
uss_coeff=load( 'myuss_coeff.txt');
global uds_coeff;
uds_coeff=load( 'myuds_coeff.txt');
global u_min;
u_min=load('u_min.txt');
for i=1:n
    myt(i,1)=myt(1,1)+(i-1)*h;
    %myphi_s(i,1)=w*myt(i,1);
end
for k=1:m
    myf(k,1)=seedg(k,myt(1,1));
end
for i=2:n
    %myt(i,1)=myt(i-1,1)+h;
    for k=1:m
        sum=seedg(k,myt(i,1));
        for l=1:m
            sum=sum+0.25*ak(k,l,myt(i,1),myt(1,1),0.01*h,vw,vt,lambda);
            for j=2:(i-1)
                sum=sum+h*ak(k,l,myt(i,1),myt(j,1),0.01*h,vw,vt,lambda)*myf(l,j);
            end
            if (k==l)
                tmp=ak(k,l,myt(i,1),myt(i,1),0.01*h,vw,vt,lambda);
                A(k,l)=1.0-0.5*h*tmp;
            else
                temp=ak(k,l,myt(i,1),myt(i,1),0.01*h,vw,vt,lambda);
                A(k,l)=-0.5*h*temp;
            end
        end
        B(k,1)=sum;
    end
    BB=A\B;
    for k=1:m
        myf(k,i)=BB(k,1);
    end
    sprintf('%s %d\n','Work on step ',i) 
end

for k=1:m
    myI(k,1)=0;
end
for i=2:n
    for k=1:m%fix me
        sum=0.25*myG(k,1,myt(i,1),myt(1,1),0.01*h,vw,vt,lambda)+0.25*myG(k,2,myt(i,1),myt(1,1),0.01*h,vw,vt,lambda)+0.5*h*(myf(1,i)*myG(k,1,myt(i,1),myt(i,1),0.01*h,vw,vt,lambda)+myf(2,i)*myG(k,2,myt(i,1),myt(i,1),0.01*h,vw,vt,lambda));
        for j=2:i-1 %fix me
            sum=sum+0.5*h*(2*myf(1,j)*myG(k,1,myt(i,1),myt(j,1),0.01*h,vw,vt,lambda)+2*myf(2,j)*myG(k,2,myt(i,1),myt(j,1),0.01*h,vw,vt,lambda));
        end
        myI(k,i)=sum;
    end          
end
for i=1:n
    myf(1,i)=myf(1,i)*w*1e-9;
    myf(2,i)=myf(2,i)*w*1e-9;
    myI(1,i)=myI(1,i)*w*1e-9;
    myI(2,i)=myI(2,i)*w*1e-9;
end
dt=h/w/1e-9;
myNp=zeros(n,1);
myNp(1,1)=1.0; %first step particle number equals to unit;
for i=2:n
    sumNp=0.5*(1-dt*(myI(1,1)+myI(2,1)));
    for j=2:i-1
       sumNp=sumNp+dt*(myf(1,j)+myf(2,j)-myI(1,j)-myI(2,j));
    end
    sumNp=sumNp+dt*0.5*(myf(1,i)+myf(2,i)-myI(1,i)-myI(2,i));
    myNp(i,1)=sumNp/myNp(1,1);
end
fid = fopen('fig8_result.txt', 'w', 'n');
for i=2:n
    tmp_z=[(i-1)*h/w*1E9,myNp(i,1)];
%     tmp_z=[(i-1)*h/w*1E9,(myf(1,i)+myf(2,i))/(myI(1,i)+myI(2,i))];
    fprintf(fid,'%12.8f %12.8f\n', tmp_z);
end
fclose(fid);

