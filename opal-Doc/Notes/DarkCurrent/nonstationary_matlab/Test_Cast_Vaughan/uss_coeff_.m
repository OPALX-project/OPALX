e=1.6E-019;
m=9.10938188E-031;
JtoeV=6.2415E+018;%eV/J
d=0.001;%m
f=1.64e9;%Hz
w=2*pi*f;
Vo=120;%V
vw=e*Vo/m/w/d;
lambda=w*d/vw;
% fid = fopen('Gds_Gss.txt', 'w', 'n');
phi_s=122/180*pi;
phi=mod(phi_s,2*pi);
vt=7.268929821E+5;
global u_min;
u_min=load('u_min.txt');
%=================================================================
%                         Define size
%=================================================================
m=2;
h=pi/180;
n=ceil(2*pi/h);
t0=0;%start phase;

myphi_s=zeros(n,1);
myphi_s(1,1)=t0;
for i=1:n
    
    myphi_s(i,1)=myphi_s(1,1)+(i-1)*h;
end
%------------------------------------------------------------------
%                      End of size defination
%------------------------------------------------------------------
fid = fopen('uss_coeff.txt', 'w', 'n');
%========================================================================

%===============================================================================
lambda=lambda*0;
for k=1:n
    phi_s=myphi_s(k,1);
    syms t;
    Y_in = (lambda-sin(phi_s)+sin(t+phi_s))/t-cos(phi_s);
    YY = diff(Y_in);
    fc = char(YY);
    fcv = symvar(fc); % Note: these may not be in the order you want.if the expression involves multiple variables
    f2 = ['@(', sprintf('%s,', fcv{1:end-1}), fcv{end}, ') ', fc];
    uds_diff = eval(f2);
    YYY= diff(Y_in,t,2);
    fc2 = char(YYY);
    fcv2 = symvar(fc2); % Note: these may not be in the order you want.if the expression involves multiple variables
    f3 = ['@(', sprintf('%s,', fcv2{1:end-1}), fcv2{end}, ') ', fc2];
    uds_diff2 = eval(f3);
    g=@(x) (lambda-sin(phi_s)+sin(x+phi_s))/x-cos(phi_s)-u_min(k,2);
    %========================================================================
    for j=0:100
        if (0.1*j*pi==0)
            continue;
        else
            t=fzero(uds_diff,0.1*j*pi);
            
%             if (t~=0 && ((lambda-sin(phi_s)+sin(t+phi_s))/t-cos(phi_s)>=0))
%                 break;
            
            if (t>0 && uds_diff2(t)<0)
                break
              
            end
        end
    end
    if (((lambda-sin(phi_s)+sin(t+phi_s))/t-cos(phi_s)<0)||u_min(k,3)==0)
        x=0;
    elseif (abs((lambda-sin(phi_s)+sin(t+phi_s))/t-cos(phi_s)-u_min(k,2))<0.000001)
        x=t;
    else
        [x,fvar,exitflag]=fzero(g,t+pi/2.0);
        assert(x>0 && exitflag==1,'%f %f %f %f %f ',x,exitflag,t,u_min(k,2),(lambda-sin(phi_s)+sin(t+phi_s))/t-cos(phi_s)) 
        
    end            
    t_tmp=[phi_s,x];
    fprintf(fid,'%12.8f %12.8f\n', t_tmp);
    sprintf('%s %d\n','Work on step ',k) 
end
fclose(fid);

            
            