%=================================================================
%                         Define constant
%=================================================================
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

%------------------------------------------------------------------
%                      End of constant defination
%------------------------------------------------------------------
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
fid = fopen('uds_coeff.txt', 'w', 'n');
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
    s=0;
    %========================================================================
    for j=0:5000
        if (0.002*j*pi==0)
            continue;
        else
            t=fzero(uds_diff,0.002*j*pi);
            if (t>0 && t<14 && uds_diff2(t)>0)
                break;
                
            end
        end
    end
    t_tmp=[phi_s,t];
    fprintf(fid,'%12.8f %12.8f ', t_tmp);
    s=s+2;
    i_s=ceil(t/pi);
    i_e=ceil(100/pi);
    for i=i_s:2:i_e
        x=i*pi;
        for j=0:20
            if (x-0.1*j*pi==0)
                continue;
            else
                p=fzero(uds_diff,x-0.1*j*pi);
                p2=uds_diff2(p);
                if (p2>0 && p<=x)
                    
                    break;
                end
            end
        end
        for j=1:20
            if (x+0.1*j*pi==0)
                continue;
            else
                pp=fzero(uds_diff,x+0.1*j*pi);
                pp2=uds_diff2(pp);
                if (pp2>0 && pp>x)
                    break;
                end
            end
        end
        t_tmp=[p,pp];
        fprintf(fid,'%12.8f %12.8f ', t_tmp);
        s=s+2;
    end
    if s<34
        for ss=s+1:34
            fprintf(fid,'%12.8f ', NaN);
        end
    end
        
    fprintf(fid,'\n');
    sprintf('%s %d\n','Work on step ',k) 
end
fclose(fid);


            
            