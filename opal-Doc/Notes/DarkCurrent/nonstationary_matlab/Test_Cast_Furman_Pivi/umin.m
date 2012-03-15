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
fid = fopen('u_min.txt', 'w', 'n');
for k=1:n
    if k<269
        phi_s=myphi_s(k,1);
        if phi_s>=0 && phi_s<=pi/2.0
            ret=0;
            tao=0;
        else
            g=@(u) (u+cos(phi_s))*(2*pi+acos(u+cos(phi_s))-phi_s)+sin(phi_s)-sin(acos(u+cos(phi_s)));
            
            %         sym u;
            %         Y_in =(2*pi+acos(u+cos(phi_s))-phi_s)+sin(phi_s)-sin(acos(u+cos(phi_s)));
            %         YY = diff(Y_in);
            %         fc = char(YY);
            %         fcv = symvar(fc); % Note: these may not be in the order you want.if the expression involves multiple variables
            %         f2 = ['@(', sprintf('%s,', fcv{1:end-1}), fcv{end}, ') ', fc];
            %         ineq_diff=eval(f2);
            
            [u,fval,exitflag]=fzero(g,-cos(phi_s));
            if (u+cos(phi_s))>=0
                
                ret=u;
                
                gg=@(t) (lambda-sin(phi_s)+sin(t+phi_s))/t-cos(phi_s)-u;
                if k<95
                    for i=10000:10000:100000000
                        [tao,fval,exitflag]=fzero(gg,1.0+i);
                        if exitflag==1
                            break;
                        end
                        
                        assert(u>=0)
                        assert(u<=1.2597)
                    end
                else
                    for i=1:10000
                        [tao,fval,exitflag]=fzero(gg,1.0+i);
                        if exitflag==1
                            break;
                        end
                        
                        assert(u>=0)
                        assert(u<=1.2597)
                    end
                end
                %
                %         else
                %             ret=-cos(phi_s);
                
            end
            assert((u+cos(phi_s))>=0)
        end
    else
        phi_s=myphi_s(k,1);
        if phi_s>=0 && phi_s<=pi/2.0
            ret=0;
            tao=0;
        else
            g=@(u) (u+cos(phi_s))*(2*pi+acos(u+cos(phi_s))-phi_s)+sin(phi_s)-sin(acos(u+cos(phi_s)));
            
            %         sym u;
            %         Y_in =(2*pi+acos(u+cos(phi_s))-phi_s)+sin(phi_s)-sin(acos(u+cos(phi_s)));
            %         YY = diff(Y_in);
            %         fc = char(YY);
            %         fcv = symvar(fc); % Note: these may not be in the order you want.if the expression involves multiple variables
            %         f2 = ['@(', sprintf('%s,', fcv{1:end-1}), fcv{end}, ') ', fc];
            %         ineq_diff=eval(f2);
            
            [u,fval,exitflag]=fzero(g,-cos(phi_s)-cos(phi_s/2.0));
            if (u+cos(phi_s))>=0
                
                ret=u;
                
                gg=@(t) (lambda-sin(phi_s)+sin(t+phi_s))/t-cos(phi_s)-u;
                 for i=1:10000
                    [tao,fval,exitflag]=fzero(gg,1.0+i);
                    if exitflag==1
                        break;
                    end
                        
                    assert(u>=0)
                    assert(u<=1.2597)
                end
                               
                %
                %         else
                %             ret=-cos(phi_s);
                
            end
            assert((u+cos(phi_s))>=0)
        end
    end
    
    z=[phi_s,ret,tao];
    fprintf(fid,'%12.8f %12.8f %12.8f\n', z);
            
            
            
    sprintf('%s %d\n','Work on step ',k) 
end
fclose(fid);


            
            