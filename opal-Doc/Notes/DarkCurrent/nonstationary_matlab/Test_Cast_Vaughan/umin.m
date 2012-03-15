%=================================================================
%                         Define constant
%=================================================================
e=1.6E-019; % charge of electron (C)
m=9.10938188E-031;% mass of electorn (kg)
JtoeV=6.2415E+018;% (eV/J)
d=0.001;% distance between two plates (m)
f=1.64e9;% frequency of voltage (Hz)
w=2*pi*f;
Vo=120;% magnitude of voltage(V)
vw=e*Vo/m/w/d; % velocity scalar defination see S. Anza, C. Vicente, J. Gil, V. E. Boria, B. Gimeno, and D. Raboso, Phys. Plasmas 17, 062110 2010, eqn(3)
lambda=w*d/vw; % distance scalar defination see S. Anza, C. Vicente, J. Gil, V. E. Boria, B. Gimeno, and D. Raboso, Phys. Plasmas 17, 062110 2010, eqn(3)
vt=7.268929821E+5; % thermal velocity relative to energy of 1.5eV see S. Anza, C. Vicente, J. Gil, V. E. Boria, B. Gimeno, and D. Raboso, Phys. Plasmas 17, 062110 2010, eqn(12)
%------------------------------------------------------------------
%                      End of constant defination
%------------------------------------------------------------------
%=================================================================
%                         Define size
%=================================================================
%m=2;
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
                  
            [u,fval,exitflag]=fzero(g,-cos(phi_s));
            if (u+cos(phi_s))>=0
                
                ret=u;
                
                gg=@(t) (lambda-sin(phi_s)+sin(t+phi_s))/t-cos(phi_s)-u;
                if k<95
                    for i=10000:10000:100000000 % large interval for better converge to right number
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
                     
            [u,fval,exitflag]=fzero(g,-cos(phi_s)-cos(phi_s/2.0)); % use -cos(phi_s)-cos(phi_s/2.0) as initial guess for better converge see N. K. Vdovicheva, A. G. Sazontov, and V. E. Semenov, Radiophys. Quantum Electron. 47, 580 2004.p583-584
            if (u+cos(phi_s))>=0
                
                ret=u;
                
                gg=@(t) (lambda-sin(phi_s)+sin(t+phi_s))/t-cos(phi_s)-u;
                 for i=1:10000
                    [tao,fval,exitflag]=fzero(gg,1.0+i);
                    if exitflag==1
                        break;
                    end
                        
                    assert(u>=0)
                    assert(u<=1.2597) % According to N. K. Vdovicheva, A. G. Sazontov, and V. E. Semenov, Radiophys. Quantum Electron. 47, 580 2004. p583

                end
                                            
            end
            assert((u+cos(phi_s))>=0)
        end
    end
    
    z=[phi_s,ret,tao];
    fprintf(fid,'%12.8f %12.8f %12.8f\n', z);
    sprintf('%s %d\n','Work on step ',k) 
end
fclose(fid);


            
            