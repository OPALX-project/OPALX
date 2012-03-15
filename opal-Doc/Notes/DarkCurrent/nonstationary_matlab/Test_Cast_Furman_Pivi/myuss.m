function Y = myuss(x,lambda,phi_s)
% The uss function
% The following will be function body
lambda=lambda*0;
% syms t;
% Y_in = (lambda-sin(phi_s)+sin(t+phi_s))/t-cos(phi_s);
% YY = diff(Y_in);
% fc = char(YY);
% fcv = symvar(fc); % Note: these may not be in the order you want.if the expression involves multiple variables
% f2 = ['@(', sprintf('%s,', fcv{1:end-1}), fcv{end}, ') ', fc];
% uds_diff = eval(f2);
% YYY= diff(Y_in,t,2);
% fc2 = char(YYY);
% fcv2 = symvar(fc2); % Note: these may not be in the order you want.if the expression involves multiple variables
% f3 = ['@(', sprintf('%s,', fcv2{1:end-1}), fcv2{end}, ') ', fc2];
% uds_diff2 = eval(f3);
%========================================================================
% for j=0:20
%     if (0.1*j*pi==0)
%         continue;
%     else
%         t=fzero(uds_diff,0.1*j*pi);
%         if (t~=0 && ((lambda-sin(phi_s)+sin(t+phi_s))/t-cos(phi_s)>0))
%             break;
%        
%         end
%     end
% end
% load uss_coeff.txt;
global uss_coeff;
global u_min;
k_e=size(uss_coeff,1);
for k=1:k_e
    
    if (abs(uss_coeff(k,1)-phi_s)<0.0001)
        if u_min(k,2)>0
            t=uss_coeff(k,2);
            if t==0
                Y=0;
                return
            else
                if (x>t)
                    
                    %                 Y=(lambda-sin(phi_s)+sin(t+phi_s))/t-cos(phi_s);
%                     Y=u_min(k,2);
                    Y=0;
                    if Y>=0
                        return
                    else
                        Y=0;
                        return
                    end
                    
                else
                    if(((lambda-sin(phi_s)+sin(x+phi_s))/x-cos(phi_s)>=0) && x~=0 )
                        if ((lambda-sin(phi_s)+sin(x+phi_s))/x-cos(phi_s)<=u_min(k,2))
                            Y=(lambda-sin(phi_s)+sin(x+phi_s))/x-cos(phi_s);
                        else
                            Y=u_min(k,2);
                        end
                        if isinf(Y) || isnan(Y)
                            y3=Y
                            lam=lambda
                        end
                        return
                    else
                        Y=0;
                        %         Y=-1;
                        return
                    end
                end
            end
            
        else
            Y=0;
            return
        end
    end
end
Y=[];
a=1;
assert(a==0,'x=%f ,phi_s=%f,',x,phi_s);
   

            
            