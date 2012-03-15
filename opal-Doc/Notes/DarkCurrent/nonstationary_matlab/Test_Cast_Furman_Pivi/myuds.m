function Y = myuds(x,lambda,phi_s)
global uds_coeff;
global u_min;

k_e=size(uds_coeff,1);
i_e=size(uds_coeff,2);
for k=1:k_e
    if (abs(uds_coeff(k,1)-phi_s)<0.00001)
        if u_min(k,2)>0
            if x>=u_min(k,3)
                Y=(lambda-sin(phi_s)+sin(u_min(k,3)+phi_s))/u_min(k,3)-cos(phi_s);
                assert(abs(Y-u_min(k,2))<0.00001)
                return
            elseif (x<=uds_coeff(k,2))
                if ( x~=0 )
                    Y=(lambda-sin(phi_s)+sin(x+phi_s))/x-cos(phi_s);
                    if Y>=0
                        return
                    else
                        Y=0;
                        return
                    end
                else
                    Y=(lambda-sin(phi_s)+sin(0.00001+phi_s))/0.00001-cos(phi_s);%Fix me:to calculate SEY near x=0;
                    if Y>=0
                        Y=0;
                        return
                    else
                        Y=0;
                        return
                    end
                end
            else
                for i=3:2:i_e-1
                    if (isnan(uds_coeff(k,i)))
                        break;
                    else
                        if (x>=uds_coeff(k,i) && x<= uds_coeff(k,i+1))
                            if ((lambda-sin(phi_s)+sin(x+phi_s))/x-cos(phi_s)<=(lambda-sin(phi_s)+sin(uds_coeff(k,i)+phi_s))/uds_coeff(k,i)-cos(phi_s))
                                ret=(lambda-sin(phi_s)+sin(x+phi_s))/x-cos(phi_s);
                            else
                                ret=(lambda-sin(phi_s)+sin(uds_coeff(k,i)+phi_s))/uds_coeff(k,i)-cos(phi_s);
                                
                            end
                            if ret>=0
                                Y=ret;
                                return
                            else
                                Y=0;
                                return
                            end
                        end
                    end
                end
            end
        else
            if (x<=uds_coeff(k,2))
                if ( x~=0 )
                    Y=(lambda-sin(phi_s)+sin(x+phi_s))/x-cos(phi_s);
                    if Y>=0
                        return
                    else
                        Y=0;
                        return
                    end
                else
                    Y=(lambda-sin(phi_s)+sin(0.00001+phi_s))/0.00001-cos(phi_s);%Fix me:to calculate SEY near x=0;
                    if Y>=0
                        return
                    else
                        Y=0;
                        return
                    end
                end
            else
                for i=3:2:i_e-1
                    if (isnan(uds_coeff(k,i)))
                        break;
                    else
                        if (x>=uds_coeff(k,i) && x<= uds_coeff(k,i+1))
                            if ((lambda-sin(phi_s)+sin(x+phi_s))/x-cos(phi_s)<=(lambda-sin(phi_s)+sin(uds_coeff(k,i)+phi_s))/uds_coeff(k,i)-cos(phi_s))
                                ret=(lambda-sin(phi_s)+sin(x+phi_s))/x-cos(phi_s);
                            else
                                ret=(lambda-sin(phi_s)+sin(uds_coeff(k,i)+phi_s))/uds_coeff(k,i)-cos(phi_s);
                                
                            end
                            if ret>=0
                                Y=ret;
                                return
                            else
                                Y=0;
                                return
                            end
                        end
                    end
                end
            end
        end
    end
end
Y=[];
a=1;
assert(a==0,'x=%f ,phi_s=%f,',x,phi_s);    

            
            