function  bfield = get_field_map(ns,radius_data, theta_data, radius, azimuth_dw, azimuth_up, bz_data, max_gap_width, add_flag,debug)

    r0 = radius_data(1);dr = radius_data(2);r1 = radius_data(3);nr = (r1 - r0) / dr + 1;
    t0 = theta_data(1); dt = theta_data(2); t1 = theta_data(3); nt = (t1 - t0) / dt;
    rad = [r0 : dr : r1]; th = [t0 : dt : t1 - dt];
    theta_up = pchip(radius,azimuth_up, rad);theta_dw = pchip(radius,azimuth_dw, rad);
    rmin = radius(1);rmax = radius(length(radius));
    
    if(debug(1) == 1)
        disp('enter to compare theta_up with theta_dw : ');
        plot(rad,theta_up,'r',rad,-theta_dw,'b');xlabel('r(cm)');ylabel('theta(deg)');legend('theta up','-theta dw');
        pause;
    end
    
    [bz,gap_width]=get_magnet_bz(bz_data,max_gap_width, rad);
    
    if(debug(2) == 1)
       disp('enter to plot Bz and gap width : '); 
       subplot(2,1,1);
        plot(bz_data(1,:),bz_data(2,:),'xb',rad,bz,'r');xlabel('r(cm)');ylabel('Bz(kGs)');
       subplot(2,1,2);
        plot(rad,gap_width,'r');xlabel('r(cm)');ylabel('gap width(cm)');
       pause;
    end
    
    swap = 360/ns;
    arc_to_rad = pi/180;
    
    for i = 1 : nr
        r = rad(i);
        if( r <= rmin || r >= rmax )
            for j = 1 : nt
                bfield(i,j) = 0.0;
            end
        else
            bfield(i,j) = 0.0;
            gap = gap_width(i);
            b = bz(i);
            for j = 1 : nt
                t = th(j);
                tmin = theta_dw(i);tmax = theta_up(i);tmid = (tmin + tmax)/2.0;
                
                if( t <= tmid )
                    
                    x = r*(tmin - t)*arc_to_rad/gap;fx = enge_empirical_fun(x);
                    if(add_flag == 1)
                       fx = fx + enge_empirical_fun(r*(t+swap-tmax)*arc_to_rad/gap); 
                    end
                    bfield(i,j) = b*fx;
                                       
                else
                    x = r*(t-tmax)*arc_to_rad/gap;fx = enge_empirical_fun(x);
                    if(add_flag == 1)
                       fx = fx + enge_empirical_fun(r*(tmin-t+swap)*arc_to_rad/gap); 
                    end
                    bfield(i,j) = b*fx;                    
                end                    
            end           
        end  
    end
    
    if(debug(3) == 1)
%add here
    end  
    
end