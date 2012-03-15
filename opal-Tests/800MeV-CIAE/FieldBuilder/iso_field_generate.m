clc;
rad_to_arc = 180/pi;
arc_to_rad = pi/180;
ns = 9;                     %Number of Sectors
ns_plot = 9;                %Number of Sectors to be plotted
sector_angle = 360/ns;sector_half_angle = sector_angle/2;
is_spiral = 1;             %is_siral = 0 : not taking spiral angle into account,i.e. add_angle = 0.0. is_spiral = 1 : taking apiral angle into account
bfield_flag = 1;           %bfield_flag = 0 : Field Map will not be calculated,i.e, just plot the accelerator structure.
                           %bfield_flag = 1 : Field Map will be calculated ang plotted.

if(is_spiral == 0)
    field_rotate_angle = 0;
else
    field_rotate_angle = 10.2;
    %Cavity Parameter
    cavity_length = 600;cavity_radius = 440;cavity_width = 60;cavity_azimuth = 31;cavity_rotate_angle = 13;
end

frf = 44.4;c = 299.792;har = 6;rinf = 100*har*c/(2*pi*frf);e0 = 938.272;
%intial_energy = [100 276 500 800];
intial_energy = [100:50:800];
intial_gama = (intial_energy+e0) / e0; intial_beta = sqrt(1 - 1.0./(intial_gama.*intial_gama));
intial_radius = rinf * intial_beta;
%magnet_angle =[13.02	14.44	16.12	17.76];
%magnet_angle at intial_energy = 100:50:800MeV;
magnet_angle = [12.92	13.374	13.826	14.243	14.622	15.032	15.417	15.785	16.12	16.434	16.736	17.02	17.286	17.532	17.76];

if(is_spiral == 0)
    add_angle = zeros(1,length(magnet_angle));
else
    %add_angle because of the spiral at intial_energy = 100:50:800MeV;
    add_angle  =  [0.000    2.0505	3.6702	5.1015	6.423	7.6596	8.8267	9.9317	10.978	11.967	12.902	13.783	14.61	15.371	16.06];
end

if( bfield_flag == 1)
    %radius_data : radial range,i.e.,230:1:580 ;  theta_data : azimuthal range,i.e., -20 : 0.2
    %bz_data: Max Field Bs used for interpolation. First row is radius,second is Field.
    radius_data = [230,1,580];theta_data = [-sector_half_angle,0.2,sector_half_angle];
    bz_data = [276.09	326.66	365.05	395.64	420.76	441.83	459.76	475.22	488.68	500.5	510.94	520.24	528.56	536.03	540 542.78
               15.520   15.872	16.159	16.452	16.750	16.999	17.267  17.525	17.790	18.071	18.372	18.673	18.972	19.290	19.4696 19.600-0.035]+0.002;
    max_gap_width = 10;add_flag = 1;debug_flag = [0,1,0];fieldmap_plot_type = 1;
end
intial_azimuth_up = magnet_angle/2 + add_angle; intial_azimuth_dw = -magnet_angle/2 + add_angle;

delta_energy = 10; energy = [100 : delta_energy : 800]; add_radius = 30;num = length(energy);
gama = (energy + e0) / e0; beta = sqrt(1 - 1.0./(gama.*gama));
radius = beta*rinf;radius_left = radius(1) - add_radius; radius_right = radius(num) + add_radius;
radius =[radius_left, radius, radius_right];

%azimuth_up = pchip(intial_radius,intial_azimuth_up,radius);
%azimuth_dw = pchip(intial_radius,intial_azimuth_dw,radius);
azimuth_up = spline(intial_radius,intial_azimuth_up,radius);
azimuth_dw = spline(intial_radius,intial_azimuth_dw,radius);
num_divide_up = 20;num_divide_dw = 20;
num = length(azimuth_dw);
azimuth_left = linspace(azimuth_dw(1),azimuth_up(1),num_divide_dw + 1);
azimuth_right = linspace(azimuth_dw(num),azimuth_up(num),num_divide_up + 1);

if( bfield_flag == 1)
    bfield = get_field_map(ns,radius_data, theta_data, radius, azimuth_dw, azimuth_up, bz_data, max_gap_width, add_flag,debug_flag);
    if(field_rotate_angle > 0)
        swap = floor(field_rotate_angle/theta_data(2));
        num = length(bfield(1,:));
        bfield_swap = bfield(:,1:swap);
        bfield = [bfield(:,swap+1:num),bfield_swap];
    elseif(field_rotate_angle < 0)
        swap = floor(field_rotate_angle/theta_data(2));
        num = length(bfield(1,:));
        bfield_swap = bfield(:,num-swap+1:num);
        bfield = [bfield_swap,bfield(:,1:num-swap)];        
    end
    
    if(fieldmap_plot_type ~= 0)
        disp('enter to plot field map picture : '); 
        plot_fieldmap(bfield,[radius_data(1):radius_data(2):radius_data(3)],[theta_data(1):theta_data(2):theta_data(3)-theta_data(2)],fieldmap_plot_type); 
        pause 
    end
    filename = 'cyclop.dat';format='%11.5f';format_num = 7;
    title = [0];
    savefield(filename, format, format_num,bfield,title,debug_flag);
    filename = 'opal_sectorfield.dat';format='%17.8e';format_num = 5;
    title = [2300,10,0,-5,200,351];
    savefield(filename, format, format_num,bfield,title,debug_flag);
end

plot(0,0,'o');axis equal;hold on;

%plot magnet sectors and cavity
if(is_spiral == 1)
    cavity_half_length = cavity_length / 2;cavity_half_width = cavity_width / 2;cavity_azimuth = cavity_azimuth + field_rotate_angle;
    cavity_angle = cavity_azimuth + cavity_rotate_angle;
    cavity_rad = cavity_radius/cos(cavity_rotate_angle*arc_to_rad);PDIS = cavity_radius*tan(cavity_rotate_angle*arc_to_rad);
    cavity_local_xy = [ - cavity_half_length, -cavity_half_width;
                        cavity_half_length, -cavity_half_width;
                        cavity_half_length,  cavity_half_width;
                        - cavity_half_length,  cavity_half_width;
                        - cavity_half_length, -cavity_half_width];
end
    
for i = 0 : ns_plot - 1
    ang = (sector_angle*i+field_rotate_angle)*arc_to_rad;
    x_left = radius_left*cos(azimuth_left*arc_to_rad + ang);y_left = radius_left*sin(azimuth_left*arc_to_rad + ang);
    x_right = radius_right*cos(azimuth_right*arc_to_rad + ang);y_right = radius_right*sin(azimuth_right*arc_to_rad + ang);
    x_dw = radius.*cos(azimuth_dw*arc_to_rad + ang);y_dw = radius.*sin(azimuth_dw*arc_to_rad + ang);
    x_up = radius.*cos(azimuth_up*arc_to_rad + ang);y_up = radius.*sin(azimuth_up*arc_to_rad + ang);
    
    plot(x_left,y_left,'b',x_right,y_right,'b',x_dw,y_dw,'b',x_up,y_up,'b');
    
    if((is_spiral == 1) && (mod(i,2) == 0))
        [cavity_azimuth,cavity_angle,PDIS];
        for k = 1 : 5
            cavity_global_xy(k,1) = cavity_rad*cos(cavity_azimuth*arc_to_rad) + cavity_local_xy(k,1)*cos(cavity_angle*arc_to_rad)-cavity_local_xy(k,2)*sin(cavity_angle*arc_to_rad);
            cavity_global_xy(k,2) = cavity_rad*sin(cavity_azimuth*arc_to_rad) + cavity_local_xy(k,1)*sin(cavity_angle*arc_to_rad)+cavity_local_xy(k,2)*cos(cavity_angle*arc_to_rad);
        end
        plot(cavity_global_xy(:,1),cavity_global_xy(:,2));
        cavity_azimuth = cavity_azimuth + 2*sector_angle;
        cavity_angle = cavity_angle + 2*sector_angle;
    end
end

is_eo_plot = 0; %if you want to plot the seo, please give the eo_data.mat
if(is_spiral == 1 && is_eo_plot == 1)
    load('eo_data.mat');
    number_of_energy = 4;number_of_data = 100;
    for i = 1:number_of_energy
        k1 = number_of_data*(i-1)+1;k2 = k1 + number_of_data - 1;
        eo_rad = eo_data(k1 : k2,2);eo_th = eo_data(k1 : k2,1);
        eo_x = [];eo_y = [];
        for j = 1 : ns_plot
            eo_x = [eo_x,eo_rad.*cos((eo_th+(j-1)*sector_angle)*arc_to_rad)];
            eo_y = [eo_y,eo_rad.*sin((eo_th+(j-1)*sector_angle)*arc_to_rad)];
        end
        if( i ~=3 )
            plot(eo_x,eo_y,'--k');
        end
    end
end

%plot line, so you can check whether the field_rotate_angle is right. 
%if you do not want to do that, just annotate it with a %
rmax = 580;
plot([0,rmax],[0,0],'--k',[0,rmax*cos(sector_angle*arc_to_rad)],[0,rmax*sin(sector_angle*arc_to_rad)],'--k');
xlabel('X[cm]');ylabel('Y[cm]');



