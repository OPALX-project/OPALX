function plot_fieldmap(bfield,rad,th,selection,plot_sector_number)
        if nargin == 4
            plot_sector_number = 9;
        end
        num = length(th);dth = th(2) - th(1);
        sector_angle = num*dth;
        field = bfield;the = th;
        for i = 2 : plot_sector_number
            field = [field,bfield];
            the = [the,th + (i - 1) * sector_angle];
        end

        [th,rad] = meshgrid(the*pi/180,rad);
        [x,y,z] = pol2cart(th,rad,field);
        if(selection == 1)
            mesh(x,y,z);title('mesh(selection = 1)');
            xlabel('X[cm]');ylabel('Y[cm]');zlabel('Bz[kGs]');
        elseif(selection == 2)
            meshc(x,y,z);title('meshc(selection = 2)');
            xlabel('X[cm]');ylabel('Y[cm]');zlabel('Bz[kGs]');
        elseif(selection == 3)
            meshz(x,y,z);title('meshz(selection = 3)');
            xlabel('X[cm]');ylabel('Y[cm]');zlabel('Bz[kGs]');
        elseif(selection == 4)
            contour(x,y,z,50);title('Equipotential magnetical field lines of one sector');
            xlabel('X[cm]');ylabel('Y[cm]');
        else
           %add 
        end
end