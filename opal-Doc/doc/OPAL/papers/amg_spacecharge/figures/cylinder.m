t = 0:pi/20:2*pi;
[X,Y,Z] = cylinder(2+.5*sin(t).*cos(.5*t));
surf(Z,X,Y)
axis off
shading interp
colormap(pink)
