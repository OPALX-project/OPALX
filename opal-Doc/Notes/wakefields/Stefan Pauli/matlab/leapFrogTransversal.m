function [x delatX] = leapFrogTransversal(t0, tend, dt, x0, dx0, a, energy)
% Implementation of the Leap Frog Algorithm
%
%Arguments:
% t0 = starting time
% tend = time of the solution : x(t=tend)
% dt = time steps
% x0 = starting point
% dx0 = starting velocity
% a = function for d^2x/d^2t (x,t)
%
% x = solution at timepoint tend
% trajektorie = [ x(t=t0) , x(t=t0+dt), ... x(t=tend)]

% Mit Steinwurf getestet:
% g=@(x,t,wake,profile) [ 0 -10 0]';
% [x, traj] =leapFrog(0, 3, .001, [0 0 0], [10 10 0],[0 0 0 ],g);
% plot(traj(1,:) , traj(2,:))


%dx1 = dx/dt(t-dt/2)
%dx2 = dx/dt(t+dt/2)
%dx = dx/dt(t)

% initial config
x = x0(:);
dx1 = dx0(:);
j=1;



% leapFrog Steps
for t = t0:dt:tend
    dx2 = dx1 + a(x-((t-t0)*dx0),t,energy*x(1))*dt;
    x = x +dx2*dt;
    dx1 = dx2;
  %  trajektorie(:,j) = x;
    j=j+1;
end

delatX = x(1)-((t-t0+dt)*dx0(1))-x0(1);
end