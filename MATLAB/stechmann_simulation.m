% stechmann_simulation.m
% A simple script that defines a Stechmann2016 simulation and runs it.  The
% script does NOT do any data processing, plotting, etc., and it is
% recommended to do perform any such activities in a separate script.

% Problem parameters
b = 25;    % Diffusion coefficient
T = @(x,y,~) 100/24*ones(1,length(x),length(y));    % 
F = @(x,y,~) -0.4*ones(1,length(x),length(y));    % 
S = @(x,y,~) 7.5*ones(1,length(x),length(y));    % Turbulence coefficient

% Generate source function and Jacobian
[SRC,DSRC] = stechmann_functions_general(T,F,S);

% Initial and boundary conditions
% BCFUN = @(x,y,t,w) BCFUN_defn(x,y,t,w);
ICFUN = @(x,y,w) ICFUN_defn(x,y,w);

% Simulation parameters
xspan = [-100,100];
yspan = xspan;
tspan = [0,100];
Nx = 201;
Ny = Nx;
Nt = ceil(tspan(2)*(Nx+1)^2 / (1*(xspan(2)-xspan(1))^2));
sf = 5;    % Noise smoothing factor.  The greater, the smoother.  0 is completely uncorrelated.

% Simulation
[x,y,t,U] = SPDE_IE_2D_PBC(b,SRC,DSRC,xspan,Nx,yspan,Ny,tspan,Nt,ICFUN,sf);
U = squeeze(U);

function f = ICFUN_defn(~,~,w)
    f = zeros(size(w));
end

% function f = BCFUN_defn(~,~,~,w)
%     f = zeros(size(w));
% end