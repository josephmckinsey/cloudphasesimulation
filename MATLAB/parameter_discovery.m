% stechmann_simulation.m
% A simple script that defines a Stechmann2016 simulation and runs it.  The
% script does NOT do any data processing, plotting, etc., and it is
% recommended to do perform any such activities in a separate script.

% increase the number of threads for isengard
parpool(1);

% time for name 
time_name = datestr(datetime('now'));
save(time_name);

% Simulation parameters
xspan = [-100,100];
yspan = xspan;
end_time = 300; %length of simulation
tspan = [0,end_time];
Nx = 201;
Ny = Nx;
Nt = ceil(tspan(2)*(Nx+1)^2 / (1*(xspan(2)-xspan(1))^2));
sf = 0;    % Noise smoothing factor.  The greater, the smoother.  0 is completely uncorrelated.

n = 85; %number of samples along F and S
tsamples = 50:10:end_time;
f = linspace(-1,1,n); %source/sink coeff
s = linspace(0,10,n); %turbulence coeff

%resulting data matricies
%first index: index of source/sink param (f)
%second index: index of turbulence coeff (s)
%thrid index: timestep (tsamples)
sigmabar = NaN(n,n,length(tsamples)); 
meanq = NaN(n,n,length(tsamples));
varq = NaN(n,n,length(tsamples));

tic
for i = 1:n
    parfor j = 1:n
        % Problem parameters
        b = 25;    % Diffusion coefficient
        T = @(x,y,~) 100/24*ones(1,length(x),length(y));    % 
        F = @(x,y,~) f(i)*ones(1,length(x),length(y));    % 
        S = @(x,y,~) s(j)*ones(1,length(x),length(y));    % Turbulence coefficient
		
        % Generate source function and Jacobian
		[SRC,DSRC] = stechmann_functions_general(T,F,S);

		% Initial and boundary conditions
		% BCFUN = @(x,y,t,w) BCFUN_defn(x,y,t,w);
		ICFUN = @(x,y,w) ICFUN_defn(x,y,w);
		
		% Simulation
        [x,y,t,U] = SPDE_IE_2D_PBC(b,SRC,DSRC,xspan,Nx,yspan,Ny,tspan,Nt,ICFUN,sf);
        U = squeeze(U);

        %samplet = tsamples(k);
        sigmabar(i,j,:) = mean(mean(U(:,:,tsamples)>0));
        meanq(i,j,:) = mean(mean(U(:,:,tsamples)));
%         A = size(U);
%         varq(1,j,:) = var(reshape(U(:,:,tsamples),1,A(1)^2));
    end
    save(time_name,'sigmabar','meanq','varq','f','s','tsamples','-append');
end
toc

function f = ICFUN_defn(~,~,w)
    f = zeros(size(w));
end

% function f = BCFUN_defn(~,~,~,w)
%     f = zeros(size(w));
% end
