function [x,y,t,U] = SPDE_IE_2D_PBC(Dvec, source, dsource, xspan, Nx, yspan, Ny, tspan, Nt, ICfun, sf)
%SPDE_IE_2D  Simulate systems of stochastic PDEs with periodic
%  boundary conditions in two spatial dimensions, implicit method.
%  [X,Y,T,U] = IE_SPDE_1D(DVEC,SRC,DSRC,XSPAN,NX,YSPAN,NY,TSPAN,NT,ICFUN,SF)
%  simulates stochastic PDEs of the form
%  du(i)/dt = DVEC(i)*(d2u(i)/dx2 + d2u(i)/dy2) + SRC(x,y,t,u,W)
%  where W is random white noise.  SPDE_IE_2D uses a second-order finite
%  difference scheme in space and the implicit Euler method for time
%  integration.
%  
%  ---INPUTS---
%  DVEC is a vector of diffusion coefficients with length(DVEC) = Neq =
%  number of equations in the system. SRC is a function handle. For vectors
%  x, y, u, and W and scalar t, SRC(x,y,t,u,W) must return an array of
%  dimensions Neq by length(x) by length(y). DSRC is a function handle,
%  analagous to the Jacobian of the source function.  DSRC(x,y,t,u,W) must
%  return an array of dimensions Neq by Neq by length(x) by length(y),
%  where the i,j,n,m-th entry is the i,j-th entry of the Jacobian of
%  source, evaluated at (x(n),y(n)). XSPAN(1) is the left boundary;
%  xspan(2) is the right boundary. YSPAN is analagous. TSPAN(1) is the
%  initial time; solutions will be simulated up to TSPAN(2). NX, NY, and NT
%  are the number of spatial points in x and y and number of times
%  at which the solution will be evaluated, respectively.  ICFUN is a
%  function handle.  For vectors x, y, and W, ICfun(x,y,W) must return an
%  array of dimensions Neq by length(x) by length(y) that represents the
%  initial condition for the system.  SF is a parameter that determines how
%  smooth the random noise will be; if  SF==0, then the noise will be
%  uncorrelated.  If SF>0, then noise will be diffused by SF
%  timesteps before being used.
%
%  ---OUTPUTS---
%  X = linspace(XSPAN(1), XSPAN(2), NX) are the (uniformly spaced)
%  spatial evaluation points in x, including the boundary. Y is analagous.
%  T = linspace(TSPAN(1), TSPAN(2), NT+1) are the (uniformly spaced)
%  temporal evaluation points, including the initial condition. U is an Neq
%  by NX by NY by NT+1 array such that U(i,n,m,j) approximates a solution to
%  the ith equation at X(n), Y(m), and T(j).

%  Author: Cooper Brown
    
    %% Initialization
    Dvec = Dvec(:);
    Neq = length(Dvec);
    x = linspace(xspan(1), xspan(2), Nx);
    y = linspace(yspan(1), yspan(2), Ny);
    t = linspace(tspan(1), tspan(2), Nt+1);
    
    % Discretization parameters
    dx = (xspan(2)-xspan(1))/(Nx-1);
    dy = (yspan(2)-yspan(1))/(Ny-1);
    dt = (tspan(2)-tspan(1))/Nt;
    Dx = dt/dx^2*Dvec;
    Dy = dt/dy^2*Dvec;
    dsto = (dt*dx*dy)^-0.5;
    
    % Generate (diffused) noise
    Wvec = NaN(Neq,Nx,Ny,Nt+1);
    for i=1:Neq
        Wvec(i,:,:,:) = dsto*NOISE_SMOOTH([Nx,Ny,Nt+1], [dx,dy,dt], sf);
    end
    
    % Initialize solution array, evaluate initial/boundary conditions
    U = NaN(Neq, Nx, Ny, Nt+1);
    U(:,:,:,1) = ICfun(x, y, Wvec(:,:,:,1));
    Uvec = permute(U(:,:,:,1), [2,3,1]);
    Uvec = reshape(Uvec, [Neq*Nx*Ny, 1]);    % For Newton iteration
    
    % Initialize block-diagonal iteration matrix
    Amat = kron(spdiags(Dx,0,Neq,Neq), kron(speye(Ny), spdiags([-1,-1,2,-1,-1].*ones(Nx,1),[-Nx+1,-1:1,Nx-1],Nx,Nx))) ...
        + kron(spdiags(Dy,0,Neq,Neq), kron(spdiags([-1,-1,2,-1,-1].*ones(Ny,1),[-Ny+1,-1:1,Ny-1],Ny,Ny), speye(Nx))) ...
        + speye(Neq*Nx*Ny);
    
    % Order to fill in Jacobian-ish matrix
    k = kron(ones(Neq^2,1), (1:Nx*Ny)') + kron((0:Neq-1)', Nx*ones(Neq*Nx*Ny, 1));
    l = kron(ones(Neq,1), (1:Nx*Ny*Neq)');
    
    %% Time integration using the implicit Euler method
    for i=2:Nt+1
        
        Upre = Uvec;    % Uvec currently holds the previous values
        
        % Newton iteration
        for j=1:1    % 3 seems to work well enough. 1 if the system is linear
            
            % Source vector fvec, with boundary conditions
            Ueval = permute(reshape(Uvec, [Nx, Ny, Neq]), [3,1,2]);
            fvec = dt*source(x, y, t(i), Ueval, Wvec(:,:,:,i));
            fvec = reshape(permute(fvec, [2,3,1]), [Neq*Nx*Ny, 1]);
            
            % Jacobian matrix term dfmat
            df = dt*dsource(x, y, t(i), Ueval, Wvec(:,:,:,i));
            df = reshape(permute(df, [3,4,1,2]), [Neq^2*Nx*Ny, 1]);
            dfmat = sparse(k, l, df);
            
            % Finalize the linear system and solve
            dfmat = Amat - dfmat;
            fvec = fvec + Upre - Amat*Uvec;
            Uvec = Uvec + dfmat\fvec;
            
        end
        
        % Add to solution array
        U(:,:,:,i) = permute(reshape(Uvec, [Nx, Ny, Neq]), [3,1,2]);
        
    end
    
end