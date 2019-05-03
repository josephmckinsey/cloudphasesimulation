function [x, y, t, U] = EE_SPDE_2D(Dvec, source, xspan, Nx, yspan, Ny, tspan, Nt, ICfun)
%EE_SPDE_2D  Simulate systems of stochastic PDEs with zero Dirichlet
%  boundary conditions in two spatial dimensions, explicit method.
%  [x,y,t,U] = EE_SPDE_1D(Dvec, source, xspan, Nx, tspan, Nt, ICfun)
%  simulates stochastic PDEs of the form
%  du(i)/dt = Dvec(i)*(d2u(i)/dx2+d2u(i)/dy2) + source(x,y,t,u,W)
%  on a rectangular domain, where W is random white noise.  The solution is
%  assumed to be 0 on the boundary.  EE_SPDE_2D uses a second-order finite
%  difference scheme in space and the explicit Euler method for time
%  integration.
%  
%  INPUTS 
%  Dvec is a vector of diffusion coefficients with length(Dvec) = Neq = 
%  number of equations in the system.
%  source is a function handle. For vectors x, y, u, Neq by length(x)
%  matrix W, and scalar t, source(x,y,t,u,W) must return an array of
%  dimensions Neq by length(x) by length(y).
%  xspan(1) is the left boundary in x; xspan(2) is the right boundary.
%  yspan works the same way, for y.
%  tspan(1) is the initial time; solutions will be simulated up to
%  tspan(2).
%  Nx, Ny, and Nt are the number of interior spatial points in x and y and
%  number of times at which the solution will be evaluated, respectively.
%  ICfun is a function handle.  For a vector x, ICfun(x,y) must return an
%  array of dimensions Neq by length(x) by length(y) that represents the
%  initial condition for the system.
%
%  OUTPUTS
%  x = linspace(xspan(1), xspan(2), Nx+2) are the (uniformly spaced)
%  spatial evaluation points in x.
%  y = linspace(yspan(1), yspan(2), Ny+2) are the (uniformly spaced)
%  spatial evaluation points in y.
%  t = linspace(tspan(1), tspan(2), Nt+1) are the (uniformly spaced)
%  temporal evaluation points.
%  U is an Neq by Nx+2 by Ny+2 by Nt+1 array such that U(i,j,k,l)
%  approximates a solution to the ith equation at x(j), y(k), and t(l).

%  Author: Cooper Brown

    Neq = length(Dvec);
    x = linspace(xspan(1), xspan(2), Nx+2);
    y = linspace(yspan(1), yspan(2), Ny+2);
    t = linspace(tspan(1), tspan(2), Nt+1);
    
    U = NaN(Neq, Nx+2, Ny+2, Nt+1);
    U(:,:,:,1) = ICfun(x,y);
    U(:,1,:,:) = 0;    % Zero Dirichlet boundary conditions
    U(:,end,:,:) = 0;
    U(:,:,1,:) = 0;
    U(:,:,end,:) = 0;
    
    dx = (xspan(2)-xspan(1))/(Nx+1);
    dy = (yspan(2)-yspan(1))/(Ny+1);
    dt = (tspan(2)-tspan(1))/Nt;
    Dx = dt/dx^2*Dvec;
    Dy = dt/dy^2*Dvec;
    dsto = (dt*dx*dy)^-0.5;
    
    % Time integration using the explicit Euler method
    for i=2:Nt+1
        U(:,2:end-1,2:end-1,i) = U(:,2:end-1,2:end-1,i-1) ...
            + Dx.*(U(:,1:end-2,2:end-1,i-1) - 2*U(:,2:end-1,2:end-1,i-1) + U(:,3:end,2:end-1,i-1)) ...
            + Dy.*(U(:,2:end-1,1:end-2,i-1) - 2*U(:,2:end-1,2:end-1,i-1) + U(:,2:end-1,3:end,i-1)) ...
            + dt*source(x(2:end-1), y(2:end-1), t(i-1), U(:,2:end-1,2:end-1,i-1), dsto*randn(Neq,Nx,Ny));
    end
    
end