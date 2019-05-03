Dvec = 25;
xspan = [-10,10];
yspan = xspan;
tspan = [0,10];

Nx = 100;
Ny = Nx;
Nt = ceil(tspan(2)*(Nx+1)^2 / (0.01*(xspan(2)-xspan(1))^2));
[x,y,t,U] = EE_SPDE_2D(Dvec,@source,xspan,Nx,yspan,Ny,tspan,Nt,@ICfun);

for i=1:10:size(U,4)
    figure(1);
    h = pcolor(x,y,squeeze(U(1,:,:,i))');
    set(h, 'EdgeColor', 'none');
%     axis([-1, 1, -1, 1, -2, 2]);
end

function f = source(~, ~, ~, u, w)
    f = -1/100*u - 0.4 + 7.5*w;
end

function f = ICfun(x,y)
    x = x(:);
    y = y(:);
    f = sin(pi()*x)*sin(pi()*y)';
%     f = zeros(size(x), size(y));
end