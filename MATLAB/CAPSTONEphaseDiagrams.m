n = 100;

F = linspace(-1,1,n);
D = linspace(0,15,n);

deltaX = 1.0;
b = 25.0;
t = 100.0/24;
L = 550.0;


sigmaBar = zeros(n);
chi = zeros(n);
for i = 1:n
     for j = 1:n
         a = (t*F(j))/D(i);
         theRoot = sqrt((2*pi*b)/log(L/deltaX));
        sigmaBar(j,i) = .5*(1+erf(a*theRoot));
        chi(j,i) = (t/D(i))*sqrt(2*b/log(L/deltaX))*exp(-1*(2*pi*b*(t*F(j))^2)/(D(i)^2*log(L/deltaX)));
     end
end

figure;
surf(D,F,sigmaBar);


figure;
surf(D,F,chi);
zlim([0 5]);

figure;
contourf(D,F,sigmaBar,12,'lineStyle','none');
