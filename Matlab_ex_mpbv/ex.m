clear all
close all

xc = 1;
xmesh = [0 0.25 0.5 0.75 xc xc 1.25 1.5 1.75 2];
yinit = [1; 1];
sol = bvpinit(xmesh,yinit);

lambda = 2;
n = 5e-2;
for kappa = 2:5
   eta = lambda^2/(n*kappa^2);
   p = [n kappa lambda eta];
   sol = bvp4c(@(x,y,r) f(x,y,r,p), @bc, sol);
   K2 = lambda*sinh(kappa/lambda)/(kappa*cosh(kappa));
   approx = 1/(1 - K2);
   computed = 1/sol.y(1,end);
   fprintf('  %2i    %10.3f     %10.3f \n',kappa,computed,approx);
end

plot(sol.x,sol.y(1,:),'--o',sol.x,sol.y(2,:),'--o')
line([1 1], [0 2], 'Color', 'k')
legend('v(x)','C(x)')
title('A Three-Point BVP Solved with bvp5c')
xlabel({'x', '\lambda = 2, \kappa = 5'})
ylabel('v(x) and C(x)')

%-------------------------------------------
function dydx = f(x,y,region,p) % equations being solved
n = p(1);
eta = p(4);

dydx = zeros(2,1);
dydx(1) = (y(2) - 1)/n;

switch region
    case 1    % x in [0 1]
        dydx(2) = (y(1)*y(2) - x)/eta;
    case 2    % x in [1 lambda]
        dydx(2) = (y(1)*y(2) - 1)/eta;
end
end
%-------------------------------------------
function res = bc(XL,XR) % boundary conditions
res = [XL(1,1)               % v(0) = 0
       XR(1,1) - XL(1,2)     % Continuity of v(x) at x=1
       XR(2,1) - XL(2,2)     % Continuity of C(x) at x=1
       XR(2,end) - 1];       % C(lambda) = 1
end
%-------------------------------------------