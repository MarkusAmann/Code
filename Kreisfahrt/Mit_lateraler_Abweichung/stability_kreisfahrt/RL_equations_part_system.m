function F = RL_equations_part_system(x,params)
% x = (v,dr,l1,l4,kappa)
v = x(1);
dr = x(2);
l1 = x(3);
l4 = x(4);
kappa = x(5);
F(1) = kappa*v - params.kapparef*v/(1-dr*params.kapparef);
F(2) =  -2*params.fy*kappa^2*v^3 - l1/(1-dr*params.kapparef) - l4*kappa + l4*params.kapparef/(1-dr*params.kapparef);
F(3) = l1 - l4*v;
F(4) = kappa + l4/(params.fy*v^3);
F(5) = 1/2*params.fy*kappa^2*v^4 + l1*v/(1-dr*params.kapparef) + l4*kappa*v - l4*v*params.kapparef/(1-dr*params.kapparef) + 1;
end