function f = fsys(X,u,p,region)
j = u(1);
kappa = u(2);
sr = X(1);
v = X(2);
a = X(3);
dr = X(4);
psir = X(5);

% Gerade mit kapparef = 0
f = [v*cos(psir)/(1-dr*p.kapparef);...
a;...
j;...
v*sin(psir);...
kappa*v - p.kapparef*v*cos(psir)/(1-dr*p.kapparef)];
  
end