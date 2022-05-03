function f = fsys(X,u,p)
j = u(1);
kappa = u(2);
sr = X(1);
v = X(2);
a = X(3);
dr = X(4);
psir = X(5);

f = [v*cos(psir)/(1-dr*p.kapparef);...
    a;...
    j;...
    v*sin(psir);...
    kappa*v - p.kapparef*v*cos(psir)/(1-dr*p.kapparef)];
end