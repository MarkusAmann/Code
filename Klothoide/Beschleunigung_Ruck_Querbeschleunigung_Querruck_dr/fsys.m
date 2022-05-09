function f = fsys(X,u,p)
j = u(1);
dkappa = u(2);
sr = X(1);
v = X(2);
a = X(3);
dr = X(4);
psir = X(5);
kappa = X(6);

f = [v*cos(psir)/(1-dr*(p.kappa_*sr + p.kapparef0));...
    a;...
    j;...
    v*sin(psir);...
    kappa*v - (p.kappa_*sr + p.kapparef0)*v*cos(psir)/(1-dr*(p.kappa_*sr + p.kapparef0));...
    dkappa];
end