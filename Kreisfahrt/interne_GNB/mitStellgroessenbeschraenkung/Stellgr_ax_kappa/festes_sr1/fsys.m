function f = fsys(X,u,p)
ax = u(1);
kappa = u(2);
sr = X(1);
v = X(2);
dr = X(3);
psir = X(4);

f = [v*cos(psir)/(1-dr*p.kapparef);...
    ax;...
    v*sin(psir);...
    kappa*v - p.kapparef*v*cos(psir)/(1-dr*p.kapparef)];
end