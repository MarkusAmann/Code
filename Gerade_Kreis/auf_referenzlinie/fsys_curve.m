function f = fsys_curve(X,u,p)
ax = u;
sr = X(1);
v = X(2);
dr = X(3);
psir = X(4);

f = [v*cos(psir)/(1-dr*p.kapparef_curve);...
    ax;...
    v*sin(psir);...
    p.kapparef_curve*v - p.kapparef_curve*v*cos(psir)/(1-dr*p.kapparef_curve)];
end