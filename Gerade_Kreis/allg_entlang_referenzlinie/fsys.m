function f = fsys(X,u,p,region)
ax = u(1);
kappa = u(2);
sr = X(1);
v = X(2);
dr = X(3);
psir = X(4);

switch region
    case 1
        f = [v*cos(psir)/(1-dr*p.kapparef_straight);...
            ax;...
            v*sin(psir);...
            kappa*v - p.kapparef_straight*v*cos(psir)/(1-dr*p.kapparef_straight)];
    case 2 
        f = [v*cos(psir)/(1-dr*p.kapparef_curve);...
            ax;...
            v*sin(psir);...
            kappa*v - p.kapparef_curve*v*cos(psir)/(1-dr*p.kapparef_curve)];
end
end