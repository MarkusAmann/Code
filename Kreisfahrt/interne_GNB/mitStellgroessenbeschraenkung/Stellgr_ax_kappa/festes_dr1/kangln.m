function f = kangln(X,u,p)

ax = u(1);
kappa = u(2);
sr = X(1);
v = X(2);
dr = X(3);
psir = X(4);
l1 = X(5);
l2 = X(6);
l3 = X(7);
l4 = X(8);

switch p.use_dr
    case 1
        % Bestrafung von dr im Gütefunktional
        f = [0;...
        -(2*p.fy*v^3*kappa^2 + l1*cos(psir)/(1-dr*p.kapparef) + l3*sin(psir) + l4*kappa - l4*p.kapparef*cos(psir)/(1-dr*p.kapparef));...
        -(p.fr*dr + p.kapparef*l1*v*cos(psir)/(1-dr*p.kapparef)^2 - p.kapparef^2*l4*v*cos(psir)/(1-dr*p.kapparef)^2);...
        -(l3*v*cos(psir) - l1*v*sin(psir)/(1-dr*p.kapparef) + p.kapparef*l4*v*sin(psir)/(1-dr*p.kapparef))];
    case 0
        % ohne Zustandsbestrafung im Gütefunktional
        f = [0;...
        -(2*p.fy*v^3*kappa^2 + l1*cos(psir)/(1-dr*p.kapparef) + l3*sin(psir) + l4*kappa - l4*p.kapparef*cos(psir)/(1-dr*p.kapparef));...
        -(p.kapparef*l1*v*cos(psir)/(1-dr*p.kapparef)^2 - p.kapparef^2*l4*v*cos(psir)/(1-dr*p.kapparef)^2);...
        -(l3*v*cos(psir) - l1*v*sin(psir)/(1-dr*p.kapparef) + p.kapparef*l4*v*sin(psir)/(1-dr*p.kapparef))];
end
end
