function f = kangln(X,u,p)

j = u(1);
kappa = u(2);
sr = X(1);
v = X(2);
a = X(3);
dr = X(4);
psir = X(5);
l1 = X(6);
l2 = X(7);
l3 = X(8);
l4 = X(9);
l5 = X(10);

switch p.use_dr
    case 1
        % Bestrafung von dr im Gütefunktional
        f = [0;...
            -(2*p.fy*v^3*kappa^2 + l1*cos(psir)/(1-dr*p.kapparef) + l4*sin(psir) + l5*kappa - l5*p.kapparef*cos(psir)/(1-dr*p.kapparef));...
            -(p.fa*a + l2);...
            -(p.fr*dr + p.kapparef*l1*v*cos(psir)/(1-dr*p.kapparef)^2 - p.kapparef^2*l5*v*cos(psir)/(1-dr*p.kapparef)^2);...
            -(l4*v*cos(psir) - l1*v*sin(psir)/(1-dr*p.kapparef) + p.kapparef*l5*v*sin(psir)/(1-dr*p.kapparef))];
    case 0
        % ohne Zustandsbestrafung im Gütefunktional
        f = [0;...
            -(2*p.fy*v^3*kappa^2 + l1*cos(psir)/(1-dr*p.kapparef) + l3*sin(psir) + l5*kappa - l5*p.kapparef*cos(psir)/(1-dr*p.kapparef));...
            -(p.fa*a + l2);...
            -(p.kapparef*l1*v*cos(psir)/(1-dr*p.kapparef)^2 - p.kapparef^2*l5*v*cos(psir)/(1-dr*p.kapparef)^2);...
            -(l3*v*cos(psir) - l1*v*sin(psir)/(1-dr*p.kapparef) + p.kapparef*l5*v*sin(psir)/(1-dr*p.kapparef))];
end
