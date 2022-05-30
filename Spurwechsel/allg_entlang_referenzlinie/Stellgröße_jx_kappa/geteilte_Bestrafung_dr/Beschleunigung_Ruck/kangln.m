function f = kangln(X,u,p,region)
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

switch region
    case 1 % Gerade vor dem Spurwechsel, mit Bestrafung von dr
        % Bestrafung von dr im Gütefunktional
                 f = [0;...
                    -(2*p.fay*v^3*kappa^2 + l1*cos(psir)/(1-dr*p.kapparef) + l4*sin(psir) + l5*kappa - l5*p.kapparef*cos(psir)/(1-dr*p.kapparef));...
                    -(p.fax*a + l2);...
                    -(p.fr*dr + p.kapparef*l1*v*cos(psir)/(1-dr*p.kapparef)^2 - p.kapparef^2*l5*v*cos(psir)/(1-dr*p.kapparef)^2);...
                    -(l4*v*cos(psir) - l1*v*sin(psir)/(1-dr*p.kapparef) + p.kapparef*l5*v*sin(psir)/(1-dr*p.kapparef))];
    case 2 % Gerade nach dem Spurwechsel, ohne Bestrafung von dr bzw. mit Bestrafung von dr-dr1
        % ohne Zustandsbestrafung im Gütefunktional
                 f = [0;...
                    -(2*p.fay*v^3*kappa^2 + l1*cos(psir)/(1-dr*p.kapparef) + l4*sin(psir) + l5*kappa - l5*p.kapparef*cos(psir)/(1-dr*p.kapparef));...
                    -(p.fax*a + l2);...
                    -(p.use_dr_min_dr1*(p.fr*dr - p.fr*p.dr1) + p.kapparef*l1*v*cos(psir)/(1-dr*p.kapparef)^2 - p.kapparef^2*l5*v*cos(psir)/(1-dr*p.kapparef)^2);...
                    -(l4*v*cos(psir) - l1*v*sin(psir)/(1-dr*p.kapparef) + p.kapparef*l5*v*sin(psir)/(1-dr*p.kapparef))];
end
end