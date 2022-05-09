function f = kangln(X,u,p,region)
j = u(1);
dkappa = u(2);
sr = X(1);
v = X(2);
a = X(3);
dr = X(4);
psir = X(5);
kappa = X(6);
l1 = X(7);
l2 = X(8);
l3 = X(9);
l4 = X(10);
l5 = X(11);
l6 = X(12);

switch region
    case 1 % Gerade vor dem Spurwechsel, mit Bestrafung von dr
        % Bestrafung von dr im Gütefunktional
                f = [0;...
            -(2*p.fay*v^3*kappa^2 + 2*p.fjy*dkappa^2*v^3 + 4*p.fjy*kappa^2*a^2*v + 6*p.fjy*kappa*a*v^2*dkappa + l1*cos(psir)/(1-dr*p.kapparef) + l4*sin(psir) + l5*kappa - l5*p.kapparef*cos(psir)/(1-dr*p.kapparef));...
            -(p.fax*a + 4*p.fjy*kappa^2*a*v^2 + 2*p.fjy*kappa*v^3*dkappa + l2);...
            -(p.fr*dr + p.kapparef*l1*v*cos(psir)/(1-dr*p.kapparef)^2 - p.kapparef^2*l5*v*cos(psir)/(1-dr*p.kapparef)^2);...
            -(l4*v*cos(psir) - l1*v*sin(psir)/(1-dr*p.kapparef) + p.kapparef*l5*v*sin(psir)/(1-dr*p.kapparef));...
            -(p.fay*kappa*v^4 + 4*p.fjy*a^2*v^2*kappa + 2*p.fjy*a*v^3*dkappa + l5*v)];
    case 2 % Gerade nach dem Spurwechsel, ohne Bestrafung von dr bzw. mit Bestrafung von dr-dr1
        % ohne Zustandsbestrafung im Gütefunktional
                f = [0;...
            -(2*p.fay*v^3*kappa^2 + 2*p.fjy*dkappa^2*v^3 + 4*p.fjy*kappa^2*a^2*v + 6*p.fjy*kappa*a*v^2*dkappa + l1*cos(psir)/(1-dr*p.kapparef) + l4*sin(psir) + l5*kappa - l5*p.kapparef*cos(psir)/(1-dr*p.kapparef));...
            -(p.fax*a + 4*p.fjy*kappa^2*a*v^2 + 2*p.fjy*kappa*v^3*dkappa + l2);...
            -(p.use_dr_min_dr1*(p.fr*dr - p.fr*p.dr1) + p.kapparef*l1*v*cos(psir)/(1-dr*p.kapparef)^2 - p.kapparef^2*l5*v*cos(psir)/(1-dr*p.kapparef)^2);...
            -(l4*v*cos(psir) - l1*v*sin(psir)/(1-dr*p.kapparef) + p.kapparef*l5*v*sin(psir)/(1-dr*p.kapparef));...
            -(p.fay*kappa*v^4 + 4*p.fjy*a^2*v^2*kappa + 2*p.fjy*a*v^3*dkappa + l5*v)];
end
end