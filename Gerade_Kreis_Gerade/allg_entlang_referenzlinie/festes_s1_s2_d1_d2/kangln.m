function f = kangln(X,u,p,region)
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

switch region
    case 1 % Gerade
        switch p.use_dr
            case 1
                % Bestrafung von dr im Gütefunktional
                f = [0;...
                -(2*p.fy*v^3*kappa^2 + l1*cos(psir)/(1-dr*p.kapparef_straight) + l3*sin(psir) + l4*kappa - l4*p.kapparef_straight*cos(psir)/(1-dr*p.kapparef_straight));...
                -(p.fr*dr + p.kapparef_straight*l1*v*cos(psir)/(1-dr*p.kapparef_straight)^2 - p.kapparef_straight^2*l4*v*cos(psir)/(1-dr*p.kapparef_straight)^2);...
                -(l3*v*cos(psir) - l1*v*sin(psir)/(1-dr*p.kapparef_straight) + p.kapparef_straight*l4*v*sin(psir)/(1-dr*p.kapparef_straight))];
            case 0
                % ohne Zustandsbestrafung im Gütefunktional
                f = [0;...
                -(2*p.fy*v^3*kappa^2 + l1*cos(psir)/(1-dr*p.kapparef_straight) + l3*sin(psir) + l4*kappa - l4*p.kapparef_straight*cos(psir)/(1-dr*p.kapparef_straight));...
                -(p.kapparef_straight*l1*v*cos(psir)/(1-dr*p.kapparef_straight)^2 - p.kapparef_straight^2*l4*v*cos(psir)/(1-dr*p.kapparef_straight)^2);...
                -(l3*v*cos(psir) - l1*v*sin(psir)/(1-dr*p.kapparef_straight) + p.kapparef_straight*l4*v*sin(psir)/(1-dr*p.kapparef_straight))];
        end
    case 2 % Kurve
        switch p.use_dr
            case 1
                 % Bestrafung von dr im Gütefunktional
                 f = [0;...
                 -(2*p.fy*v^3*kappa^2 + l1*cos(psir)/(1-dr*p.kapparef_curve) + l3*sin(psir) + l4*kappa - l4*p.kapparef_curve*cos(psir)/(1-dr*p.kapparef_curve));...
                 -(p.fr*dr + p.kapparef_curve*l1*v*cos(psir)/(1-dr*p.kapparef_curve)^2 - p.kapparef_curve^2*l4*v*cos(psir)/(1-dr*p.kapparef_curve)^2);...
                 -(l3*v*cos(psir) - l1*v*sin(psir)/(1-dr*p.kapparef_curve) + p.kapparef_curve*l4*v*sin(psir)/(1-dr*p.kapparef_curve))];
            case 0
                % ohne Zustandsbestrafung im Gütefunktional
                f = [0;...
                -(2*p.fy*v^3*kappa^2 + l1*cos(psir)/(1-dr*p.kapparef_curve) + l3*sin(psir) + l4*kappa - l4*p.kapparef_curve*cos(psir)/(1-dr*p.kapparef_curve));...
                -(p.kapparef_curve*l1*v*cos(psir)/(1-dr*p.kapparef_curve)^2 - p.kapparef_curve^2*l4*v*cos(psir)/(1-dr*p.kapparef_curve)^2);...
                -(l3*v*cos(psir) - l1*v*sin(psir)/(1-dr*p.kapparef_curve) + p.kapparef_curve*l4*v*sin(psir)/(1-dr*p.kapparef_curve))];
        end
    case 3 % Gerade
        switch p.use_dr
            case 1
                 % Bestrafung von dr im Gütefunktional
                 f = [0;...
                 -(2*p.fy*v^3*kappa^2 + l1*cos(psir)/(1-dr*p.kapparef_straight) + l3*sin(psir) + l4*kappa - l4*p.kapparef_straight*cos(psir)/(1-dr*p.kapparef_straight));...
                 -(p.fr*dr + p.kapparef_straight*l1*v*cos(psir)/(1-dr*p.kapparef_straight)^2 - p.kapparef_straight^2*l4*v*cos(psir)/(1-dr*p.kapparef_straight)^2);...
                 -(l3*v*cos(psir) - l1*v*sin(psir)/(1-dr*p.kapparef_straight) + p.kapparef_straight*l4*v*sin(psir)/(1-dr*p.kapparef_straight))];
            case 0
                % ohne Zustandsbestrafung im Gütefunktional
                f = [0;...
                -(2*p.fy*v^3*kappa^2 + l1*cos(psir)/(1-dr*p.kapparef_straight) + l3*sin(psir) + l4*kappa - l4*p.kapparef_straight*cos(psir)/(1-dr*p.kapparef_straight));...
                -(p.kapparef_straight*l1*v*cos(psir)/(1-dr*p.kapparef_straight)^2 - p.kapparef_straight^2*l4*v*cos(psir)/(1-dr*p.kapparef_straight)^2);...
                -(l3*v*cos(psir) - l1*v*sin(psir)/(1-dr*p.kapparef_straight) + p.kapparef_straight*l4*v*sin(psir)/(1-dr*p.kapparef_straight))];
        end
end
end