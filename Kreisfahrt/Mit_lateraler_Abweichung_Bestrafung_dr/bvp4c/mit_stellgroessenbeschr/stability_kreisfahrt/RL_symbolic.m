% J = int(1/2*fr*dr^2 + 1/2*fx*ax^2 + 1/2*fy*kappa^2*v^4))

syms v dr psir ax kappa l1 l2 l3 l4 fx fy fr kappar
% ohne Vereinfachung dr*kappar << 1
eqns = [kappa*v - kappar*v/(1-dr*kappar);...
    -(2*fy*kappa^2*v^3 + l1/(1-dr*kappar) + l3*psir + l4*kappa - l4*kappar/(1-dr*kappar));...
    -(fr*dr + kappar*l1*v/(1-dr*kappar)^2 - kappar^2*l4*v/(1-dr*kappar)^2);...
    kappa + l4/(fy*v^3);...
    1/2*fr*dr^2 + 1/2*fy*kappa^2*v^4 + l1*v/(1-dr*kappar) + l4*kappa*v - l4*kappar*v/(1-dr*kappar) + 1];
eqns_subs = subs(eqns, [psir; ax; fx; fy; fr; kappar], [0;0;1;1;1;0.01]);
unknown = [v; dr; l1; l4; kappa];
S = solve(eqns_subs == 0, unknown);
v = double(S.v);
dr = double(S.dr);
l1 = double(S.l1);
l4 = double(S.l4);
kappa = double(S.kappa);

% mit Vereinfachung dr*kappar << 1
% eqns = [-(2*fy*kappar^2*v^3 + l1);...
%     -(fr*dr + kappar*l1*v - kappar^2*l4*v);...
%     kappar + l4/(fy*v^3);...
%     1/2*fr*dr^2 + 1/2*fy*kappar^2*v^4 + l1*v + 1];
% eqns_subs = subs(eqns, [psir; ax; fx; fy; fr; kappar], [0;0;1;1;1;0.01]);
% unknown = [v; dr; l1; l4];
% S = solve(eqns_subs == 0, unknown);
%    Gleichung für dotl4 --> l3 = 0: -(l3*v - l1*v*psir/(1-dr*kappar) + l4*v*psir/(1-dr*kappar));...
