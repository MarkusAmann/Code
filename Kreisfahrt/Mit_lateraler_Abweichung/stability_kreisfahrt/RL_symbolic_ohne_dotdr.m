% J = int(1/2*fx*ax^2 + 1/2*fy*kappa^2*v^4))
% Gleichungen ohne dotdr = 0, reduzierte Gleichungen
syms v dr psir ax kappa l1 l2 l3 l4 fx fy fp kappar
% ohne Vereinfachung dr*kappar << 1
eqns = [-l4/(fy*v^2) - kappar*v/(1-dr*kappar);...
    -(l4^2/(fy*v^3) + l3*psir);...
    -(l3*v + l4*v*psir/(1-dr*kappar)*(1-kappar^2));...
    -1/2*l4^2/(fy*v^2) + l3*v*psir + 1];
eqns_subs = subs(eqns, [fx; fy; kappar], [1;1;0.01]);
unknown = [v; dr; psir; l4; l3; kappa];
S = solve(eqns_subs == 0, unknown);


% J = int(1/2*fx*ax^2 + 1/2*fy*kappa^2*v^4))
% Gleichungen ohne dotdr = 0
syms v dr psir ax kappa l1 l2 l3 l4 fx fy fp kappar
% ohne Vereinfachung dr*kappar << 1
eqns = [kappa*v - kappar*v/(1-dr*kappar);...
    -(2*fy*kappa^2*v^3 + l1/(1-dr*kappar) + l3*psir + l4*kappa - l4*kappar/(1-dr*kappar));...
    -(l1 - kappar*l4);...
    -(l3*v + l4*v*psir/(1-dr*kappar) - l1*v*psir*kappar/(1-dr*kappar));...
    kappa + l4/(fy*v^3);...
    1/2*fy*kappa^2*v^4 + l1*v/(1-dr*kappar) + l3*v*psir + l4*kappa*v - l4*kappar*v/(1-dr*kappar) + 1];
eqns_subs = subs(eqns, [ax; fx; fy; fp; kappar], [0;1;1;1;0.01]);
unknown = [v; dr; psir; l1; l4; l3; kappa];
S = solve(eqns_subs == 0, unknown);
v = double(S.v);
idx_v_real_pos = find(imag(v) == 0 & v>0);
v_real_pos = v(idx_v_real_pos);
dr = double(S.dr);
dr_match_v = dr(idx_v_real_pos);
psir = double(S.psir);
psir_match_v = psir(idx_v_real_pos);
l1 = double(S.l1);
l1_match_v = l1(idx_v_real_pos);
l3 = double(S.l3);
l3_match_v = l3(idx_v_real_pos);
l4 = double(S.l4);
l4_match_v = l4(idx_v_real_pos);
kappa = double(S.kappa);
kappa_match_v = kappa(idx_v_real_pos);

syms v dr psir ax kappa l1 l2 l3 l4 fx fy fr kappar
% mit Vereinfachung dr*kappar << 1
eqns = [kappa*v - kappar*v;...
    -(2*fy*kappar^2*v^3 + l1 + l3*psir + l4*kappa - l4*kappar);...
    -(l1 - kappar*l4);...
    -(l3*v + l4*v*psir - l1*v*psir*kappar);...
    kappar + l4/(fy*v^3);...
    1/2*fy*kappa^2*v^4 + l1*v + l3*v*psir + l4*kappa*v - l4*kappar*v + 1];
eqns_subs = subs(eqns, [psir; ax; fx; fy; fr; kappar], [0;0;1;1;1;0.01]);
unknown = [v; dr; psir; l1; l4; kappa];
S = solve(eqns_subs == 0, unknown);
kappa = double(S.kappa);
v = double(S.v);
dr = double(S.dr);
l1 = double(S.l1);
l4 = double(S.l4);
