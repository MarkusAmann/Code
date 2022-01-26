%% Ruhelage
syms vsym drsym l1sym l5sym fysym frsym kapparefsym imag_entry
eqns = [1/2*fysym*kapparefsym^2*vsym^4 + l1sym*vsym + 1/2*frsym*drsym^2 + 1;...
    l1sym + 2*fysym*kapparefsym^2*vsym^3;...
   kapparefsym*vsym*(kapparefsym*l5sym - l1sym) - frsym*drsym;...
    l5sym + fysym*kapparefsym*vsym^3];
RLs = solve(eqns == 0, [vsym, drsym, l1sym, l5sym]);
v_RL = RLs.vsym;
dr_RL = RLs.drsym;
l1_RL = RLs.l1sym;
l5_RL = RLs.l5sym;
RL_mat = [v_RL, dr_RL, l1_RL, l5_RL];
complex_check = double(subs(RL_mat,[kapparefsym, frsym, fysym],[p.kapparef, p.fr, p.fy]));
X_RL = complex_check(1,:);

%% Stabilität
syms ax j v psir dr kappa l1 l2 l3 l4 l5
j = -l3/p.fj; kappa = -l5/(p.fy*v^3);
dXdt = [ax;...
        j;...
        v*sin(psir);...
        kappa*v - p.kapparef*v*cos(psir)/(1-dr*p.kapparef);...
        0;...
        -(2*p.fy*kappa^2*v^3 + l1*cos(psir)/(1-dr*p.kapparef) + l4*sin(psir) + l5*kappa - l5*p.kapparef*cos(psir)/(1-dr*p.kapparef));...
        -l2;...
        -(l1*v*p.kapparef*cos(psir)/(1-dr*p.kapparef)^2 - l5*v*p.kapparef^2*cos(psir)/(1-dr*p.kapparef)^2 + p.fr*dr);...
        -(l4*v*cos(psir) - l1*v*sin(psir)/(1-dr*p.kapparef) + l5*v*p.kapparef*sin(psir)/(1-dr*p.kapparef))
        ];

A_sym = jacobian(dXdt,[v;ax;dr;psir;l1;l2;l3;l4;l5]);
A = subs(A_sym, [v;ax;dr;psir;l1;l2;l3;l4;l5], [X_RL(1); 0; X_RL(2); 0; X_RL(3); 0; 0; 0; X_RL(4)]);
A = double(A);
eigs = eig(A);
real_eig = real(eigs);
unstable_eig = real_eig > 0;
if sum(unstable_eig) > 0
    display("Die Systemmatrix des linearen Teils ist instabil.")
end