%% Ruhelage
fun = @(x)RL_equations_whole_system(x,p);
% x = (v,dr,psir,l1,l2,l3,l4,ax,kappa) für RL_equations_whole_system
x0 = [10 10 0 -1 0 0 -5 0 0.01];

% fun = @(x)RL_equations_part_system(x,p);
% x = (v,dr,l1,l4,kappa) für RL_equations_part_system
% x0 = [10 10 -1 -1 0.01];

options = optimoptions('fsolve','MaxIterations',1e5,'MaxFunctionEvaluations',1e5);
[x, fval] = fsolve(fun,x0,options);
v_RL = x(1);
dr_RL = x(2);
psir_RL = x(3);
l1_RL = x(4);
l2_RL = x(5);
l3_RL = x(6);
l4_RL = x(7);
ax_RL = x(8);
kappa_RL = x(9);
ax_RL = x(8);
kappa_RL = x(9);
X_RL = [v_RL; dr_RL; psir_RL; l1_RL; l2_RL; l3_RL; l4_RL];

%% Stabilität
syms v dr psir l1 l2 l3 l4 ax kappa
dXdt = [ax;...
        v*sin(psir);...
        kappa*v - p.kapparef*v*cos(psir)/(1-dr*p.kapparef);...
        0;...
        -(2*p.fy*kappa^2*v^3 + l1*cos(psir)/(1-dr*p.kapparef) + l3*sin(psir) + l4*kappa - l4*p.kapparef*cos(psir)/(1-dr*p.kapparef));...
        -(l1*v*p.kapparef*cos(psir)/(1-dr*p.kapparef)^2 - l4*v*p.kapparef^2*cos(psir)/(1-dr*p.kapparef)^2);...
        -(l3*v*cos(psir) - l1*v*sin(psir)/(1-dr*p.kapparef) + l4*v*p.kapparef*sin(psir)/(1-dr*p.kapparef))
        ];

A_sym = jacobian(dXdt,[v;dr;psir;l1;l2;l3;l4]);
A = subs(A_sym, [v;dr;psir;l1;l2;l3;l4;ax;kappa], x.');
A = double(A);
eigs = eig(A);
real_eig = real(eigs);
unstable_eig = real_eig > 0;
if sum(unstable_eig) > 0
    display("Die Systemmatrix des linearen Teils ist instabil.")
end