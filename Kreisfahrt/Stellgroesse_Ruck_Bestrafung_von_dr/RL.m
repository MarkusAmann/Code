syms s v a dr psir j kappa l1 l2 l3 l4 l5 fa fj fy fr kappar
% ohne Vereinfachung dr*kappar << 1
x_sym = [s v a dr psir l1 l2 l3 l4 l5];
u_sym = [j kappa];
f_sym = [fa fj fy fr];
fa_num = 1; fj_num = 1; fy_num = 1; fr_num = 10;
f_num = [fa_num fj_num fy_num fr_num];
kappar_num = 0.01;

eqn_RL = 1/2*fy^2/fr*kappar^6*v^8-3/2*fy*kappar^2*v^4+1;
v_RL_all_sym = solve(eqn_RL == 0, v);
v_RL_all_sym_subs = double(subs(v_RL_all_sym,[f_sym kappar],[f_num kappar_num]));
eqn_RL_num = subs(eqn_RL,f_sym,f_num);
v_RL_all = solve(eqn_RL_num == 0, v);
v_RL_all_num = subs(v_RL_all,kappar,kappar_num);
v_RL_all_num = double(v_RL_all_num);
v_RL = min(v_RL_all_num(imag(v_RL_all_num) == 0 & real(v_RL_all_num) >= 0));
dr_RL = kappar_num^3*v_RL^4*fy_num/fr_num;
l1_RL = -2*fy_num*kappar_num^2*v_RL^3;
l5_RL = -fy_num*kappar_num*v_RL^3;
kappa_RL = kappar_num;
kappa_RL_corrected = kappar_num/(1-dr_RL*kappar_num);
ay_RL = kappa_RL_corrected*v_RL^2;
s_RL = 10;
a_RL = 0;
psir_RL = 0;
l2_RL = 0; l3_RL = 0; l4_RL = 0; 
j_RL = 0; 
X_RL = [s_RL v_RL a_RL dr_RL psir_RL l1_RL l2_RL l3_RL l4_RL l5_RL];
eqns_sym = [v*cos(psir)/(1-dr*kappar);...
    a;...
    j;...
    v*sin(psir);...
    kappa*v - kappar*v*cos(psir)/(1-dr*kappar);...
    0;...
    -(2*fy*v^3*kappa^2 + l1*cos(psir)/(1-dr*kappar) + l4*sin(psir) + l5*kappa - l5*kappar*cos(psir)/(1-dr*kappar));...
    -(fa*a + l2);...
    -(fr*dr + kappar*l1*v*cos(psir)/(1-dr*kappar)^2 - kappar^2*l5*v*cos(psir)/(1-dr*kappar)^2);...
    -(l4*v*cos(psir) - l1*v*sin(psir)/(1-dr*kappar) + kappar*l5*v*sin(psir)/(1-dr*kappar))];

eqns_subs = subs(eqns_sym,[f_sym kappar],[f_num kappar_num]);

A = jacobian(eqns_subs,x_sym);
B = jacobian(eqns_subs,u_sym);
A_num = double(subs(A,[x_sym u_sym],[X_RL j_RL kappa_RL_corrected]));
B_num = double(subs(B,[x_sym u_sym],[X_RL j_RL kappa_RL_corrected]));
ew = eig(A_num);
ew(5:6)