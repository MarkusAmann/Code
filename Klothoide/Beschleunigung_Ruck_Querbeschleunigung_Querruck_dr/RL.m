syms s v a dr psir kappa j dkappa l1 l2 l3 l4 l5 l6 fax fjx fay fjy fr kappar
% ohne Vereinfachung dr*kappar << 1
x_sym = [s v a dr psir kappa l1 l2 l3 l4 l5 l6];
u_sym = [j dkappa];
f_sym = [fax fjx fay fjy fr];
fax_num = 1; fjx_num = 1; fay_num = 1; fjy_num = 1; fr_num = 10;
f_num = [fax_num fjx_num fay_num fjy_num fr_num];
kappar_num = 0.01;

eqn_RL = 1/2*fay^2/fr*kappar^6*v^8-3/2*fay*kappar^2*v^4+1;
v_RL_all_sym = solve(eqn_RL == 0, v);
v_RL_all_sym_subs = double(subs(v_RL_all_sym,[f_sym kappar],[f_num kappar_num]));
eqn_RL_num = subs(eqn_RL,f_sym,f_num);
v_RL_all = solve(eqn_RL_num == 0, v);
v_RL_all_num = subs(v_RL_all,kappar,kappar_num);
v_RL_all_num = double(v_RL_all_num);
v_RL = min(v_RL_all_num(imag(v_RL_all_num) == 0 & real(v_RL_all_num) >= 0));
dr_RL = kappar_num^3*v_RL^4*fay_num/fr_num;
l1_RL = -2*fay_num*kappar_num^2*v_RL^3;
l5_RL = -fay_num*kappar_num*v_RL^3;
kappa_RL = kappar_num;
kappa_RL_corrected = kappar_num/(1-dr_RL*kappar_num);
ay_RL = kappa_RL_corrected*v_RL^2;
s_RL = 10;
a_RL = 0;
psir_RL = 0;
l2_RL = 0; l3_RL = 0; l4_RL = 0; l6_RL = 0;
j_RL = 0; 
X_RL = [s_RL v_RL a_RL dr_RL psir_RL kappa_RL l1_RL l2_RL l3_RL l4_RL l5_RL l6_RL];
eqns_sym = [v*cos(psir)/(1-dr*kappar);...
    a;...
    j;...
    v*sin(psir);...
    kappa*v - kappar*v*cos(psir)/(1-dr*kappar);...
    dkappa;...
    0;...
    -(2*fay*v^3*kappa^2 + 2*fjy*dkappa^2*v^3 + 4*fjy*kappa*a^2*v + 6*fjy*kappa*a*v^2*dkappa + l1*cos(psir)/(1-dr*kappar) + l4*sin(psir) + l5*kappa - l5*kappar*cos(psir)/(1-dr*kappar));...
    -(fax*a + 4*fjy*kappa*a*v^2 + 2*fjy*kappa*v^3*dkappa + l2);...
    -(fr*dr + kappar*l1*v*cos(psir)/(1-dr*kappar)^2 - kappar^2*l5*v*cos(psir)/(1-dr*kappar)^2);...
    -(l4*v*cos(psir) - l1*v*sin(psir)/(1-dr*kappar) + kappar*l5*v*sin(psir)/(1-dr*kappar));...
    -(fay*kappa*v^4 + 2*fjy*a^2*v^2 + 2*fjy*a*v^3*dkappa + l5*v)];

eqns_subs = subs(eqns_sym,[f_sym kappar],[f_num kappar_num]);

A = jacobian(eqns_subs,x_sym);
B = jacobian(eqns_subs,u_sym);
A_num = double(subs(A,[x_sym u_sym],[X_RL j_RL kappa_RL_corrected]));
B_num = double(subs(B,[x_sym u_sym],[X_RL j_RL kappa_RL_corrected]));
ew = eig(A_num);
ew(5:6)