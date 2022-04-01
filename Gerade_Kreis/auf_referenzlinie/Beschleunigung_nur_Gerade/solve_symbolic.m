syms fx fy kapparef s0 s1 sf v0 % bekannte symbolische Variablen
syms t1 tf nu_1 nu_2 l1l a0 % unbekannte symbolische Variablen
known_vars = [fx fy kapparef s0 s1 sf v0];
known_vars_num = [1 1 0.1 0 40 100 5];
unknwon_vars = [t1 tf nu_1 nu_2 l1l a0];

vr = (2/(3*fy*kapparef^2))^(1/4);
l1r = -2*fy*kapparef^2*vr^3;
ax_t1 = l1l*t1/fx + a0;
v_t1 = l1l*t1^2/(2*fx) + a0*t1 + v0;
l2l_t1 = -ax_t1*fx;

eqns = [l1l - l1r - 2*nu_1;...
    l1l*t1 + fx*a0 + 2*nu_2;...
    l1l*t1^3/(6*fx) + a0/2*t1^2 + v0*t1 + s0 - s1;...
    v_t1 - vr;...
    tf - t1 + (s1 - sf)/vr;...
    1/2*fx*ax_t1^2 + l1l*v_t1 + l2l_t1*ax_t1 - 1/2*fy*kapparef^2*vr^4 - l1r*vr];
eqns_num = subs(eqns,known_vars,known_vars_num);

S = solve(eqns_num, unknwon_vars);
t1 = vpa(S.t1,10);
tf = vpa(S.tf,10);
nu_1 = vpa(S.nu_1,10);
nu_2 = vpa(S.nu_2,10);
l1l = vpa(S.l1l,10);
a0 = vpa(S.a0,10);
t1_val = double(t1(1));
tf_val = double(tf(1));
nu_1_val = double(nu_1(1));
nu_2_val = double(nu_2(1));
l1l_val = double(l1l(1));
a0_val = double(a0(1));

t1_vec = linspace(0,t1_val,1000);
t2_vec = linspace(t1_val,tf_val,1000);
t_vec = [t1_vec t2_vec];
l1l_vec = l1l_val*ones(size(t1_vec));
l1r_vec = l1l_val*ones(size(t2_vec)) - 2*nu_1_val;
l1_vec = [l1l_vec l1r_vec];
l2l_vec = -l1l_vec.*t1_vec - known_vars_num(1)*a0_val;
l2r_vec = zeros(size(t2_vec));
l2_vec = [l2l_vec l2r_vec];
ax_vec = - l2_vec/known_vars_num(1);
v_vec = cumtrapz(t_vec,ax_vec) + known_vars_num(7);
s_vec = cumtrapz(t_vec,v_vec) + known_vars_num(4);
kappa_vec = [zeros(size(t1_vec)) known_vars_num(3)*ones(size(t2_vec))];

J_fun = 1/2*known_vars_num(1)*ax_vec.^2 + 1/2*known_vars_num(2)*kappa_vec.^2.*v_vec.^4;
J = trapz(t_vec,J_fun) + tf_val