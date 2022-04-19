clear;

syms fj fy kapparef s0 s1 sf v0 a0 a1 % bekannte symbolische Variablen
syms t1 c1 c2 c3 % unbekannte symbolische Variablen
known_vars = [fj fy kapparef s0 s1 sf v0 a0 a1];
known_vars_num = [1 2 0.1 0 40 100 5 0 0];
% known_vars_num(7) = ((2/(3*known_vars_num(2)*known_vars_num(3)^2))^(1/4));
% known_vars_num(7) = 5;

unknwon_vars = [t1 c1 c2 c3];

l1 = c1;
l2_t1 = -c1*t1 + c2;
j_t1 = -c1/(2*fj)*t1^2 + c2/fj*t1 - c3/fj;
ax_t1 = -c1/(6*fj)*t1^3 + c2/(2*fj)*t1^2 - c3/fj*t1 + a0;
v_t1 = -c1/(24*fj)*t1^4 + c2/(6*fj)*t1^3 - c3/(2*fj)*t1^2 + a0*t1 + v0;
s_t1 = -c1/(120*fj)*t1^5 + c2/(24*fj)*t1^4 - c3/(6*fj)*t1^3 + a0/2*t1^2 + v0*t1 + s0;

eqns = [s_t1 - s1;...
    ax_t1 - a1;...
    l2_t1 - (sf-s1)*3/2*fy*v_t1^2*kapparef^2 + (sf-s1)/v_t1^2;
    -1/2*fj*j_t1^2 + l1*v_t1 + l2_t1*ax_t1 + 1];
eqns_num = subs(eqns,known_vars,known_vars_num);

% dd = vpasolve(eqns, unknwon_vars, [15,0,0,0]);
S = solve(eqns_num, unknwon_vars);
t1 = vpa(S.t1,10);
c1 = vpa(S.c1,10);
c2 = vpa(S.c2,10);
c3 = vpa(S.c3,10);
t1_val = double(t1);
c1_val = double(c1);
c2_val = double(c2);
c3_val = double(c3);

pos = find((imag(t1_val)==0)&(real(t1_val)>0));
if length(pos)>1
    max   = inf;
    p_max = 0; 
    for i = 1:length(pos)
        if t1_val(pos(i)) < max
            max   = t1_val(pos(i));
            p_max = pos(i);
        end
    end
    pos = p_max;
end