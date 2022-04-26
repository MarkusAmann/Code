% close all; 
clear;

syms c1_sym k1_sym fa_sym fj_sym sf_sym s0_sym v0_sym vf_sym a0_sym af_sym tf_sym;
var_sym = [fa_sym fj_sym sf_sym s0_sym v0_sym vf_sym a0_sym af_sym tf_sym];

%% hier parameter vorgeben
fa = 1;            % Gewichtung Längsbeschleunigung
fj = 1;            % Gewichtung Ruck.

sf = 200;   % Länge der Gerade  

s0 = 0; v0 = 10; vf = 0; a0 = 0; af = 0; tf = 40;
var_num = [fa fj sf s0 v0 vf a0 af tf];

c2_sym = c1_sym*tf_sym;
cs_sym = s0_sym - a0_sym*fj_sym/fa_sym - c1_sym*tf_sym*fj_sym/fa_sym^2;
cv_sym = v0_sym + a0_sym*sqrt(fj_sym/fa_sym) + c1_sym*tf_sym*sqrt(fj_sym/fa_sym^3);
eqns = [k1_sym*sqrt(fa_sym/fj_sym)*exp(sqrt(fa_sym/fj_sym)*tf_sym) + c1_sym/fa_sym;...
    k1_sym*fj_sym/fa_sym*exp(sqrt(fa_sym/fj_sym)*tf_sym) + c1_sym*tf_sym^3/(6*fa_sym) - c2_sym*tf_sym^2/(2*fa_sym) + cv_sym*tf_sym + cs_sym - sf_sym];
S = solve(eqns == 0,[c1_sym, k1_sym]);
c1 = double(subs(S.c1_sym,var_sym,var_num));
k1 = double(subs(S.k1_sym,var_sym,var_num));
c2 = double(subs(c2_sym,[var_sym c1_sym k1_sym],[var_num c1 k1]));
k2 = a0 + c2/fa;

tmax = -sqrt(fj/fa)*log(c1/k2*sqrt(fj/fa^3));
tmin = sqrt(fj/fa)*log(-c1/k1*sqrt(fj/fa^3));
amax = k2*exp(-sqrt(fa/fj)*tmax) + c1/fa*tmax - c2/fa;
amin = k1*exp(sqrt(fa/fj)*tmin) + c1/fa*tmin - c2/fa;