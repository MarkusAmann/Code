syms c1 k1 k2 tf fa fj t sf a0 s0 v0
c1 = (sf+fj/fa*a0-s0-a0*tf*sqrt(fj/fa)-v0*tf)/(tf^2*sqrt(fj/fa^3)-tf*fj/fa^2-2*tf^3/(6*fa)-sqrt(fj^3/fa^5));
k1 = (-c1*sqrt(fj/fa^3))/(exp(sqrt(fa/fj)*tf));
c2 = -af*fa + c1*tf + k1*fa*exp(sqrt(fa/fj)*tf);
k2 = a0 + c1*tf/fa;
j = k1*sqrt(fa/fj)*exp(sqrt(fa/fj)*t) - k2*sqrt(fa/fj)*exp(-sqrt(fa/fj)*t) + c1/fa;
tf_num = 25; fa_num = 0.01; fj_num = 1; s0_num = 0; v0_num = 0.1; a0_num = 0; sf_num = 1000;
c1_num = double(subs(c1,[fa fj s0 v0 a0 sf tf],[fa_num fj_num s0_num v0_num a0_num sf_num tf_num]));
c2_num = double(subs(c2,[fa fj s0 v0 a0 sf tf],[fa_num fj_num s0_num v0_num a0_num sf_num tf_num]));
k1_num = double(subs(k1,[fa fj s0 v0 a0 sf tf],[fa_num fj_num s0_num v0_num a0_num sf_num tf_num]));
k2_num = double(subs(k2,[fa fj s0 v0 a0 sf tf],[fa_num fj_num s0_num v0_num a0_num sf_num tf_num]));

j_of_t = subs(j,[fa fj s0 v0 a0 sf tf],[fa_num fj_num s0_num v0_num a0_num sf_num tf_num]);
j_vec = subs(j_of_t,t,[0:0.01:tf_num]);