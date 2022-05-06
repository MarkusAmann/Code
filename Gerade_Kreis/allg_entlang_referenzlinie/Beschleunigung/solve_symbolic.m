syms sr(t) v(t) dr(t) psir(t) l1(t) l2(t) l3(t) l4(t) ax(t) kappa(t)
x0 = [0 1 0 0].'; l0 = [0 0 0 0].'; %l0 = 0.1*randn(4,1);
alim = 1.06*1000; kappalim = 1/4*1000; use_umax = 0;
umax = [alim;kappalim]; umin = -[alim;kappalim];
t0 = 0; tf = 10; N = 102; fx = 1; fy = 1; fr = 1; kapparef = 0.1; sf = pi/2*1/kapparef; drf = 0; psirf = 0; dr1 = 0.5; % dr1 ist der seitliche Versatz im Scheitelpunkt der Kurve
tf_free = 0; t1 = tf/3; 
p.use_umax = use_umax; p.umax = umax; p.umin = umin; p.fx = fx; p.fy = fy; p.fr = fr; p.kapparef = kapparef; p.sf = sf; p.drf = drf; p.psirf = psirf; p.dr1 = dr1;
p.x0 = x0; p.l0 = l0; p.t0 = t0; p.tf = tf; p.tf_free = tf_free; p.t1 = t1; p.N = N; 

ax = -l2/p.fx; kappa = -l4/(p.fy*v^3);
f_sys = [diff(sr,t) == v*cos(psir)/(1-dr*p.kapparef);...
    diff(v,t) == ax;...
    diff(dr,t) == v*sin(psir);...
    diff(psir,t) == kappa*v - p.kapparef*v*cos(psir)/(1-dr*p.kapparef)];
f_kan = [diff(l1,t) == 0;...
    diff(l2,t) == -(2*p.fy*v^3*kappa^2 + l1*cos(psir)/(1-dr*p.kapparef) + l3*sin(psir) + l4*kappa - l4*p.kapparef*cos(psir)/(1-dr*p.kapparef));...
    diff(l3,t) == -(p.fr*dr + p.kapparef*l1*v*cos(psir)/(1-dr*p.kapparef)^2 - p.kapparef^2*l4*v*cos(psir)/(1-dr*p.kapparef)^2);...
    diff(l4,t) == -(l3*v*cos(psir) - l1*v*sin(psir)/(1-dr*p.kapparef) + p.kapparef*l4*v*sin(psir)/(1-dr*p.kapparef))];
f = [f_sys; f_kan];

S = dsolve(f);