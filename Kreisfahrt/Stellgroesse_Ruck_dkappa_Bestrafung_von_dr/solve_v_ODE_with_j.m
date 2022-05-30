t0 = 0; tf = 1; N = 100; fjy = 1; fjx = 1; fax = 1; fay = 1; fr = 1; kapparef = 0.01; sf = 1000; drf = 0; psirf = 0;
x0 = [0 5 0 0 0 kapparef].'; l0 = [0 0 0 0 0 0].'; %l0 = 0.1*randn(4,1);
p.fjx = fjx; p.fjy = fjy; p.fax = fax; p.fay = fay; p.fr = fr; p.kapparef = kapparef; p.sf = sf; p.drf = drf; p.psirf = psirf;
p.x0 = x0; p.l0 = l0; p.t0 = t0; p.tf = tf; p.N = N;  

%%
[t,ode_sol] = ode45(@(t,X) v_Sol_ODE(t,X,p),[0 10],[x0(2); x0(3); 0.766098652067641; -0.844631416170913]);


function dXdt = v_Sol_ODE(t,X,p)
v = X(1);
dv = X(2);
d2v = X(3);
d3v = X(4);
v_RL = (2/(3*p.fay*p.kapparef^2))^(1/4);
l1_RL = -2*p.fay*p.kapparef^2*v_RL^3;

dx1 = dv;
dx2 = d2v;
dx3 = d3v;
dx4 = p.fax/p.fjx*d2v - 2*p.fay/p.fjx*p.kapparef^2*v^3 - l1_RL/p.fjx;

dXdt = [dx1; dx2; dx3; dx4];
end