x0 = [0 0]; l0 = [0.1 0.5]; % not used right now, l0 is chosen to be zeros
t0 = 0; tf = 60; N = 500; fx = 1; fy = 1; kapparef = 0.02; sf = 200;
p.fx = fx; p.fy = fy; p.kappa = kapparef; p.sf = sf; p.x0 = x0; p.l0 = l0;
p.t0 = t0; p.tf = tf; p.N = N;
t = linspace(t0, tf, N+1);
[Xopt,fval,exitflag,output] = indirect_discretization(p);
sopt = Xopt(1:4:end);
vopt = Xopt(2:4:end);
l1opt = Xopt(3:4:end);
l2opt = Xopt(4:4:end);
axopt = -l2opt/(2*fx);
dot_psiopt = vopt*kapparef;
ayopt = dot_psiopt*kapparef;
psiopt = cumtrapz(dot_psiopt);
dx_global_opt = vopt.*cos(psiopt);
dy_global_opt = vopt.*sin(psiopt);
x_global_opt = cumtrapz(dx_global_opt);
y_global_opt = cumtrapz(dy_global_opt);

figure
subplot(3,1,1)
plot(t,sopt)
grid on
hold on
subplot(3,1,2)
plot(t,vopt)
grid on
hold on
subplot(3,1,3)
plot(t,axopt)
grid on
hold on

figure 
plot(t,ayopt)
grid on
hold on

figure 
plot(x_global_opt,y_global_opt)
grid on
hold on
