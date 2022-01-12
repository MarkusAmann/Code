%%
x0 = [0 0 0 0]; l0 = [0.1 0.5 0.1 0.5]; % not used right now, l0 is chosen to be zeros
t0 = 0; tf = 60; N = 100; fx = 1; fy = 1; kapparef = 0.02; sf = 50; drf = 0; psirf = 0;
maxIter = 400;
maxFunEval = 1e5;
p.fx = fx; p.fy = fy; p.kapparef = kapparef; p.sf = sf; p.drf = drf; p.psirf = psirf;
p.x0 = x0; p.l0 = l0; p.t0 = t0; p.tf = tf; p.N = N;  
p.maxIter = maxIter; p.maxFunEval = maxFunEval;
t = linspace(t0, tf, N+1);
deltat = mean(diff(t));
p.deltat = deltat;

[Xopt,fval,exitflag,output] = indirect_discretization(p);

%%
sopt = Xopt(1:8:end);
vopt = Xopt(2:8:end);
dropt = Xopt(3:8:end);
psiropt = Xopt(4:8:end);
l1opt = Xopt(5:8:end);
l2opt = Xopt(6:8:end);
l3opt = Xopt(7:8:end);
l4opt = Xopt(8:8:end);
axopt = -l2opt/fx;
kappaopt = -l4opt./(fy*vopt.^3);
% kappaopt(1) = interp1(t(2:end),kappaopt(2:end),0,'pchip','extrap');

dot_psiref = kapparef*vopt.*cos(psiropt)./(1-dropt*kapparef);
dot_psi = kappaopt.*vopt;
dot_psir = dot_psi - dot_psiref;
dot_psir_grad = gradient(psiropt, deltat);

ayopt = vopt.^2.*kappaopt;
psiopt = cumtrapz(dot_psi);
psiref = cumtrapz(dot_psiref);
dx_global_opt = vopt.*cos(psiopt);
dy_global_opt = vopt.*sin(psiopt);
x_global_opt = cumtrapz(dx_global_opt);
y_global_opt = cumtrapz(dy_global_opt);

%%
figure
subplot(3,1,1)
plot(t,sopt)
ylabel('s_r [m]')
grid on
hold on
subplot(3,1,2)
plot(t,vopt)
ylabel('v [m/s]')
grid on
hold on
subplot(3,1,3)
plot(t,axopt)
ylabel('a_x_{opt} [m/s^2]')
xlabel('t [s]')
grid on
hold on

figure
subplot(2,1,1)
plot(t, dropt)
hold on
grid on
ylabel('d_r_{opt} [m]')
subplot(2,1,2)
plot(t, psiropt)
hold on
grid on
ylabel('psi_r_{opt} [rad]')
xlabel('t [s]')


figure 
plot(t,ayopt)
ylabel('a_y [m/s^2]')
xlabel('t [s]')
grid on
hold on

figure 
plot(t,psiopt)
ylabel('\psi_{opt} [rad]')
xlabel('t [s]')
grid on
hold on

figure 
plot(t,kappaopt)
ylabel('\kappa_{opt} [1/m]')
xlabel('t [s]')
grid on
hold on

figure 
plot(x_global_opt,y_global_opt)
grid on
hold on
dx_ref = vopt.*cos(psiref);
dy_ref = vopt.*sin(psiref);
x_ref = cumtrapz(dx_ref);
y_ref = cumtrapz(dy_ref);
plot(x_ref,y_ref,'k--')
ylabel('y position [m]')
xlabel('x position [m]')
legend('trajectory', 'reference')

figure
subplot(2,1,1)
plot(t, x_global_opt)
hold on
grid on
ylabel('x position [m]')
subplot(2,1,2)
plot(t, y_global_opt)
hold on
grid on
ylabel('y position [m]')
xlabel('t [s]')





