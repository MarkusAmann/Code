%%
x0 = [0 1 0 0]; l0 = [0.1 0.5 0.1 0.5]; % not used right now, l0 is chosen to be zeros
t0 = 0; tf = 100; N = 100; fx = 1; fy = 1; fr = 1; kapparef = 0.01; sf = 50; drf = 0; psirf = 0;
tf_free = 1;
maxIter = 400;
maxFunEval = 1e5;
p.fx = fx; p.fy = fy; p.fr = fr; p.kapparef = kapparef; p.sf = sf; p.drf = drf; p.psirf = psirf;
p.x0 = x0; p.l0 = l0; p.t0 = t0; p.tf = tf; p.tf_free = tf_free; p.N = N;  
p.maxIter = maxIter; p.maxFunEval = maxFunEval;
t = linspace(t0, tf, N+1);
deltat = mean(diff(t));
p.deltat = deltat;

switch tf_free
    case 0
        [Xopt,fval,exitflag,output] = indirect_discretization(p);
        % optimal states
        sopt = Xopt(1:8:end);
        vopt = Xopt(2:8:end);
        dropt = Xopt(3:8:end);
        psiropt = Xopt(4:8:end);
        l1opt = Xopt(5:8:end);
        l2opt = Xopt(6:8:end);
        l3opt = Xopt(7:8:end);
        l4opt = Xopt(8:8:end);

    case 1
        p.tf = 100;
        [Xopt,fval,exitflag,output] = indirect_discretization(p, tf_free);
        % optimal states
        sopt = Xopt(1:8:end-1);
        vopt = Xopt(2:8:end);
        dropt = Xopt(3:8:end);
        psiropt = Xopt(4:8:end);
        l1opt = Xopt(5:8:end);
        l2opt = Xopt(6:8:end);
        l3opt = Xopt(7:8:end);
        l4opt = Xopt(8:8:end);
        delta_opt = Xopt(end);
        tf_opt = delta_opt*p.tf
        t = linspace(t0, tf_opt, N+1);
end

%%
% optimal control inputs
axopt = -l2opt/fx;
kappaopt = -l4opt./(fy*vopt.^3);
% kappaopt(1) = interp1(t(2:end),kappaopt(2:end),0,'pchip','extrap');

% time-derivatives of optimal values and values of reference curve
dot_sopt = vopt.*cos(psiropt)./(1-dropt*kapparef);
dot_psi = kappaopt.*vopt;
dot_psiref = kapparef*dot_sopt;
psiref = cumtrapz(t,dot_psiref);
dot_psir = dot_psi - dot_psiref;
vref = dot_sopt;

% lateral acc., heading angle of car
ayopt = vopt.^2.*kappaopt;
psiopt = cumtrapz(t,dot_psi); 

% coordinate transformation of car movement to global coordinates
dx_global_opt = vopt.*cos(psiopt);
dy_global_opt = vopt.*sin(psiopt);
x_global_opt = cumtrapz(t,dx_global_opt);
y_global_opt = cumtrapz(t,dy_global_opt);

% coordinate transformation of reference curve to global coordinates
dx_ref = vref.*cos(psiref);
dy_ref = vref.*sin(psiref);
x_ref = cumtrapz(t,dx_ref);
y_ref = cumtrapz(t,dy_ref);

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
subplot(2,2,1)
plot(t,l1opt)
ylabel('l_{1,opt}')
grid on
hold on
subplot(2,2,2)
plot(t,l2opt)
ylabel('l_{2,opt}')
grid on
hold on
subplot(2,2,3)
plot(t,l3opt)
ylabel('l_{3,opt}')
xlabel('t [s]')
grid on
hold on
subplot(2,2,4)
plot(t,l4opt)
ylabel('l_{4,opt}')
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
plot(x_global_opt,y_global_opt)
grid on
hold on
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





