clc
clear all
close all

%%
x0 = [0 1 0 0].'; l0 = [0.1 0 0 0].';
t0 = 0; tf = 100; N = 100; fx = 1; fy = 1; fr = 1; kapparef = 0.01; sf = 1000; drf = 0; psirf = 0;
tf_free = 0;
maxIter = 400;
maxFunEval = 1e5;
p.fx = fx; p.fy = fy; p.fr = fr; p.kapparef = kapparef; p.sf = sf; p.drf = drf; p.psirf = psirf;
p.x0 = x0; p.l0 = l0; p.t0 = t0; p.tf = tf; p.tf_free = tf_free; p.N = N;  
p.maxIter = maxIter; p.maxFunEval = maxFunEval;

bvpoptions = bvpset('RelTol',1e-5,'Stats','on');

switch tf_free
    case 0
        bvpoptions = bvpset(bvpoptions,'FJacobian',@jac_fix_tf);
        t = linspace(p.t0, p.tf, p.N+1);
        deltat = mean(diff(t));
        p.deltat = deltat;
        init_guess = @(x)guess_fix_tf(t,p);
        solinit = bvpinit(t,init_guess);
        sol = bvp4c(@bvpfcn_fix_tf, @bcfcn_fix_tf, solinit, bvpoptions, p);
        % optimal states
        sol_mesh = sol.x;
        sopt = sol.y(1,:);
        vopt = sol.y(2,:);
        dropt = sol.y(3,:);
        psiropt = sol.y(4,:);
        l1opt = sol.y(5,:);
        l2opt = sol.y(6,:);
        l3opt = sol.y(7,:);
        l4opt = sol.y(8,:);

    case 1
        bvpoptions = bvpset(bvpoptions,'FJacobian',@jac_free_tf);
        t = linspace(p.t0, p.tf, p.N+1);
        deltat = mean(diff(t));
        p.deltat = deltat;
        init_guess = @(x)guess_free_tf(t,p);
        solinit = bvpinit(t,init_guess,0.1);
        sol = bvp4c(@bvpfcn_free_tf, @bcfcn_free_tf, solinit, bvpoptions, p);
        % optimal states
        sol_mesh = sol.x;
        sopt = sol.y(1,:);
        vopt = sol.y(2,:);
        dropt = sol.y(3,:);
        psiropt = sol.y(4,:);
        l1opt = sol.y(5,:);
        l2opt = sol.y(6,:);
        l3opt = sol.y(7,:);
        l4opt = sol.y(8,:);
        delta_opt = sol.parameters;
        tf_opt = delta_opt*p.tf;
        sol_mesh = sol_mesh*delta_opt;
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
psiref = cumtrapz(sol_mesh,dot_psiref);
dot_psir = dot_psi - dot_psiref;
vref = dot_sopt;

% lateral acc., heading angle of car
ayopt = vopt.^2.*kappaopt;
psiopt = cumtrapz(sol_mesh,dot_psi); 

% coordinate transformation of car movement to global coordinates
dx_global_opt = vopt.*cos(psiopt);
dy_global_opt = vopt.*sin(psiopt);
x_global_opt = cumtrapz(sol_mesh,dx_global_opt);
y_global_opt = cumtrapz(sol_mesh,dy_global_opt);

% coordinate transformation of reference curve to global coordinates
dx_ref = vref.*cos(psiref);
dy_ref = vref.*sin(psiref);
x_ref = cumtrapz(sol_mesh,dx_ref);
y_ref = cumtrapz(sol_mesh,dy_ref);

%%
figure
subplot(3,1,1)
plot(sol_mesh,sopt)
ylabel('s_r [m]')
grid on
hold on
subplot(3,1,2)
plot(sol_mesh,vopt)
ylabel('v [m/s]')
grid on
hold on
subplot(3,1,3)
plot(sol_mesh,axopt)
ylabel('a_x_{opt} [m/s^2]')
xlabel('t [s]')
grid on
hold on

figure
subplot(2,1,1)
plot(sol_mesh, dropt)
hold on
grid on
ylabel('d_r_{opt} [m]')
subplot(2,1,2)
plot(sol_mesh, psiropt)
hold on
grid on
ylabel('psi_r_{opt} [rad]')
xlabel('t [s]')

figure
subplot(2,2,1)
plot(sol_mesh,l1opt)
ylabel('l_{1,opt}')
grid on
hold on
subplot(2,2,2)
plot(sol_mesh,l2opt)
ylabel('l_{2,opt}')
grid on
hold on
subplot(2,2,3)
plot(sol_mesh,l3opt)
ylabel('l_{3,opt}')
xlabel('t [s]')
grid on
hold on
subplot(2,2,4)
plot(sol_mesh,l4opt)
ylabel('l_{4,opt}')
xlabel('t [s]')
grid on
hold on

figure 
plot(sol_mesh,kappaopt)
ylabel('\kappa_{opt} [1/m]')
xlabel('t [s]')
grid on
hold on

figure 
plot(sol_mesh,ayopt)
ylabel('a_y [m/s^2]')
xlabel('t [s]')
grid on
hold on

figure 
plot(sol_mesh,psiopt)
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
plot(sol_mesh, x_global_opt)
hold on
grid on
ylabel('x position [m]')
subplot(2,1,2)
plot(sol_mesh, y_global_opt)
hold on
grid on
ylabel('y position [m]')
xlabel('t [s]')





