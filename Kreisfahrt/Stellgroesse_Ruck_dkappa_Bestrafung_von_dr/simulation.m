% clc
clear all
close all

%% Parameter
jlim = 1.06*1000; dkappalim = 1/4*1000; use_umax = 0; use_dr = 1;
umax = [jlim;dkappalim]; umin = -[jlim;dkappalim];
t0 = 0; tf = 1; N = 100; fjy = 1; fjx = 1; fax = 1; fay = 1; fr = 1; kapparef = 0.01; sf = 500; drf = 0; psirf = 0;
x0 = [0 5 0 0.0067 0 kapparef].'; l0 = [0 0 0 0 0 0].'; %l0 = 0.1*randn(4,1);

p.use_umax = use_umax; p.umax = umax; p.umin = umin; p.use_dr = use_dr; p.fjx = fjx; p.fjy = fjy; p.fax = fax; p.fay = fay; p.fr = fr; p.kapparef = kapparef; p.sf = sf; p.drf = drf; p.psirf = psirf;
p.x0 = x0; p.l0 = l0; p.t0 = t0; p.tf = tf; p.N = N;  

%% Approximierte Lösung der Geschwindigkeit für das System ohne Querabweichung mit (v')^2 = fy/fx*kapparef^2*v^4-2*l1*v+cv \approx fy/fx*kapparef^2*v^4
syms fx fy kappar v(t) t C1 C2 l1
ode = diff(v,t)^2 == fy/fx*kappar^2*v^4 + 2*l1*v + C1; 
dsol = dsolve(ode);
dsol_num = subs(dsol,[fx fy kappar],[fjx, fay, kapparef]);
C1_num = -1/x0(2);
C2_num = -1/x0(2);
dsol_num = subs(dsol_num,[C1 C2],[C1_num C2_num]);

%% Ruhelage der Kreisfahrt mit Bestrafung von dr und 1-dr*kappar = 1
syms vsym
v_RL_all = solve(1/2*p.fay^2/p.fr*p.kapparef^6*vsym^8-3/2*p.fay*kapparef^2*vsym^4+1 == 0, vsym);
v_RL_all = double(v_RL_all);
v_RL = min(v_RL_all(imag(v_RL_all) == 0 & real(v_RL_all) >= 0));
v_RL_approx_auf_ref = (2/(3*p.fay*p.kapparef^2))^(1/4);
dr_RL = p.kapparef^3*v_RL^4*p.fay/p.fr;
l1_RL = -2*p.fay*p.kapparef^2*v_RL^3;
l5_RL = -p.fay*p.kapparef*v_RL^3;
% kappa_RL = p.kapparef;
% kappa_RL_corrected = p.kapparef/(1-dr_RL*p.kapparef);
% ay_RL = kappa_RL_corrected*v_RL^2;
%p.x0 = [0 v_RL 0 dr_RL 0].'; p.l0 = [l1_RL 0 0 0 l5_RL].';

%% Optimierung
bvpoptions = bvpset('RelTol',1e-5,'Stats','on');

t = linspace(p.t0, p.tf, p.N+1);
deltat = mean(diff(t));
p.deltat = deltat;
init_guess = @(x)guess_free_tf(t,p);
start_inits = [10];
inits = start_inits;
init_order = floor((log10(start_inits)>0).*log10(start_inits));
inits = 10.^(floor((log10(start_inits)>0).*log10(start_inits))).*abs(randn(size(start_inits))).*sign(start_inits);
error_flag = 1;
while error_flag
    try
        solinit = bvpinit(t,init_guess,inits); % [nu_tilde, delta_t1, delta_t2]
        sol = bvp4c(@sys_gesamt_free_tf, @bcfcn_free_tf, solinit, bvpoptions, p);
        error_flag = 0;
    catch ME
        if strcmp(ME.identifier,'MATLAB:bvp4c:SingJac')
            warn_message = strcat(ME.message, ' Reinitilization necessary.');
            warning(warn_message);
            error_flag = 1;
            init_order = floor((log10(start_inits)>0).*log10(start_inits));
            inits = 10.^(floor((log10(start_inits)>0).*log10(start_inits))).*abs(randn(size(start_inits))).*sign(start_inits);
        else
            error(ME.message)
        end
    end
end
% optimal states
sol_mesh = sol.x;
sopt = sol.y(1,:);
vopt = sol.y(2,:);
aopt = sol.y(3,:);
dropt = sol.y(4,:);
psiropt = sol.y(5,:);
kappaopt = sol.y(6,:);
l1opt = sol.y(7,:);
l2opt = sol.y(8,:);
l3opt = sol.y(9,:);
l4opt = sol.y(10,:);
l5opt = sol.y(11,:);
l6opt = sol.y(12,:);
delta_opt = sol.parameters;
tf_opt = delta_opt*p.tf
sol_mesh = sol_mesh*delta_opt;

%%
% optimal control inputs
for i=1:length(sol_mesh)
    u(:,i) = uopt(sol.y(:,i),p); % Steuerung
end
jopt = u(1,:);
dkappaopt = u(2,:);
jyopt = dkappaopt.*vopt.^2 + 2*kappaopt.*aopt.*vopt;
J_fun = 1/2*p.fjx*jopt.^2 + 1/2*p.fax*aopt.^2 + 1/2*p.fjy*jyopt.^2 + 1/2*p.fay*kappaopt.^2.*vopt.^4 + 1/2*p.fr*dropt.^2*p.use_dr;
J = trapz(sol_mesh,J_fun) + tf_opt 

dotjopt_numeric = gradient(jopt,sol_mesh);

%%
% time-derivatives of optimal values and values of reference curve
dot_sopt = vopt.*cos(psiropt)./(1-dropt*kapparef);
dot_psi = kappaopt.*vopt;
dot_psiref = kapparef*dot_sopt;
psiref = cumtrapz(sol_mesh,dot_psiref);
dot_psir = dot_psi - dot_psiref;
vref = dot_sopt;

% lateral acc., heading angle of car
ayopt = vopt.^2.*kappaopt;
jyopt = dkappaopt.*vopt.^2 + 2*kappaopt.*aopt.*vopt;
psiopt = cumtrapz(sol_mesh,dot_psi); 

% coordinate transformation of car movement to global coordinates
dx_global_opt = vopt.*cos(psiopt);
dy_global_opt = vopt.*sin(psiopt);
x_global_opt = cumtrapz(sol_mesh,dx_global_opt);
y_global_opt = cumtrapz(sol_mesh,dy_global_opt) + p.x0(4);

% coordinate transformation of reference curve to global coordinates
dx_ref = vref.*cos(psiref);
dy_ref = vref.*sin(psiref);
x_ref = cumtrapz(sol_mesh,dx_ref);
y_ref = cumtrapz(sol_mesh,dy_ref);

%%
set(gcf,'renderer','Painters')
figure(1)
subplot(2,2,1)
plot(sol_mesh,sopt,'-','Linewidth',2)
ylabel('$s_r\,[m]$','Interpreter','Latex')
xlabel('$t\,[s]$','Interpreter','Latex')
grid on
hold on
subplot(2,2,2)
plot(sol_mesh,vopt,'-','Linewidth',2)
ylabel('$v\,[m/s]$','Interpreter','Latex')
xlabel('$t\,[s]$','Interpreter','Latex')
grid on
hold on
subplot(2,2,3)
plot(sol_mesh,aopt,'-','Linewidth',2)
ylabel('$a_x\,[m/s^2]$','Interpreter','Latex')
xlabel('$t\,[s]$','Interpreter','Latex')
grid on
hold on
subplot(2,2,4)
plot(sol_mesh,jopt,'-','Linewidth',2)
ylabel('$j_x\,[m/s^3]$','Interpreter','Latex')
xlabel('$t\,[s]$','Interpreter','Latex')
grid on
hold on

figure(11)
subplot(2,2,1)
plot(sol_mesh, dropt,'-.','Linewidth',2)
hold on
grid on
xlabel('$t\,[s]$','Interpreter','Latex')
ylabel('$d_r [m]$','Interpreter','Latex')
subplot(2,2,2)
plot(sol_mesh, psiropt,'-.','Linewidth',2)
hold on
grid on
ylabel('$\psi_r [rad]$','Interpreter','Latex')
xlabel('$t\,[s]$','Interpreter','Latex')
subplot(2,2,3)
plot(sol_mesh,kappaopt,'-.','Linewidth',2)
ylabel('$\kappa [1/m]$','Interpreter','Latex')
xlabel('$t\,[s]$','Interpreter','Latex')
grid on
hold on
subplot(2,2,4)
plot(sol_mesh,dkappaopt,'-.','Linewidth',2)
ylabel('$\dot{\kappa} [1/(ms)]$','Interpreter','Latex')
xlabel('$t\,[s]$','Interpreter','Latex')
grid on
hold on

figure(22)
subplot(2,2,1)
plot(sol_mesh, dropt,'-.','Linewidth',2)
hold on
grid on
xlabel('$t\,[s]$','Interpreter','Latex')
ylabel('$d_r [m]$','Interpreter','Latex')
subplot(2,2,2)
plot(sol_mesh, psiropt,'-.','Linewidth',2)
hold on
grid on
ylabel('$\psi_r [rad]$','Interpreter','Latex')
xlabel('$t\,[s]$','Interpreter','Latex')
subplot(2,2,3)
plot(sol_mesh,kappaopt,'-.','Linewidth',2)
ylabel('$\kappa [1/m]$','Interpreter','Latex')
xlabel('$t\,[s]$','Interpreter','Latex')
grid on
hold on
subplot(2,2,4)
plot(sol_mesh,ayopt,'-.','Linewidth',2)
ylabel('$a_y [m/s^2]$','Interpreter','Latex')
xlabel('$t\,[s]$','Interpreter','Latex')
grid on
hold on


figure(2)
subplot(2,1,1)
plot(sol_mesh, dropt,'-.','Linewidth',2)
hold on
grid on
ylabel('d_r_{opt} [m]')
subplot(2,1,2)
plot(sol_mesh, psiropt,'-.','Linewidth',2)
hold on
grid on
ylabel('psi_r_{opt} [rad]')
xlabel('t [s]')

figure(3)
subplot(2,3,1)
plot(sol_mesh,l1opt)
ylabel('l_{1,opt}')
grid on
hold on
subplot(2,3,2)
plot(sol_mesh,l2opt)
ylabel('l_{2,opt}')
grid on
hold on
subplot(2,3,3)
plot(sol_mesh,l3opt)
ylabel('l_{3,opt}')
xlabel('t [s]')
grid on
hold on
subplot(2,3,4)
plot(sol_mesh,l4opt)
ylabel('l_{4,opt}')
xlabel('t [s]')
grid on
hold on
subplot(2,3,5)
plot(sol_mesh,l5opt)
ylabel('l_{5,opt}')
xlabel('t [s]')
grid on
hold on
subplot(2,3,6)
plot(sol_mesh,l6opt)
ylabel('l_{6,opt}')
xlabel('t [s]')
grid on
hold on

figure(4)
subplot(2,1,1)
plot(sol_mesh,kappaopt,'-.','Linewidth',2)
ylabel('\kappa_{opt} [1/m]')
xlabel('t [s]')
grid on
hold on
subplot(2,1,2)
plot(sol_mesh,dkappaopt,'-.','Linewidth',2)
ylabel('dkappa_{opt} [1/(m*s)]')
xlabel('t [s]')
grid on
hold on

figure(5)
subplot(2,1,1)
plot(sol_mesh,ayopt,'-','Linewidth',2)
ylabel('a_y [m/s^2]')
xlabel('t [s]')
grid on
hold on
subplot(2,1,2)
plot(sol_mesh,jyopt,'-','Linewidth',2)
ylabel('j_y [m/s^3)]')
xlabel('t [s]')
grid on
hold on

figure(44)
subplot(2,2,1)
plot(sol_mesh,kappaopt,'-','Linewidth',2)
ylabel('$\kappa\,[1/m]$','Interpreter','Latex')
xlabel('$t\,[s]$','Interpreter','Latex')
grid on
hold on
subplot(2,2,2)
plot(sol_mesh,gradient(kappaopt,sol_mesh),'-','Linewidth',2)
ylabel('$\dot{\kappa}\,[1/(ms)]$','Interpreter','Latex')
xlabel('$t\,[s]$','Interpreter','Latex')
grid on
hold on
subplot(2,2,3)
plot(sol_mesh,ayopt,'-','Linewidth',2)
ylabel('$a_y\,[m/s^2]$','Interpreter','Latex')
xlabel('$t\,[s]$','Interpreter','Latex')
grid on
hold on
subplot(2,2,4)
plot(sol_mesh,gradient(ayopt,sol_mesh),'-','Linewidth',2)
ylabel('$j_y\,[m/s^3)]$','Interpreter','Latex')
xlabel('$t\,[s]$','Interpreter','Latex')
grid on
hold on

figure(55)
subplot(3,1,1)
plot(sol_mesh,kappaopt,'-','Linewidth',2)
ylabel('$\kappa\,[1/m]$','Interpreter','Latex')
xlabel('$t\,[s]$','Interpreter','Latex')
grid on
hold on
subplot(3,1,2)
plot(sol_mesh,ayopt,'-','Linewidth',2)
ylabel('$a_y\,[m/s^2]$','Interpreter','Latex')
xlabel('$t\,[s]$','Interpreter','Latex')
grid on
hold on
subplot(3,1,3)
plot(sol_mesh,gradient(ayopt,sol_mesh),'-','Linewidth',2)
ylabel('$j_y\,[m/s^3)]$','Interpreter','Latex')
xlabel('$t\,[s]$','Interpreter','Latex')
grid on
hold on

figure(6)
plot(sol_mesh,psiopt)
ylabel('\psi_{opt} [rad]')
xlabel('t [s]')
grid on
hold on

figure(7)
plot(x_global_opt,y_global_opt,'-','Linewidth',2)
grid on
hold on
% plot(x_ref,y_ref,'k--','Linewidth',1.5)
ylabel('y_g-Position [m]')
xlabel('x_g-Position [m]')
% legend('trajectory', 'reference')

figure(8)
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