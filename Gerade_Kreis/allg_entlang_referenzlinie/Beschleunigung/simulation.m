% clc
clear all
% close all

%% Parameter
x0 = [0 1 0 0].'; l0 = [0 0 0 0].'; %l0 = 0.1*randn(4,1);
alim = 2; kappalim = 0.2*100; use_umax = 0;
umax = [alim;kappalim]; umin = -[alim;kappalim]; use_dr = 1;
t0 = 0; tf = 2; N = 100; fx = 1; fy = 1; fr = 1; kapparef_straight = 0.0; kapparef_curve = 0.01; sf = 700; drf = 0; psirf = 0; t1 = tf/2; s1 = 200; %s1 = sf-2*pi/4*1/kapparef_curve; % Strecke, nach der von Gerade auf Kreis umgeschaltet wird

p.use_umax = use_umax; p.umax = umax; p.umin = umin; p.fx = fx; p.fy = fy; p.fr = fr; p.kapparef_straight = kapparef_straight; p.kapparef_curve = kapparef_curve; 
p.sf = sf; p.drf = drf; p.psirf = psirf; p.s1 = s1; p.t1 = t1;
p.x0 = x0; p.l0 = l0; p.t0 = t0; p.tf = tf; p.N = N; p.use_dr = use_dr;

%% Ruhelage für Kreisfahrt mit Querdynamik und Bestrafung von dr
v_RL = (2/(3*p.fy*p.kapparef_curve^2))^(1/4);
dr_RL = (p.fy*p.kapparef_curve^3*v_RL^4)/(p.fr);
l1_RL = -2*p.fy*p.kapparef_curve^2*v_RL;
l4_RL = -p.fy*p.kapparef_curve*v_RL^3;
kappa_RL = p.kapparef_curve/(1-dr_RL*p.kapparef_curve);
ay_RL = kappa_RL*v_RL^2;

%% Optimierung
bvpoptions = bvpset('RelTol',1e-5,'Stats','on','Nmax',5e4);

t0_1 = linspace(p.t0, p.t1, p.N/(tf/t1));
t1_f = linspace(p.t1, p.tf, p.N*(1-(t1/tf)));
t = [t0_1 t1_f];
deltat = mean(diff(t));
p.deltat = deltat;
init_guess = @(x,region)guess_free_tf(x,region,p);
start_inits = [-0.1 12 13];
inits = start_inits;
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
dropt = sol.y(3,:);
psiropt = sol.y(4,:);
l1opt = sol.y(5,:);
l2opt = sol.y(6,:);
l3opt = sol.y(7,:);
l4opt = sol.y(8,:);
nu_tilde = sol.parameters(1);
delta_t1_opt = sol.parameters(2)
delta_t2_opt = sol.parameters(3)
t1_opt = delta_t1_opt*p.t1
tf_opt = delta_t1_opt*p.t1 + delta_t2_opt*(p.tf - p.t1)
split_idx = [find(diff(sol_mesh)==0) find(diff(sol_mesh)==0)+1];
sol_mesh_1 = sol_mesh(1:split_idx(1))*delta_t1_opt;
sol_mesh_2 = delta_t1_opt*p.t1 + delta_t2_opt*(sol_mesh(split_idx(2):end) - p.t1);
sol_mesh = [sol_mesh_1 sol_mesh_2];

%%
% optimal control inputs
for i=1:length(sol_mesh)
    u(:,i) = uopt(sol.y(:,i),p); % Steuerung
end
axopt = u(1,:);
kappaopt = u(2,:);
J_fun = 1/2*p.fx*axopt.^2 + 1/2*p.fy*kappaopt.^2.*vopt.^4 + 1/2*p.fr*dropt.^2*p.use_dr;
J = trapz(sol_mesh,J_fun) + tf_opt 

%%
% time-derivatives of optimal values and values of reference curve
kapparef_vec = [ones(size(sol_mesh_1))*p.kapparef_straight ones(size(sol_mesh_2))*p.kapparef_curve];
dot_sopt = vopt.*cos(psiropt)./(1-dropt.*kapparef_vec);
dot_psi = kappaopt.*vopt;
dot_psiref = kapparef_vec.*dot_sopt;
psiref = cumtrapz(sol_mesh,dot_psiref);
dot_psir = dot_psi - dot_psiref;
vref = dot_sopt;
ayref = vref.^2.*kapparef_vec;

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
figure(1)
subplot(3,1,1)
plot(sol_mesh,sopt)
ylabel('s_r [m]')
grid on
hold on
subplot(3,1,2)
plot(sol_mesh,vopt,sol_mesh,dot_sopt)
ylabel('v [m/s]')
legend('v_{opt}','ds_r')
grid on
hold on
subplot(3,1,3)
plot(sol_mesh,axopt)
ylabel('a_x_{opt} [m/s^2]')
xlabel('t [s]')
grid on
hold on

figure(2)
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

figure(3)
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

figure(4)
plot(sol_mesh,kappaopt,sol_mesh,kapparef_vec)
ylabel('\kappa [1/m]')
xlabel('t [s]')
legend('\kappa_{opt}','\kappa_{ref}')
grid on
hold on

figure(5)
plot(sol_mesh,ayopt,sol_mesh,ayref)
ylabel('a_y [m/s^2]')
xlabel('t [s]')
legend('ay_{opt}','ay_{ref}')
grid on
hold on

figure(6) 
plot(sol_mesh,psiopt)
ylabel('\psi_{opt} [rad]')
xlabel('t [s]')
grid on
hold on

figure(7)
plot(x_global_opt,y_global_opt)
grid on
hold on
plot(x_ref,y_ref,'k--')
ylabel('y position [m]')
xlabel('x position [m]')
legend('trajectory', 'reference')

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





