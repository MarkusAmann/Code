% clc
clear all
close all

%% Parameter
x0 = [0 5].'; l0 = [0 0].'; %l0 = 0.1*randn(4,1);
alim = 1; use_umax = 0;
umax = alim; umin = -alim;
t0 = 0; t1 = 1; tf = 2; N = 100; fx = 1; fy = 1; fr = 1; kapparef_straight = 0.0; kapparef_curve = 0.1; sf = 100; drf = 0; psirf = 0; s1 = 40; %s1 = sf-2*pi/4*1/kapparef_curve; % Strecke, nach der von Gerade auf Kreis umgeschaltet wird

p.use_umax = use_umax; p.umax = umax; p.umin = umin; p.fx = fx; p.fy = fy; p.fr = fr; p.kapparef_straight = kapparef_straight; p.kapparef_curve = kapparef_curve; 
p.sf = sf; p.drf = drf; p.psirf = psirf; p.s1 = s1; p.t1 = t1;
p.x0 = x0; p.l0 = l0; p.t0 = t0; p.tf = tf; p.N = N;  

%% Ruhelage 
v_RL = ((2/(3*p.fy*p.kapparef_curve^2))^(1/4));
l1_RL = -2*p.fy*p.kapparef_curve^2*v_RL^3;
ay_RL = v_RL^2*p.kapparef_curve;

%% Optimierung
bvpoptions = bvpset('RelTol',1e-5,'Stats','on','Nmax',5e4);

t0_1 = linspace(p.t0, p.t1, p.N/(tf/t1));
t1_f = linspace(p.t1, p.tf, p.N*(1-(t1/tf)));
t = [t0_1 t1_f];
deltat = mean(diff(t));
p.deltat = deltat;
init_guess = @(x,region)guess_free_tf(x,region,p);
start_inits = [0.1 0.5 12 13];
inits = start_inits;
error_flag = 1;
while error_flag
    try
        solinit = bvpinit(t,init_guess,inits); % [nu_tilde1, delta_t1, delta_t2, nu_tilde2]
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
l1opt = sol.y(3,:);
l2opt = sol.y(4,:);
nu_tilde1 = sol.parameters(1);
nu_tilde2 = sol.parameters(2);
delta_t1_opt = sol.parameters(3)
delta_t2_opt = sol.parameters(4)
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
axopt = u;
kappaopt = [ones(size(sol_mesh_1))*p.kapparef_straight ones(size(sol_mesh_2))*p.kapparef_curve];
J_fun = 1/2*p.fx*axopt.^2+1/2*p.fy*kappaopt.^2.*vopt.^4;
J = trapz(sol_mesh,J_fun) + tf_opt 

%%
% time-derivatives of optimal values and values of reference curve
kapparef_vec = [ones(size(sol_mesh_1))*p.kapparef_straight ones(size(sol_mesh_2))*p.kapparef_curve];
dot_sopt = vopt;
dot_psi = kappaopt.*vopt;
dot_psiref = kapparef_vec.*dot_sopt;
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
figure(1)
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

figure(2)
subplot(2,1,1)
plot(sol_mesh,l1opt)
ylabel('l_{1,opt}')
grid on
hold on
subplot(2,1,2)
plot(sol_mesh,l2opt)
ylabel('l_{2,opt}')
grid on
hold on

figure(3) 
plot(sol_mesh,kappaopt)
ylabel('\kappa_{opt} [1/m]')
xlabel('t [s]')
grid on
hold on

figure(4) 
plot(sol_mesh,ayopt)
ylabel('a_y [m/s^2]')
xlabel('t [s]')
grid on
hold on

figure(5)
plot(sol_mesh,psiopt)
ylabel('\psi_{opt} [rad]')
xlabel('t [s]')
grid on
hold on

figure(6)
plot(x_global_opt,y_global_opt)
grid on
hold on
plot(x_ref,y_ref,'k--')
ylabel('y position [m]')
xlabel('x position [m]')
legend('trajectory', 'reference')

figure(7)
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


%% Interpolation von ax_min und vmax s1 = [10 20 30 40 50 80 100 120 140 160 180 190] bei sf = 200, kapparef = 0.1, v0 = 2.5, fx = fy = 1
t1_opt_vec = [3.2162 5.45145 7.21096 8.72289 10.0552 13.4143 15.3119 17.0335 18.6164 20.0908 21.4764 22.1404];
axmin = [-0.692955 -1.02408 -1.17919 -1.27068 -1.33321 -1.4427 -1.48559 -1.51701 -1.5412 -1.56058 -1.57653 -1.58352];
t_vmax = [1.8269 2.8926 3.7507 4.50979 5.13018 6.84405 7.76848 8.6133 9.41375 10.1137 10.8509 11.1455];
vmax = [3.3369 4.1691 4.9007 5.54186 6.12171 7.60788 8.45722 9.22925 9.94209 10.6072 11.2335 11.5338];
s1_vec = [10 20 30 40 50 80 100 120 140 160 180 190];
