% clc
use_solution_as_init = 1;
close_figs = 0;
if close_figs 
    close all
end
if ~use_solution_as_init
    clearvars -except use_solution_as_init close_figs
end

%% Parameter
x0 = [0 1 0].'; l0 = [0 0 0].'; %l0 = 0.1*randn(4,1);
jlim = 3; use_umax = 0;
umax = jlim; umin = -jlim;
t0 = 0; t1 = 1; tf = 2; N = 100; fa = 1; fy = 1; fj = 1; kapparef_straight = 0.0; kapparef_curve = 1/(75.8778/2); s1 = 200; sf = s1+2*pi/kapparef_curve; %s1 = sf-2*pi/4*1/kapparef_curve; % Strecke, nach der von Gerade auf Kreis umgeschaltet wird
%sf = 400; %kapparef_curve = 0.01;
p.use_umax = use_umax; p.umax = umax; p.umin = umin; p.fa = fa; p.fy = fy; p.fj = fj; p.kapparef_straight = kapparef_straight; p.kapparef_curve = kapparef_curve; 
p.sf = sf; p.s1 = s1; p.t1 = t1;
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
start_inits = [-0.1 12 13];
inits = start_inits;
init_order = floor((log10(start_inits)>0).*log10(start_inits));
inits = 10.^(floor((log10(start_inits)>0).*log10(start_inits))).*abs(randn(size(start_inits))).*sign(start_inits);
error_flag = 1;
while error_flag
    try
        if use_solution_as_init
            solinit = bvpinit(sol,[p.t0 p.tf]);
        else
            solinit = bvpinit(t,init_guess,inits); % [nu_tilde, delta_t1, delta_t2]
            end
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
axopt = sol.y(3,:);
l1opt = sol.y(4,:);
l2opt = sol.y(5,:);
l3opt = sol.y(6,:);
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
if exist('u','var')
    if ~isempty(u)
        clear u
    end
end
for i=1:length(sol_mesh)
    u(:,i) = uopt(sol.y(:,i),p); % Steuerung
end
jopt = u;
kappaopt = [ones(size(sol_mesh_1))*p.kapparef_straight ones(size(sol_mesh_2))*p.kapparef_curve];
J_fun = 1/2*p.fj*jopt.^2+1/2*p.fa*axopt.^2+1/2*p.fy*kappaopt.^2.*vopt.^4;
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
subplot(4,1,1)
plot(sol_mesh,sopt)
ylabel('s_r [m]')
grid on
hold on
subplot(4,1,2)
plot(sol_mesh,vopt)
ylabel('v [m/s]')
grid on
hold on
subplot(4,1,3)
plot(sol_mesh,axopt)
ylabel('a_x_{opt} [m/s^2]')
xlabel('t [s]')
grid on
hold on
subplot(4,1,4)
plot(sol_mesh,jopt)
ylabel('j_x_{opt} [m/s^2]')
xlabel('t [s]')
grid on
hold on

figure(2)
subplot(3,1,1)
plot(sol_mesh,l1opt)
ylabel('l_{1,opt}')
grid on
hold on
subplot(3,1,2)
plot(sol_mesh,l2opt)
ylabel('l_{2,opt}')
grid on
hold on
subplot(3,1,3)
plot(sol_mesh,l3opt)
ylabel('l_{3,opt}')
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

figure(8)
plot(sopt,vopt)
hold on
grid on
ylabel('v [m/s]')
xlabel('s [m]')




