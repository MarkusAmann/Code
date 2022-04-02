% clc
Doppelintegrator_fallunterscheidung_l1
clear all
% close all

%% Parameter
x0 = [0 10].'; l0 = [-0.5 0].'; %l0 = 0.1*randn(3,1);
alim = 0.2; 
umax = alim; umin = -alim; use_umax = 0;
t0 = 0; tf = 10; N = 1000; fa = 1; sf = 100; t1 = tf/8; s1 = 20;
tf_free = 1; % ACHTUNG: ich habe den Eindruck, dass die Optimierung bei aktiver Stellgrößenbeschränkung und freier Endzeit Probleme hat
% aufgrund des instabilen Teils der Lösung
p.use_umax = use_umax; p.umax = umax; p.umin = umin; p.fa = fa; p.sf = sf; p.s1 = s1;
p.x0 = x0; p.l0 = l0; p.t0 = t0; p.tf = tf; p.tf_free = tf_free; p.N = N; p.t1 = t1;

%% Optimierung
bvpoptions = bvpset('RelTol',1e-5,'Stats','on','NMax',4000);
switch tf_free
    case 0
        t = linspace(p.t0, p.tf, p.N+1);
        deltat = mean(diff(t));
        p.deltat = deltat;
        init_guess = @(x)guess_fix_tf(t,p);
        solinit = bvpinit(t,init_guess);
        sol = bvp4c(@sys_gesamt_fix_tf, @bcfcn_fix_tf, solinit, bvpoptions, p);
        % optimal states
        sol_mesh = sol.x;
        sopt = sol.y(1,:);
        vopt = sol.y(2,:);
        l1opt = sol.y(3,:);
        l2opt = sol.y(4,:);
        
    case 1
        t0_1 = linspace(p.t0, p.t1, p.N/(tf/t1));
        t1_f = linspace(p.t1, p.tf, p.N*(1-(t1/tf)));
        t = [t0_1 t1_f];
        deltat = mean(diff(t));
        p.deltat = deltat;
        init_guess = @(x,region)guess_free_tf(x,region,p);
        solinit = bvpinit(t,init_guess,[-0.1 0.8 0.8]); % [nu_tilde, delta_t1, delta_t2]
        sol = bvp4c(@sys_gesamt_free_tf, @bcfcn_free_tf, solinit, bvpoptions, p);
        % optimal states
        sol_mesh = sol.x;
        sopt = sol.y(1,:);
        vopt = sol.y(2,:);
        l1opt = sol.y(3,:);
        l2opt = sol.y(4,:);
        nu_tilde = sol.parameters(1);
        delta_t1_opt = sol.parameters(2)
        delta_t2_opt = sol.parameters(3)
        t1_opt = delta_t1_opt*p.t1
        tf_opt = delta_t1_opt*p.t1 + delta_t2_opt*(p.tf - p.t1)
        split_idx = [find(diff(sol_mesh)==0) find(diff(sol_mesh)==0)+1];
        sol_mesh_1 = sol_mesh(1:split_idx(1))*delta_t1_opt;
        sol_mesh_2 = delta_t1_opt*p.t1 + delta_t2_opt*(sol_mesh(split_idx(2):end) - p.t1);
        sol_mesh = [sol_mesh_1 sol_mesh_2];
        
end

% optimal control inputs
for i=1:length(sol_mesh)
   aopt(:,i) = uopt(sol.y(:,i),p); % Steuerung
end

J_fun = 1/2*p.fa*aopt.^2;
J = trapz(sol_mesh,J_fun) + tf_opt + t1_opt

%%
figure('Name','sva')
subplot(3,1,1)
plot(sol_mesh,sopt,':','LineWidth',2)
ylabel('s_r [m]')
grid on
hold on
subplot(3,1,2)
plot(sol_mesh,vopt,':','LineWidth',2)
ylabel('v [m/s]')
grid on
hold on
subplot(3,1,3)
plot(sol_mesh,aopt,':','LineWidth',2)
ylabel('a_x_{opt} [m/s^2]')
xlabel('t [s]')
grid on
hold on

figure('Name','adj')
subplot(2,1,1)
plot(sol_mesh,l1opt,':','LineWidth',2)
ylabel('l_{1,opt}')
grid on
hold on
subplot(2,1,2)
plot(sol_mesh,l2opt,':','LineWidth',2)
ylabel('l_{2,opt}')
grid on
hold on
