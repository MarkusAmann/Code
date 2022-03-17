% clc
clear all
close all

%% Parameter
x0 = [0 1 0].'; l0 = [0 0 0].'; %l0 = 0.1*randn(3,1);
jlim = 0.5; 
umax = jlim; umin = -jlim; use_umax = 0;
t0 = 0; tf = 30; N = 100; fa = 10; fj = 10; sf = 100; t1 = tf/3; s1 = 20;
tf_free = 1; % ACHTUNG: ich habe den Eindruck, dass die Optimierung bei aktiver Stellgrößenbeschränkung und freier Endzeit Probleme hat
% aufgrund des instabilen Teils der Lösung
p.use_umax = use_umax; p.umax = umax; p.umin = umin; p.fa = fa; p.fj = fj; p.sf = sf; p.s1 = s1;
p.x0 = x0; p.l0 = l0; p.t0 = t0; p.tf = tf; p.tf_free = tf_free; p.N = N; p.t1 = t1;


%% Optimierung
bvpoptions = bvpset('RelTol',1e-5,'Stats','on');
switch tf_free
    case 0
        t0_1 = linspace(p.t0, p.t1, p.N/(tf/t1));
        t1_f = linspace(p.t1, p.tf, p.N*(1-(t1/tf)));
        t = [t0_1 t1_f];
        deltat = mean(diff(t));
        p.deltat = deltat;
        init_guess = @(x,region)guess_fix_tf(t,p);
        solinit = bvpinit(t,init_guess,[0.1 1]); % [nu_tilde, delta_t1]
        sol = bvp4c(@sys_gesamt_fix_tf, @bcfcn_fix_tf, solinit, bvpoptions, p);

        % optimal states
        nu_tilde = sol.parameters(1);
        delta_t1_opt = sol.parameters(2)
        t1_opt = delta_t1_opt*p.t1
        sol_mesh = sol.x;
        split_idx = [find(diff(sol_mesh)==0) find(diff(sol_mesh)==0)+1];
        sol_mesh_t1 = delta_t1_opt*sol_mesh(1:split_idx(1));
        sol_mesh_t2 = sol_mesh(split_idx(2):end);
        sol_mesh_t2_zero_shift = sol_mesh_t2 - sol_mesh_t2(1);
        delta_2 = (tf-sol_mesh_t1(end))/sol_mesh_t2_zero_shift(end);
        sol_mesh_t2_zero_shift_transformed = delta_2*sol_mesh_t2_zero_shift;
        sol_mesh_t2 = sol_mesh_t2_zero_shift_transformed + sol_mesh_t1(end);
        t1_for_interp = linspace(p.t0, t1_opt, 1000);
        t2_for_interp = linspace(sol_mesh_t2(1), sol_mesh_t2(end), 1000);
        y_t1_interp = interp1(sol_mesh_t1.',sol.y(:,1:split_idx(1)).',t1_for_interp,'spline');
        y_t2_interp = interp1(sol_mesh_t2.',sol.y(:,split_idx(2):end).',t2_for_interp,'spline');
        sol_mesh = [t1_for_interp t2_for_interp];
        sol.x = sol_mesh;
        sol.y = [y_t1_interp; y_t2_interp].';

        sopt = sol.y(1,:);
        vopt = sol.y(2,:);
        aopt = sol.y(3,:);
        l1opt = sol.y(4,:);
        l2opt = sol.y(5,:);
        l3opt = sol.y(6,:);
        
        % optimal control inputs
        for i=1:length(sol_mesh)
           jopt(:,i) = uopt(sol.y(:,i),p); % Steuerung
        end
        
        J_fun = 1/2*p.fa*aopt.^2 + 1/2*p.fj*jopt.^2;
        J = trapz(sol_mesh,J_fun)

    case 1
        t0_1 = linspace(p.t0, p.t1, p.N/(tf/t1));
        t1_f = linspace(p.t1, p.tf, p.N*(1-(t1/tf)));
        t = [t0_1 t1_f];
        deltat = mean(diff(t));
        p.deltat = deltat;
        init_guess = @(x,region)guess_free_tf(t,p);
        solinit = bvpinit(t,init_guess,[0.1 1]); % [nu_tilde, delta_t1, delta_t2]
        sol = bvp4c(@sys_gesamt_free_tf, @bcfcn_free_tf, solinit, bvpoptions, p);
        % optimal states
        nu_tilde = sol.parameters(1);
        delta_tf_opt = sol.parameters(2)
        tf_opt = delta_tf_opt*p.tf
        sol_mesh = sol.x;
        sopt = sol.y(1,:);
        vopt = sol.y(2,:);
        aopt = sol.y(3,:);
        l1opt = sol.y(4,:);
        l2opt = sol.y(5,:);
        l3opt = sol.y(6,:);
        sol_mesh = sol_mesh*delta_tf_opt;
        
        % optimal control inputs
        for i=1:length(sol_mesh)
           jopt(:,i) = uopt(sol.y(:,i),p); % Steuerung
        end

        J_fun = 1/2*p.fa*aopt.^2 + 1/2*p.fj*jopt.^2;
        J = trapz(sol_mesh,J_fun) + tf_opt
end 

%%
figure('Name','sva')
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
plot(sol_mesh,aopt)
ylabel('a_x_{opt} [m/s^2]')
xlabel('t [s]')
grid on
hold on

figure('Name','adj')
subplot(3,1,1)
plot(sol_mesh,l1opt,'*')
ylabel('l_{1,opt}')
grid on
hold on
subplot(3,1,2)
plot(sol_mesh,l2opt,'*')
ylabel('l_{2,opt}')
grid on
hold on
subplot(3,1,3)
plot(sol_mesh,l3opt,'*')
ylabel('l_{3,opt}')
xlabel('t [s]')
grid on
hold on

figure('Name','j')
plot(sol_mesh,jopt)
ylabel('j_{opt} [m/s^3]')
xlabel('t [s]')
grid on
hold on

