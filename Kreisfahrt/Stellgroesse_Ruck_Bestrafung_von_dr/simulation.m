% clc
clear all
close all

%% Parameter
x0 = [0 5 0 0 0].'; l0 = [0 0 0 0 0].'; %l0 = 0.1*randn(4,1);
jlim = 1.06*1000; kappalim = 1/4*1000; use_umax = 0; use_dr = 1;
umax = [jlim;kappalim]; umin = -[jlim;kappalim];
t0 = 0; tf = 1; N = 100; fj = 1; fa = 1; fy = 1; fr = 0.1; kapparef = 0.1; sf =100; drf = 0; psirf = 0;
tf_free = 1;
% l4_init = -0.125*0;
% x0 = [0;5;0;-pi/2]; l0=[kapparef*l4_init;0;-l4_init^2/(fy*5^3);l4_init];
p.use_umax = use_umax; p.umax = umax; p.umin = umin; p.use_dr = use_dr; p.fj = fj; p.fa = fa; p.fy = fy; p.fr = fr; p.kapparef = kapparef; p.sf = sf; p.drf = drf; p.psirf = psirf;
p.x0 = x0; p.l0 = l0; p.t0 = t0; p.tf = tf; p.tf_free = tf_free; p.N = N;  


%% Ruhelage der Kreisfahrt mir Bestrafung von dr und 1-dr*kappar = 1
syms vsym
v_RL_all = solve(1/2*p.fy^2/p.fr*p.kapparef^6*vsym^8-3/2*p.fy*kapparef^2*vsym^4+1 == 0, vsym);
v_RL_all = double(v_RL_all);
v_RL = min(v_RL_all(imag(v_RL_all) == 0 & real(v_RL_all) >= 0));
dr_RL = p.kapparef^3*v_RL^4*p.fy/p.fr;
l1_RL = -2*p.fy*p.kapparef^2*v_RL^3;
l5_RL = -p.fy*p.kapparef*v_RL^3;
kappa_RL = p.kapparef;
kappa_RL_corrected = p.kapparef/(1-dr_RL*p.kapparef);
ay_RL = kappa_RL*v_RL^2;
p.x0 = [0 v_RL 0 dr_RL 0].'; p.l0 = [l1_RL 0 0 0 l5_RL].';

%% Optimierung
bvpoptions = bvpset('RelTol',1e-5,'Stats','on');
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
        aopt = sol.y(3,:);
        dropt = sol.y(4,:);
        psiropt = sol.y(5,:);
        l1opt = sol.y(6,:);
        l2opt = sol.y(7,:);
        l3opt = sol.y(8,:);
        l4opt = sol.y(9,:);
        l5opt = sol.y(10,:);

    case 1
        t = linspace(p.t0, p.tf, p.N+1);
        deltat = mean(diff(t));
        p.deltat = deltat;
        init_guess = @(x)guess_free_tf(t,p);
        start_inits = [20];
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
        aopt = sol.y(3,:);
        dropt = sol.y(4,:);
        psiropt = sol.y(5,:);
        l1opt = sol.y(6,:);
        l2opt = sol.y(7,:);
        l3opt = sol.y(8,:);
        l4opt = sol.y(9,:);
        l5opt = sol.y(10,:);
        delta_opt = sol.parameters;
        tf_opt = delta_opt*p.tf;
        sol_mesh = sol_mesh*delta_opt;
end

%%
% optimal control inputs
for i=1:length(sol_mesh)
    u(:,i) = uopt(sol.y(:,i),p); % Steuerung
end
jopt = u(1,:);
kappaopt = u(2,:);
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
y_global_opt = cumtrapz(sol_mesh,dy_global_opt) + p.x0(4);

% coordinate transformation of reference curve to global coordinates
dx_ref = vref.*cos(psiref);
dy_ref = vref.*sin(psiref);
x_ref = cumtrapz(sol_mesh,dx_ref);
y_ref = cumtrapz(sol_mesh,dy_ref);

%%
figure
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
plot(sol_mesh,aopt)
ylabel('a_x_{opt} [m/s^2]')
xlabel('t [s]')
grid on
hold on
subplot(4,1,4)
plot(sol_mesh,jopt)
ylabel('j_x_{opt} [m/s^3]')
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





