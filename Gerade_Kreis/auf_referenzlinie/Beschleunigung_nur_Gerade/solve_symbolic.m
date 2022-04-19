clear all

%%
syms fx fy kapparef s0 s1 sf v0 % bekannte symbolische Variablen
syms t1 tf nu_1 nu_2 l1l a0 % unbekannte symbolische Variablen
known_vars = [fx fy kapparef s0 s1 sf v0];
known_vars_num = [1 1 0.1 0 40 100 5];
% known_vars_num(7) = ((2/(3*known_vars_num(2)*known_vars_num(3)^2))^(1/4));
unknwon_vars = [t1 tf nu_1 nu_2 l1l a0];

vr = (2/(3*fy*kapparef^2))^(1/4);
l1r = -2*fy*kapparef^2*vr^3;
ax_t1 = l1l*t1/fx + a0;
v_t1 = l1l*t1^2/(2*fx) + a0*t1 + v0;
l2l_t1 = -ax_t1*fx;

eqns = [l1l - l1r - 2*nu_1;...
    l1l*t1 + fx*a0 + 2*nu_2;...
    l1l*t1^3/(6*fx) + a0/2*t1^2 + v0*t1 + s0 - s1;...
    v_t1 - vr;...
    tf - t1 + (s1 - sf)/vr;...
    1/2*fx*ax_t1^2 + l1l*v_t1 + l2l_t1*ax_t1 - 1/2*fy*kapparef^2*vr^4 - l1r*vr];
eqns_num = subs(eqns,known_vars,known_vars_num);

S = solve(eqns_num, unknwon_vars);
t1 = vpa(S.t1,10);
tf = vpa(S.tf,10);
nu_1 = vpa(S.nu_1,10);
nu_2 = vpa(S.nu_2,10);
l1l = vpa(S.l1l,10);
a0 = vpa(S.a0,10);

pos = find((imag(t1)==0)&(real(t1)>0));
if length(pos)>1
    max   = inf;
    p_max = 0; 
    for i = 1:length(pos)
        if t1(pos(i)) < max
            max   = t1(pos(i));
            p_max = pos(i);
        end
    end
    pos = p_max;
end

t1_val = double(t1(pos));
tf_val = double(tf(pos));
nu_1_val = double(nu_1(pos));
nu_2_val = double(nu_2(pos));
l1l_val = double(l1l(pos));
a0_val = double(a0(pos));

t1_vec = linspace(0,t1_val,1000);
t2_vec = linspace(t1_val,tf_val,1000);
t_vec = [t1_vec t2_vec];
l1l_vec = l1l_val*ones(size(t1_vec));
l1r_vec = l1l_val*ones(size(t2_vec)) - 2*nu_1_val;
l1_vec = [l1l_vec l1r_vec];
l2l_vec = -l1l_vec.*t1_vec - known_vars_num(1)*a0_val;
l2r_vec = zeros(size(t2_vec));
l2_vec = [l2l_vec l2r_vec];
ax_vec = - l2_vec/known_vars_num(1);
v_vec = cumtrapz(t_vec,ax_vec) + known_vars_num(7);
s_vec = cumtrapz(t_vec,v_vec) + known_vars_num(4);
kappa_vec = [zeros(size(t1_vec)) known_vars_num(3)*ones(size(t2_vec))];
ay_vec = v_vec.^2.*kappa_vec;
dot_psi_vec = kappa_vec.*v_vec;
psi_vec = cumtrapz(t_vec,dot_psi_vec); 

% coordinate transformation of car movement to global coordinates
dx_global_opt_vec = v_vec.*cos(psi_vec);
dy_global_opt_vec = v_vec.*sin(psi_vec);
x_global_opt_vec = cumtrapz(t_vec,dx_global_opt_vec);
y_global_opt_vec = cumtrapz(t_vec,dy_global_opt_vec);

J_fun = 1/2*known_vars_num(1)*ax_vec.^2 + 1/2*known_vars_num(2)*kappa_vec.^2.*v_vec.^4;
J = trapz(t_vec,J_fun) + tf_val;

fprintf('l1_1: %f\nl1_2: %f\nnu_tilde_1: %f\nnu_tilde_2: %f\nt1: %f\ntf: %f\nJ: %f\n',l1_vec(1),l1_vec(end),nu_1_val,nu_2_val,t1_val,tf_val,J)

%% Ruhelagenwerte
v_RL = ((2/(3*known_vars_num(2)*known_vars_num(3)^2))^(1/4));
l1_RL = -2*known_vars_num(2)*known_vars_num(3)^2*v_RL^3;
ay_RL = v_RL^2*known_vars_num(3);

%% Visualization
figure
subplot(3,1,1)
% plot(t_vec,s_vec,'--','LineWidth',2)
plot(t_vec,s_vec)
ylabel('s [m]')
grid on
hold on
subplot(3,1,2)
% plot(t_vec,v_vec,'--','LineWidth',2)
plot(t_vec,v_vec)
ylabel('v [m/s]')
grid on
hold on
subplot(3,1,3)
% plot(t_vec,ax_vec,'--','LineWidth',2)
plot(t_vec,ax_vec)
ylabel('a [m/s^2]')
xlabel('t [s]')
grid on
hold on

figure('Name','adj')
subplot(2,1,1)
% plot(t_vec,l1_vec,'--','LineWidth',2)
plot(t_vec,l1_vec)
ylabel('l_{1,opt}')
grid on
hold on
subplot(2,1,2)
% plot(t_vec,l2_vec,'--','LineWidth',2)
plot(t_vec,l2_vec)
ylabel('l_{2,opt}')
grid on
hold on

figure('Name','kapparef')
plot(t_vec,kappa_vec)
ylabel('\kappa_{opt} [1/m]')
xlabel('t [s]')
grid on
hold on

figure('Name','ay')
plot(t_vec,ay_vec)
ylabel('a_y [m/s^2]')
xlabel('t [s]')
grid on
hold on

figure('Name','psi')
plot(t_vec,psi_vec)
ylabel('\psi_{opt} [rad]')
xlabel('t [s]')
grid on
hold on

figure('Name','xy-Position')
plot(x_global_opt_vec,y_global_opt_vec)
grid on
hold on
ylabel('y position [m]')
xlabel('x position [m]')
legend('trajectory')

figure('Name','xy-Time')
subplot(2,1,1)
plot(t_vec, x_global_opt_vec)
hold on
grid on
ylabel('x position [m]')
subplot(2,1,2)
plot(t_vec, y_global_opt_vec)
hold on
grid on
ylabel('y position [m]')
xlabel('t [s]')

figure('Name', 'v-s')
plot(s_vec,v_vec)
hold on
grid on
ylabel('v [m/s]')
xlabel('s [m]')
