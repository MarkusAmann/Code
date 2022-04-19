clear all;

syms fx fy kapparef s0 s1 sf v0 % bekannte symbolische Variablen
syms t1 l1 a0 % unbekannte symbolische Variablen
known_vars = [fx fy kapparef s0 s1 sf v0];
known_vars_num = [1 1 0.1 0 80 200 5];
% known_vars_num(7) = ((2/(3*known_vars_num(2)*known_vars_num(3)^2))^(1/4));
% known_vars_num(7) = 5;

unknwon_vars = [t1 l1 a0];

ax_t1 = l1/fx*t1 + a0;
v_t1 = l1*t1^2/(2*fx) + a0*t1 + v0;

eqns = [l1*t1^3/(6*fx) + a0/2*t1^2 + v0*t1 + s0 - s1;...
    ax_t1 - (sf - s1)/(fx*v_t1^2) + 3/2*(sf - s1)*fy/fx*kapparef^2*v_t1^2;...
    -1/2*fx*ax_t1^2 + l1*v_t1 + 1];
eqns_num = subs(eqns,known_vars,known_vars_num);

S = solve(eqns_num, unknwon_vars);
t1 = vpa(S.t1,10);
l1 = vpa(S.l1,10);
a0 = vpa(S.a0,10);
t1_val = double(t1);
l1_val = double(l1);
a0_val = double(a0);

pos = find((imag(t1_val)==0)&(real(t1_val)>0));
if length(pos)>1
    max   = inf;
    p_max = 0; 
    for i = 1:length(pos)
        if t1_val(pos(i)) < max
            max   = t1_val(pos(i));
            p_max = pos(i);
        end
    end
    pos = p_max;
end

fx = known_vars_num(1); fy = known_vars_num(2); kapparef = known_vars_num(3); s0 = known_vars_num(4); s1 = known_vars_num(5); sf = known_vars_num(6); v0 = known_vars_num(7);

a0_sol = double(a0_val(pos));
l1_sol = double(l1_val(pos));
t1_sol = double(t1_val(pos));
v1 = l1_sol*t1_sol.^2/(2*fx) + a0_sol*t1_sol + v0;
tf_sol = (sf - s1)/v1 + t1_sol;

%% Zeitlösungen ausrechnen
t1_vec = linspace(0,t1_sol,1000);
t2_vec = linspace(t1_sol,tf_sol,1000);
t_vec = [t1_vec t2_vec];

l2l_vec = -l1_sol*t1_vec - fx*a0_sol;
l2r_vec = zeros(size(t2_vec));
l2_vec = [l2l_vec l2r_vec];

axl_vec = -l2l_vec/fx;
ax_vec = -l2_vec/fx;

vl_vec = l1_sol*t1_vec.^2/(2*fx) + a0_sol*t1_vec + v0;
vr_vec = vl_vec(end)*ones(size(t2_vec));
v_vec = [vl_vec vr_vec];
v_vec_int = cumtrapz(t_vec,ax_vec) + v0;

sl_vec = l1_sol*t1_vec.^3/(6*fx) + a0_sol/2*t1_vec.^2 + v0.*t1_vec + s0;
sr_vec = vl_vec(end).*t2_vec + s1 - v1*t1_sol;
s_vec = [sl_vec sr_vec];
s_vec_int = cumtrapz(t_vec,v_vec) + s0;

kappa_vec = [zeros(size(t1_vec)) kapparef*ones(size(t2_vec))];
ay_vec = v_vec.^2.*kappa_vec;
dot_psi_vec = kappa_vec.*v_vec;
psi_vec = cumtrapz(t_vec,dot_psi_vec); 

% coordinate transformation of car movement to global coordinates
dx_global_opt_vec = v_vec.*cos(psi_vec);
dy_global_opt_vec = v_vec.*sin(psi_vec);
x_global_opt_vec = cumtrapz(t_vec,dx_global_opt_vec);
y_global_opt_vec = cumtrapz(t_vec,dy_global_opt_vec);

J_fun = 1/2*fx*ax_vec.^2 + 1/2*fy*kappa_vec.^2.*v_vec.^4;
J = trapz(t_vec,J_fun) + tf_sol;

J_alt_fun = 1/2*fx*axl_vec.^2;
J_alt = trapz(t1_vec,J_alt_fun) + t1_sol + (sf - s1)/(v1) + 1/2*(sf - s1)/(v1)*fy*kapparef^2*v1^4;

fprintf('l1_l: %f\nt1: %f\ntf: %f\nJ: %f\n',l1_sol,t1_sol,tf_sol,J)

%% Ruhelagenwerte
v_RL = ((2/(3*fy*kapparef^2))^(1/4));
l1_RL = -2*fy*kapparef^2*v_RL^3;
ay_RL = v_RL^2*kapparef;

%%
figure(1)
subplot(4,1,1)
plot(t_vec,s_vec)
ylabel('s_r [m]')
grid on
hold on
subplot(4,1,2)
plot(t_vec,v_vec)
ylabel('v [m/s]')
grid on
hold on
subplot(4,1,3)
plot(t_vec,ax_vec)
ylabel('a_x_{opt} [m/s^2]')
xlabel('t [s]')
grid on
hold on
subplot(4,1,4)
plot(t_vec,gradient(ax_vec,t_vec))
ylabel('j_x_{opt} [m/s^3]')
xlabel('t [s]')
grid on
hold on

figure(2)
subplot(2,1,1)
plot(t_vec,[l1_sol*ones(size(t1_vec)) zeros(size(t2_vec))])
ylabel('l_{1,opt}')
grid on
hold on
subplot(2,1,2)
plot(t_vec,l2_vec)
ylabel('l_{2,opt}')
grid on
hold on

figure(3) 
plot(t_vec,kappa_vec)
ylabel('\kappa_{opt} [1/m]')
xlabel('t [s]')
grid on
hold on

figure(4) 
plot(t_vec,ay_vec)
ylabel('a_y [m/s^2]')
xlabel('t [s]')
grid on
hold on

figure(5)
plot(t_vec,psi_vec)
ylabel('\psi_{opt} [rad]')
xlabel('t [s]')
grid on
hold on

figure(6)
plot(x_global_opt_vec,y_global_opt_vec)
grid on
hold on
ylabel('y position [m]')
xlabel('x position [m]')
legend('trajectory')

figure(7)
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

figure(8)
plot(s_vec,v_vec)
hold on
grid on
ylabel('v [m/s]')
xlabel('s [m]')




