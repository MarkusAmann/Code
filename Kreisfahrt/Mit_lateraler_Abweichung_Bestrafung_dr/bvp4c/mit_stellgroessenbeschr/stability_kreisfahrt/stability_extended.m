close all
clear all
% clc

%% Simulation of extended system
fj = 1; fy = 1; fr = 1; kapparef = 0.01;
p.fj = fj; p.fy = fy; p.fr = fr; p.kapparef = kapparef; 
ruhelage_extended;

X0 = X_RL;
% X0 = 1.05*X_RL;

X_init = [X0(1); 0; X0(2); 0; X0(3); 0; 0; 0; X0(4)];
[t,X] = ode45(@(t,X) myodefun_extended(t,X,p), 0:0.001:100, X_init);
v_traj = X(:,1); ax_traj = X(:,2); dr_traj = X(:,3); psir_traj = X(:,4);...
    l1_traj = X(:,5); l2_traj = X(:,6); l3_traj = X(:,7); l4_traj = X(:,8); l5_traj = X(:,9); 
j_traj = -l3_traj/p.fj;
s_traj = cumtrapz(t,v_traj);
dot_psi_traj = p.kapparef*v_traj;
psi_traj = cumtrapz(t,dot_psi_traj);

% coordinate transformation of car movement to global coordinates
dx_global = v_traj.*cos(psi_traj);
dy_global = v_traj.*sin(psi_traj);
x_global = cumtrapz(t,dx_global);
y_global = cumtrapz(t,dy_global);

%%
figure
subplot(3,1,1)
plot(t,s_traj)
ylabel('s [m]')
grid on
hold on
subplot(3,1,2)
plot(t,v_traj)
ylabel('v [m/s]')
grid on
hold on
subplot(3,1,3)
plot(t,ax_traj)
ylabel('a_x [m/s^2]')
xlabel('t [s]')
grid on
hold on

figure
plot(x_global,y_global)
grid on
hold on
ylabel('y position [m]')
xlabel('x position [m]')

figure
subplot(2,1,1)
plot(t,l1_traj)
ylabel('l_1')
grid on
hold on
subplot(2,1,2)
plot(t,l2_traj)
ylabel('l_2')
grid on
hold on

