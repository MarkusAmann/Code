syms r1 r2 r3 r4 r5 r6 r7 r8 r9
syms k fx fy 
a = fy*k^2;
A = [0 0 -1/fx; 0 0 0; -sqrt(12*a) -1 0];
R = [r1 r2 r3; r4 r5 r6; r7 r8 r9];
eqn = A.'*R + R*A;
Q = eye(3);
S = solve(eqn == -Q, [r1 r2 r3 r4 r5 r6 r7 r8 r9]);

%% Simulation of reduced system
fx = 1; fy = 1; kapparef = 0.005;
p.fx = fx; p.fy = fy; p.kapparef = kapparef; 

X_RL = [(4/(3*p.fy*p.kapparef^2))^(1/4);...
    -p.fy*p.kapparef^2*(4/(3*p.fy*p.kapparef^2))^(3/4);...
    0];
X0 = X_RL;
% X0 = 1.05*X_RL;
[t,X] = ode45(@(t,X) myodefun(t,X,p), 0:0.01:100, X0);
v_traj = X(:,1); l1_traj = X(:,2); l2_traj = X(:,3); 
ax_traj = -l2_traj/p.fx;
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

