% clc
clear all
close all

%% Parameter
x0 = [0 5 0 0].'; l0 = 0*[0.1 0.1 0.1 0].'; %l0 = 0.1*randn(4,1);
alim = 1.06*1000; kappalim = 1/4*1000;
umax = [alim;kappalim]; umin = -[alim;kappalim]; use_umax = 0;
t0 = 0; tf = 1; N = 100; fx = 1; fy = 1; fr = 1; fp = 1; kapparef = 0.001; sf = 100; drf = 0; psirf = 0;
tf_free = 1;
p.use_umax = use_umax; p.umax = umax; p.umin = umin; p.fx = fx; p.fy = fy; p.fr = fr; p.fp = fp; p.kapparef = kapparef; p.sf = sf; p.drf = drf; p.psirf = psirf;
p.x0 = x0; p.l0 = l0; p.t0 = t0; p.tf = tf; p.tf_free = tf_free; p.N = N;  

if kapparef == 0.01
    X_init = [0;5;0;0;-0.1067;-0.6478;-0.0876;-5.0446]; tf = 11.5188; %Anfangswert für kapparef = 0.01; sf = 100; ohne Bestrafung von Zuständen, vf, drf, psirf sind frei
elseif kapparef == 0
    X_init = [0;5;0;0;-0.0923;-1.0379;0;0]; tf = 11.2472; %Anfangswert für kapparef = 0; sf = 100; ohne Bestrafung von Zuständen, vf, drf, psirf sind frei
end
X_init = [0;5;0;0;-0.0923;-1.0379;0;0]; tf = 11.2472; % Anfangswert von optimaler Geradeausfahrt
l4_init = -1.25*0.5;
X_init = [0;5;0;pi/2;kapparef*l4_init;0;-l4_init^2/(fy*5^3);l4_init]; tf = 11.2472*10; 

[t,X] = ode45(@(t,X) odefun_kreisfahrt(t,X,p), 0:0.01:tf, X_init);
X=X.';
sr_traj = X(1,:); v_traj = X(2,:); dr_traj = X(3,:); psir_traj = X(4,:);...
    l1_traj = X(5,:); l2_traj = X(6,:); l3_traj = X(7,:); l4_traj = X(8,:);

%%
% optimal control inputs
for i=1:length(t)
    u(:,i) = uopt(X(:,i),p); % Steuerung
end
axopt = u(1,:);
kappaopt = u(2,:);
% kappaopt(1) = interp1(t(2:end),kappaopt(2:end),0,'pchip','extrap');

% time-derivatives of optimal values and values of reference curve
dot_sopt = v_traj.*cos(psir_traj)./(1-dr_traj*kapparef);
dot_psi = kappaopt.*v_traj;
dot_psiref = kapparef*dot_sopt;
psiref = cumtrapz(t,dot_psiref);
dot_psir = dot_psi - dot_psiref;
vref = dot_sopt;

% lateral acc., heading angle of car
ayopt = v_traj.^2.*kappaopt;
psiopt = cumtrapz(t,dot_psi); 

% coordinate transformation of car movement to global coordinates
dx_global_opt = v_traj.*cos(psiopt);
dy_global_opt = v_traj.*sin(psiopt);
x_global_opt = cumtrapz(t,dx_global_opt);
y_global_opt = cumtrapz(t,dy_global_opt);

% coordinate transformation of reference curve to global coordinates
dx_ref = vref.*cos(psiref);
dy_ref = vref.*sin(psiref);
x_ref = cumtrapz(t,dx_ref);
y_ref = cumtrapz(t,dy_ref);

%%
figure
subplot(3,1,1)
plot(t,sr_traj)
ylabel('s_r [m]')
grid on
hold on
subplot(3,1,2)
plot(t,v_traj)
ylabel('v [m/s]')
grid on
hold on
plot(t,gradient(sr_traj,t),'-.')
legend('vx','dot\_sr')
subplot(3,1,3)
plot(t,axopt)
ylabel('a_x [m/s^2]')
xlabel('t [s]')
grid on
hold on

figure
subplot(2,1,1)
plot(t, dr_traj)
hold on
grid on
ylabel('d_r [m]')
subplot(2,1,2)
plot(t, psir_traj)
hold on
grid on
ylabel('psi_r [rad]')
xlabel('t [s]')

figure
subplot(2,2,1)
plot(t,l1_traj)
ylabel('l_{1}')
grid on
hold on
subplot(2,2,2)
plot(t,l2_traj)
ylabel('l_{2}')
grid on
hold on
subplot(2,2,3)
plot(t,l3_traj)
ylabel('l_{3}')
xlabel('t [s]')
grid on
hold on
subplot(2,2,4)
plot(t,l4_traj)
ylabel('l_{4}')
xlabel('t [s]')
grid on
hold on

figure 
plot(t,kappaopt)
ylabel('\kappa [1/m]')
xlabel('t [s]')
grid on
hold on

figure 
plot(t,ayopt)
ylabel('a_y [m/s^2]')
xlabel('t [s]')
grid on
hold on

figure 
plot(t,psiopt)
ylabel('\psi [rad]')
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
plot(t, x_global_opt)
hold on
grid on
ylabel('x position [m]')
subplot(2,1,2)
plot(t, y_global_opt)
hold on
grid on
ylabel('y position [m]')
xlabel('t [s]')

function dXdt = odefun_kreisfahrt(t,X,p)
    sr = X(1); v = X(2); dr = X(3); psir = X(4); l1 = X(5); l2 = X(6); l3 = X(7); l4 = X(8); 
    u = uopt(X,p);
    ax = u(1);
    kappa = u(2);
    % ohne Zustandsbestrafung im Gütefunktional mit J=tf+int(1/2*fx*ax^2+1/2*fy*kappa^2*v^4)
    dXdt = [v*cos(psir)/(1-dr*p.kapparef);...
        ax;...
        v*sin(psir);...
        kappa*v - p.kapparef*v*cos(psir)/(1-dr*p.kapparef);...
        0;...
        -(2*p.fy*v^3*kappa^2 + l1*cos(psir)/(1-dr*p.kapparef) + l3*sin(psir) + l4*kappa - l4*p.kapparef*cos(psir)/(1-dr*p.kapparef));...
        -(p.kapparef*l1*v*cos(psir)/(1-dr*p.kapparef)^2 - p.kapparef^2*l4*v*cos(psir)/(1-dr*p.kapparef)^2);...
        -(l3*v*cos(psir) - l1*v*sin(psir)/(1-dr*p.kapparef) + p.kapparef*l4*v*sin(psir)/(1-dr*p.kapparef))];
end

