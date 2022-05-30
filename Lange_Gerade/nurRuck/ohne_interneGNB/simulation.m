% clc
clear all
close all

%% Parameter
x0 = [0 0.1 0].'; l0 = [0 0 0].'; %l0 = 0.1*randn(3,1);
jlim = 0.5; 
umax = jlim; umin = -jlim; use_umax = 0;
t0 = 0; tf = 1; N = 100; fj = 1; sf = 1000; 
tf_free = 1; % ACHTUNG: ich habe den Eindruck, dass die Optimierung bei aktiver Stellgrößenbeschränkung und freier Endzeit Probleme hat
% aufgrund des instabilen Teils der Lösung
p.use_umax = use_umax; p.umax = umax; p.umin = umin; p.fj = fj; p.sf = sf; 
p.x0 = x0; p.l0 = l0; p.t0 = t0; p.tf = tf; p.tf_free = tf_free; p.N = N;  

%% Optimierung
bvpoptions = bvpset('RelTol',1e-5,'Stats','on');
t = linspace(p.t0, p.tf, p.N+1);
deltat = mean(diff(t));
p.deltat = deltat;
init_guess = @(x)guess_free_tf(t,p);
solinit = bvpinit(t,init_guess,0.1);
sol = bvp4c(@sys_gesamt_free_tf, @bcfcn_free_tf, solinit, bvpoptions, p);
% optimal states
sol_mesh = sol.x;
sopt = sol.y(1,:);
vopt = sol.y(2,:);
aopt = sol.y(3,:);
l1opt = sol.y(4,:);
l2opt = sol.y(5,:);
l3opt = sol.y(6,:);
delta_opt = sol.parameters
tf_opt = delta_opt*p.tf
sol_mesh = sol_mesh*delta_opt;

% optimal control inputs
for i=1:length(sol_mesh)
   jopt(:,i) = uopt(sol.y(:,i),p); % Steuerung
end

J_fun_2 = 1/2*fj*jopt.^2;
J_2 = trapz(sol_mesh,J_fun_2) + tf_opt*p.tf_free;
fprintf('l1: %f\ntf: %f\nJ: %f\n',l1opt(1),tf_opt,J_2)

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
xlabel('t [s]')
grid on
hold on

figure('Name','j')
plot(sol_mesh,jopt)
ylabel('j_{opt} [m/s^3]')
xlabel('t [s]')
grid on
hold on


%%
function F = fun(x,p)
% fa = p.fa; fj = p.fj;
k1 = x(1); k2 = x(2); tf = x(3); l1 = x(4); c2  = x(5); cv  = x(6); cs  = x(7);
x0 = p.x0; s0 = x0(1); v0 = x0(2); a0 = x0(3); sf = p.sf;
F(1) = cs + k1 + k2 - s0;
F(2) = cv + k1 - k2 - v0;
F(3) = k1 - c2 + k2 - a0;
F(4) = cs + k2*exp(-tf) + cv*tf - (c2*tf^2)/2 + (l1*tf^3)/6 + k1*exp(tf) - sf;
F(5) = c2 - l1*tf;
F(6) = k2*exp(-tf) - l1 - k1*exp(tf);
F(7) = l1*(cv - k2*exp(-tf) - c2*tf + (l1*tf^2)/2 + k1*exp(tf)) + (k2*exp(-tf) - c2 + l1*tf + k1*exp(tf))^2/2 + 1;
end