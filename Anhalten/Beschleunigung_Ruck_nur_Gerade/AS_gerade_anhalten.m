% close all; 
clear;

syms c1 c2 cv cs k1 k2;
syms t1;

%% hier parameter vorgeben
fa = 1;            % Gewichtung Beschl.
fj = 1;            % Gewichtung Ruck

s1 = 200;   % Länge der Gerade vor der Ampel
s2 = 500;    % Länge der Gerade nach der Ampel
sf = s1 + s2; % Gesamtlänge
t1 = 40;     % Zeit bis zum Erreichen der Ampel
kappa_ref = 0; % keine Kurve

s0 = 0; v0 = 10; vf = 0; a0 = 0; a1 = 0;

% J = (integral 0 bis t1) fj/2*j^2 + fa/2*ax^2 dt + t1 + (sf-s1)/v(t1)

%% randwertproblem

lambda1 = c1;
lambda2 = @(t) c2 - c1*t;
      j = @(t) k1*sqrt(fa/fj)*exp(sqrt(fa/fj)*t) - k2*sqrt(fa/fj)*exp(-sqrt(fa/fj)*t) + c1/fa;
      a = @(t) k1*exp(sqrt(fa/fj)*t) + k2*exp(-sqrt(fa/fj)*t) + c1/fa*t - c2/fa;
      v = @(t) k1*sqrt(fj/fa)*exp(sqrt(fa/fj)*t) - k2*sqrt(fj/fa)*exp(-sqrt(fa/fj)*t) + c1/(2*fa)*t^2 - c2/fa*t + cv;
      s = @(t) k1*(fj/fa)*exp(sqrt(fa/fj)*t) + k2*(fj/fa)*exp(-sqrt(fa/fj)*t) + c1/(6*fa)*t^3 - c2/(2*fa)*t^2 + cv*t + cs;

dd = vpasolve([s(0) == s0;...
    v(0) == v0;...
    a(0) == a0;...
    s(t1) == s1;...
%     a(t1) == a1;...
    j(t1) == 0;...
%     v(t1) == v1;...
    lambda2(t1) == -(sf-s1)/v(t1)^2], [c1 c2 cv cs k1 k2], [0,0,0,0,1e-40,0]);

% dd = vpasolve([s(0) == s0;...
%     v(0) == v0;...
%     a(0) == a0;...
%     s(t1) == s1;...
%     lambda2(t1) == -(sf-s1)/v(t1)^2;...
%     a(t1) == 0;...
%     1/2*fa*a(t1)^2-1/2*fj*j(t1)^2+lambda1*v(t1)+lambda2(t1)*a(t1)+1==0], [c1 c2 cv cs k1 k2 t1], [0,0,0,0,1e-40,0,15]);


c1_sol = double(dd.c1);
c2_sol = double(dd.c2);
cs_sol = double(dd.cs);
cv_sol = double(dd.cv);
k1_sol = double(dd.k1);
k2_sol = double(dd.k2);
if isfield(dd, 't1')
    t1_sol = double(dd.t1);
else
    t1_sol = t1;
end

lambda1_sol = matlabFunction(subs(lambda1, c1, c1_sol));
lambda2_sol = matlabFunction(subs(lambda2,{c1,c2},{c1_sol,c2_sol}));
      j_sol = matlabFunction(subs(j,{c1,k1,k2},{c1_sol,k1_sol,k2_sol}));
      a_sol = matlabFunction(subs(a,{c1,c2,k1,k2},{c1_sol,c2_sol,k1_sol,k2_sol}));
      v_sol = matlabFunction(subs(v,{c1,c2,k1,k2,cv},{c1_sol,c2_sol,k1_sol,k2_sol,cv_sol}));
      s_sol = matlabFunction(subs(s,{c1,c2,k1,k2,cv,cs},{c1_sol,c2_sol,k1_sol,k2_sol,cv_sol,cs_sol}));

v1_sol = double(v_sol(t1_sol));
tf_sol = (sf-s1)./v1_sol + t1_sol;

pos = find((imag(tf_sol)==0)&(real(tf_sol)>0));
if length(pos)>1
    pos = find(tf_sol > 0);
end

%% Wiederholung mit korrekter Lösung
c1_sol = double(dd.c1(pos));
c2_sol = double(dd.c2(pos));
cs_sol = double(dd.cs(pos));
cv_sol = double(dd.cv(pos));
k1_sol = double(dd.k1(pos));
k2_sol = double(dd.k2(pos));
if isfield(dd, 't1')
    t1_sol = double(dd.t1(pos));
else
    t1_sol = t1;
end

lambda1_sol = matlabFunction(subs(lambda1, c1, c1_sol));
lambda2_sol = matlabFunction(subs(lambda2,{c1,c2},{c1_sol,c2_sol}));
      j_sol = matlabFunction(subs(j,{c1,k1,k2},{c1_sol,k1_sol,k2_sol}));
      a_sol = matlabFunction(subs(a,{c1,c2,k1,k2},{c1_sol,c2_sol,k1_sol,k2_sol}));
      v_sol = matlabFunction(subs(v,{c1,c2,k1,k2,cv},{c1_sol,c2_sol,k1_sol,k2_sol,cv_sol}));
      s_sol = matlabFunction(subs(s,{c1,c2,k1,k2,cv,cs},{c1_sol,c2_sol,k1_sol,k2_sol,cv_sol,cs_sol}));

v1_sol = double(v_sol(t1_sol));
tf_sol = (sf-s1)./v1_sol + t1_sol;
l1_sol = double(lambda1_sol());

%% Lösungsvektoren
t1_vec = linspace(0,t1_sol,1000);
t2_vec = linspace(t1_sol,tf_sol,1000);
t_vec = [t1_vec t2_vec];

j_vec = [j_sol(t1_vec) zeros(size(t2_vec))];
a_vec = [a_sol(t1_vec) zeros(size(t2_vec))];
v_vec = [v_sol(t1_vec) v1_sol*ones(size(t2_vec))];
s_vec = [s_sol(t1_vec) linspace(s1, sf, length(t2_vec))];

l1_vec = [l1_sol*ones(size(t1_vec)) zeros(size(t2_vec))];
l2_vec = [lambda2_sol(t1_vec) zeros(size(t2_vec))];
l3_vec = -j_vec/fj;

kappa_vec = [zeros(size(t1_vec)) kappa_ref*ones(size(t2_vec))];
ay_vec = kappa_vec.*v_vec.^2;
dot_psi_vec = kappa_vec.*v_vec;
psi_vec = cumtrapz(t_vec,dot_psi_vec); 

% coordinate transformation of car movement to global coordinates
dx_global_opt_vec = v_vec.*cos(psi_vec);
dy_global_opt_vec = v_vec.*sin(psi_vec);
x_global_opt_vec = cumtrapz(t_vec,dx_global_opt_vec);
y_global_opt_vec = cumtrapz(t_vec,dy_global_opt_vec);

J_fun_1 = 1/2*fj*j_vec.^2 + 1/2*fa*a_vec.^2;
J_1 = trapz(t_vec,J_fun_1) + tf_sol;

J_fun_2 = 1/2*fj*j_sol(t1_vec).^2 + 1/2*fa*a_sol(t1_vec).^2;
J_2 = trapz(t1_vec,J_fun_2) + t1_sol + (sf-s1)/v1_sol;

fprintf('l1_l: %f\nt1: %f\ntf: %f\nJ: %f\n',c1_sol,t1_sol,tf_sol,J_2)
tmax = -sqrt(fj/fa)*log(l1_sol/k2_sol*sqrt(fj/fa^3));
amax = l1_sol*sqrt(fj/fa^3)-l1_sol*sqrt(fj/fa^3)*log(l1_sol/k2_sol*sqrt(fj/fa^3))-c2_sol/fa;

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
plot(t_vec,a_vec)
ylabel('a_x_{opt} [m/s^2]')
xlabel('t [s]')
grid on
hold on
subplot(4,1,4)
plot(t_vec,j_vec)
ylabel('j_x_{opt} [m/s^3]')
xlabel('t [s]')
grid on
hold on

figure(2)
subplot(3,1,1)
plot(t_vec,l1_vec)
ylabel('l_{1,opt}')
grid on
hold on
subplot(3,1,2)
plot(t_vec,l2_vec)
ylabel('l_{2,opt}')
grid on
hold on
subplot(3,1,3)
plot(t_vec,l3_vec)
ylabel('l_{3,opt}')
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





