% close all; 
clear;

syms c1 c2;
syms tf;

%% hier parameter vorgeben
fa = 1;            % Gewichtung Längsbeschleunigung

sf = 100;   % Länge der Gerade  

s0 = 0; v0 = 10; vf = 0;

%% randwertproblem

lambda1 = c1;
      a = @(t) c1/fa*t - c2/fa;
      v = @(t) c1/(2*fa)*t^2 - c2/fa*t + v0;
      s = @(t) c1/(6*fa)*t^3 - c2/(2*fa)*t^2 + v0*t + s0;

tf = 25;
dd = vpasolve([s(tf) == sf;...
%     a(tf) == 0;...
    v(tf) == vf;...
    ], [0,0]);

% dd = vpasolve([s(tf) == sf;...
%     a(tf) == 0;...
% %     v(tf) == vf;...
%     -1/2*fa*a(tf)^2+lambda1*v(tf)+1==0], [0,0,10]);

c1_sol = double(dd.c1);
c2_sol = double(dd.c2);
if isfield(dd, 'tf')
    tf_sol = double(dd.tf);
else
    tf_sol = tf;
end

lambda1_sol = matlabFunction(subs(lambda1, c1, c1_sol));
      a_sol = matlabFunction(subs(a,{c1,c2},{c1_sol,c2_sol}));
      v_sol = matlabFunction(subs(v,{c1,c2},{c1_sol,c2_sol}));
      s_sol = matlabFunction(subs(s,{c1,c2},{c1_sol,c2_sol}));

vf_sol = double(v_sol(tf_sol));

pos = find((imag(tf_sol)==0)&(real(tf_sol)>0));
if length(pos)>1
    pos = find(tf_sol > 0);
end

%% Wiederholung mit korrekter Lösung
c1_sol = double(dd.c1(pos));
c2_sol = double(dd.c2(pos));
if isfield(dd, 'tf')
    tf_sol = double(dd.tf(pos));
else
    tf_sol = tf;
end

lambda1_sol = matlabFunction(subs(lambda1, c1, c1_sol));
      a_sol = matlabFunction(subs(a,{c1,c2},{c1_sol,c2_sol}));
      v_sol = matlabFunction(subs(v,{c1,c2},{c1_sol,c2_sol}));
      s_sol = matlabFunction(subs(s,{c1,c2},{c1_sol,c2_sol}));

vf_sol = double(v_sol(tf_sol));
l1_sol = double(lambda1_sol());

%% Lösungsvektoren
t_vec = linspace(0,tf_sol,1000);

a_vec = a_sol(t_vec);
v_vec = v_sol(t_vec);
s_vec = s_sol(t_vec);
j_vec = gradient(a_vec,t_vec);

l1_vec = l1_sol*ones(size(t_vec));
l2_vec = -a_vec/fa;

kappa_vec = zeros(size(t_vec));
ay_vec = kappa_vec.*v_vec.^2;
dot_psi_vec = kappa_vec.*v_vec;
psi_vec = cumtrapz(t_vec,dot_psi_vec); 

% coordinate transformation of car movement to global coordinates
dx_global_opt_vec = v_vec.*cos(psi_vec);
dy_global_opt_vec = v_vec.*sin(psi_vec);
x_global_opt_vec = cumtrapz(t_vec,dx_global_opt_vec);
y_global_opt_vec = cumtrapz(t_vec,dy_global_opt_vec);

J_fun_1 = 1/2*fa*a_vec.^2;
J_1 = trapz(t_vec,J_fun_1) + tf_sol;

fprintf('l1_l: %f\ntf: %f\nJ: %f\n',c1_sol,tf_sol,J_1)

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
subplot(2,1,1)
plot(t_vec,l1_vec)
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





