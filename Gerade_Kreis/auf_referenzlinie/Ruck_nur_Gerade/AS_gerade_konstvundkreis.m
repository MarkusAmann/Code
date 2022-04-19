% close all; 
clear;

syms c1 c2 c3;
syms t1;

%% hier parameter vorgeben
fj = 1;            % Gewichtung Ruck.
fy = 1;            % Gewichtung Querbeschl.
kappa_ref = 0.1;    % Krümmung

s1 = 100;   % Länge der Gerade
s2 = 100;    % Länge der Kurve
sf = s1 + s2; % Gesamtlänge

s0 = 0; v0 = 1; a0 = 0; a1 = 0;

% J = (integral 0 bis t1) fj/2*j^2 dt + t1 + (sf-s1)/v(t1) + 1/2*(sf-s1)/v(t1)*(fy*v(t1)^4*kappa_ref^2)

%% randwertproblem

lambda1 = c1;
lambda2 = @(t) c2 - c1*t;
      j = @(t) -c1/(2*fj)*t^2 + c2/fj*t - c3/fj;
      a = @(t) -c1/(6*fj)*t^3 + c2/(2*fj)*t^2 - c3/fj*t + a0;
      v = @(t) -c1/(24*fj)*t^4 + c2/(6*fj)*t^3 - c3/(2*fj)*t^2 + a0*t + v0;
      s = @(t) -c1/(120*fj)*t^5 + c2/(24*fj)*t^4 - c3/(6*fj)*t^3 + a0/2*t^2 + v0*t + s0;

dd = vpasolve([s(t1) == s1;...
               a(t1) == a1;...
               lambda2(t1) == (sf-s1)*3/2*fy*v(t1)^2*kappa_ref^2 - (sf-s1)/v(t1)^2;...
               -1/2*fj*j(t1)^2+lambda1*v(t1)+lambda2(t1)*a(t1)+1==0], [0,0,0,15]);

pos = find((imag(dd.t1)==0)&(real(dd.t1)>0));
if length(pos)>1
    max   = inf;
    p_max = 0; 
    for i = 1:length(pos)
        if dd.t1(pos(i)) < max
            max   = dd.t1(pos(i));
            p_max = pos(i);
        end
    end
    pos = p_max;
end

c1_sol = double(dd.c1(pos));
c2_sol = double(dd.c2(pos));
c3_sol = double(dd.c3(pos));
t1_sol = double(dd.t1(pos));

lambda1_sol = matlabFunction(subs(lambda1, c1, c1_sol));
lambda2_sol = matlabFunction(subs(lambda2,{c1,c2},{c1_sol,c2_sol}));
      j_sol = matlabFunction(subs(j,{c1,c2,c3},{c1_sol,c2_sol,c3_sol}));
      a_sol = matlabFunction(subs(a,{c1,c2,c3},{c1_sol,c2_sol,c3_sol}));
      v_sol = matlabFunction(subs(v,{c1,c2,c3},{c1_sol,c2_sol,c3_sol}));
      s_sol = matlabFunction(subs(s,{c1,c2,c3},{c1_sol,c2_sol,c3_sol}));

v1_sol = double(v_sol(t1_sol));
tf_sol = (sf-s1)/v1_sol + t1_sol;
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

J_fun_1 = 1/2*fj*j_vec.^2 + 1/2*fy*kappa_vec.^2.*v_vec.^4;
J_1 = trapz(t_vec,J_fun_1) + tf_sol;

J_fun_2 = 1/2*fj*j_sol(t1_vec).^2;
J_2 = trapz(t1_vec,J_fun_2) + t1_sol + (sf-s1)/v1_sol + 1/2*(sf-s1)/v1_sol*(fy*v1_sol^4*kappa_ref^2);

fprintf('l1_l: %f\nt1: %f\ntf: %f\nJ: %f\n',c1_sol,t1_sol,tf_sol,J_2)

%% Ruhelagenwerte
v_RL = ((2/(3*fy*kappa_ref^2))^(1/4));
l1_RL = -2*fy*kappa_ref^2*v_RL^3;
ay_RL = v_RL^2*kappa_ref;

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





