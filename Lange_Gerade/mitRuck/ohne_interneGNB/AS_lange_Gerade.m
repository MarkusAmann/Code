% close all; 
clear;

syms c1 c2 cv cs k1 k2;
syms tf;

%% hier parameter vorgeben
fa = 0.001;            % Gewichtung Längsbeschleunigung
fj = 1;            % Gewichtung Ruck.

sf = 1000;   % Länge der Gerade  

s0 = 0; v0 = 10; a0 = 0; vf = 0; af = 0;

%% randwertproblem

lambda1 = c1;
lambda2 = @(t) c2 - c1*t;
      j = @(t) k1*sqrt(fa/fj)*exp(sqrt(fa/fj)*t) - k2*sqrt(fa/fj)*exp(-sqrt(fa/fj)*t) + c1/fa;
      a = @(t) k1*exp(sqrt(fa/fj)*t) + k2*exp(-sqrt(fa/fj)*t) + c1/fa*t - c2/fa;
      v = @(t) k1*sqrt(fj/fa)*exp(sqrt(fa/fj)*t) - k2*sqrt(fj/fa)*exp(-sqrt(fa/fj)*t) + c1/(2*fa)*t^2 - c2/fa*t + cv;
      s = @(t) k1*(fj/fa)*exp(sqrt(fa/fj)*t) + k2*(fj/fa)*exp(-sqrt(fa/fj)*t) + c1/(6*fa)*t^3 - c2/(2*fa)*t^2 + cv*t + cs;

% tf = 40;
% dd = vpasolve([s(0) == s0;...
%     v(0) == v0;...
%     a(0) == a0;...
%     s(tf) == sf;...
% %     a(tf) == af;...
%     j(tf) == 0;...
% %     v(tf) == vf;...
%     lambda2(tf) == 0
%     ], [0,0,0,0,0,0]);

dd = vpasolve([s(0) == s0;...
    v(0) == v0;...
    a(0) == a0;...
    s(tf) == sf;...
%     a(tf) == af;...
    j(tf) == 0;...
%     v(tf) == vf;...
    lambda2(tf) == 0;...
    1/2*fa*a(tf)^2-1/2*fj*j(tf)^2+lambda1*v(tf)+lambda2(tf)*a(tf)+1==0], [0,0,0,0,1e-40,0,19]);

c1_sol = double(dd.c1);
c2_sol = double(dd.c2);
cs_sol = double(dd.cs);
cv_sol = double(dd.cv);
k1_sol = double(dd.k1);
k2_sol = double(dd.k2);
if isfield(dd, 'tf')
    tf_sol = double(dd.tf);
else
    tf_sol = tf;
end

lambda1_sol = matlabFunction(subs(lambda1, c1, c1_sol));
lambda2_sol = matlabFunction(subs(lambda2,{c1,c2},{c1_sol,c2_sol}));
      j_sol = matlabFunction(subs(j,{c1,k1,k2},{c1_sol,k1_sol,k2_sol}));
      a_sol = matlabFunction(subs(a,{c1,c2,k1,k2},{c1_sol,c2_sol,k1_sol,k2_sol}));
      v_sol = matlabFunction(subs(v,{c1,c2,k1,k2,cv},{c1_sol,c2_sol,k1_sol,k2_sol,cv_sol}));
      s_sol = matlabFunction(subs(s,{c1,c2,k1,k2,cv,cs},{c1_sol,c2_sol,k1_sol,k2_sol,cv_sol,cs_sol}));

vf_sol = double(v_sol(tf_sol));

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
if isfield(dd, 'tf')
    tf_sol = double(dd.tf(pos));
else
    tf_sol = tf;
end

lambda1_sol = matlabFunction(subs(lambda1, c1, c1_sol));
lambda2_sol = matlabFunction(subs(lambda2,{c1,c2},{c1_sol,c2_sol}));
      j_sol = matlabFunction(subs(j,{c1,k1,k2},{c1_sol,k1_sol,k2_sol}));
      a_sol = matlabFunction(subs(a,{c1,c2,k1,k2},{c1_sol,c2_sol,k1_sol,k2_sol}));
      v_sol = matlabFunction(subs(v,{c1,c2,k1,k2,cv},{c1_sol,c2_sol,k1_sol,k2_sol,cv_sol}));
      s_sol = matlabFunction(subs(s,{c1,c2,k1,k2,cv,cs},{c1_sol,c2_sol,k1_sol,k2_sol,cv_sol,cs_sol}));

vf_sol = double(v_sol(tf_sol));
l1_sol = double(lambda1_sol());

%% Lösungsvektoren
t_vec = linspace(0,tf_sol,1000);

j_vec = j_sol(t_vec);
a_vec = a_sol(t_vec);
v_vec = v_sol(t_vec);
s_vec = s_sol(t_vec);

l1_vec = l1_sol*ones(size(t_vec));
l2_vec = lambda2_sol(t_vec);
l3_vec = -j_vec/fj;

kappa_vec = zeros(size(t_vec));
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

fprintf('l1_l: %f\ntf: %f\nJ: %f\n',c1_sol,tf_sol,J_1)

c1_calc = (sf + a0*fj/fa - s0 - a0*tf_sol*sqrt(fj/fa) - v0*tf_sol)/(tf_sol^3/(6*fa) - tf_sol^3/(2*fa) + tf_sol^2*sqrt(fj/fa^3) - tf_sol*fj/fa^2 - sqrt(fj^3/fa^5));
tmax_calc = -sqrt(fj/fa)*log(c1_calc/(a0+c1_calc*tf_sol/fa)*sqrt(fj/fa^3));
amax_calc = c1_calc*sqrt(fj/fa^3) - c1_calc*sqrt(fj/fa^3)*log(c1_calc/(a0+c1_calc*tf_sol/fa)*sqrt(fj/fa^3)) - c1_calc*tf_sol/fa;
[amax,idmax] = max(a_vec);
t_vec(idmax);

tmin_calc = sqrt(fj/fa)*(log(sqrt(fj/fa^3)) + sqrt(fj/fa)*tf_sol);
amin_calc = -c1_calc*sqrt(fj/fa^3) - c1_calc*tf_sol/fa + c1_calc*sqrt(fj/fa^3)*(log(sqrt(fj/fa^3)) + sqrt(fj/fa)*tf_sol);
[amin,idmin] = min(a_vec);
t_vec(idmin);

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
xlabel('t [s]')

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





