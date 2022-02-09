% clc
clear all
close all

%% Parameter
x0 = [0 1 0].'; l0 = [0 0 0].'; %l0 = 0.1*randn(4,1);
jlim = 0.8; 
umax = jlim; umin = -jlim; use_umax = 0;
t0 = 0; tf = 104.1238; N = 100; fa = 1; fj = 1; sf = 5000; 
tf_free = 1; % ACHTUNG: ich habe den Eindruck, dass die Optimierung bei aktiver Stellgrößenbeschränkung und freier Endzeit Probleme hat
% aufgrund des instabilen Teils der Lösung
p.use_umax = use_umax; p.umax = umax; p.umin = umin; p.fa = fa; p.fj = fj; p.sf = sf; 
p.x0 = x0; p.l0 = l0; p.t0 = t0; p.tf = tf; p.tf_free = tf_free; p.N = N;  


%% System
% x=[s;v;a], u=j
% J=tf+int(1/2*fa*a^2+1/2*fj*j^2)
A = [0 1 0 0 0 0;...
    0 0 1 0 0 0;...
    0 0 0 0 0 1;...
    0 0 0 0 0 0 ;...
    0 0 0 -1 0 0;...
    0 0 fa/fj 0 1/fj 0];

%% analytische Lösung für a(t)
% dj/dt = d^2a/dt^2 = fa*a/fj + l2/fj, mit l2 = -l1*t + c2
% --> d^2a/dt^2 = fa*a/fj - l1*t/fj + c2/fj
% --> a(t) = k1*e^(-sqrt(fa/fj)*t) + k2*e^(sqrt(fa/fj)*t) + l1*t/fa - c2/fa
s0 = x0(1);
v0 = x0(2);
a0 = x0(3);
syms a(tsym) fasym fjsym c2sym l1sym cvsym cssym tfsym C1 C2
eqn = diff(a,tsym,2) == fasym/fjsym*a - l1sym/fjsym*tsym + c2sym/fjsym;
% cond = a(0) == a0;
S = dsolve(eqn);
a_t = S;
j_t = diff(a_t,tsym);
v_t = int(a_t,tsym) + cvsym;
s_t = int(v_t,tsym) + cssym;
l1_t = l1sym;
l2_t = int(-l1_t,tsym) + c2sym;
l3_t = -fjsym*j_t;
% H_tf = 1/2*fasym*a_t^2 + 1/2*fjsym*j_t^2 + l1_t*v_t + l2_t*a_t + l3_t*j_t;
H_tf = 1/2*fasym*a_t^2 + l1_t*v_t;


% syms k1 k2 t tf l1 c2 cv cs
% j_t = k1*sqrt(fa/fj)*exp(sqrt(fa/fj)*t) - k2*sqrt(fa/fj)*exp(-sqrt(fa/fj)*t) + l1/fa;
% a_t = k1*exp(sqrt(fa/fj)*t) + k2*exp(-sqrt(fa/fj)*t) + l1/fa*t - c2/fa;
% v_t = k1*sqrt(fj/fa)*exp(sqrt(fa/fj)*t) - k2*sqrt(fj/fa)*exp(-sqrt(fa/fj)*t) + l1/(2*fa)*t^2 - c2/fa*t + cv;
% s_t = k1*fj/fa*exp(sqrt(fa/fj)*t) + k2*fj/fa*exp(-sqrt(fa/fj)*t) + l1/(6*fa)*t^3 - c2/(2*fa)*t^2 + cv*t + cs;
% l1_t = l1;
% l2_t = -l1_t*t + c2;
% l3_t = -fj*j_t;
% H_t = 1/2*fa*a_t^2 + 1/2*fj*j_t^2 + l1_t*v_t + l2_t*a_t + l3_t*j_t;
% H_tf = 1/2*fa*a_t^2 + l1_t*v_t;

% parameter_eqns = [subs(s_t,t,0) - s0;...
%     subs(v_t,t,0) - v0;...
%     subs(a_t,t,0) - a0;...
%     subs(s_t,t,tf) - sf;...
%     subs(l2_t,t,tf);...
%     subs(l3_t,t,tf);...
%     subs(H_tf,t,tf) + 1];
% sol = solve(parameter_eqns==0,[k1 k2 tf l1 c2 cv cs]);
% funhan = @(x)fun(x,p);
% [paramsol,fsol] = fsolve(funhan,[-0.1 -1 100 -0.01 -1.33 4 -1.33]);

%% Optimierung
bvpoptions = bvpset('RelTol',1e-5,'Stats','on');
switch tf_free
    case 0
        t = linspace(p.t0, p.tf, p.N+1);
        deltat = mean(diff(t));
        p.deltat = deltat;
        init_guess = @(x)guess_fix_tf(t,p);
        solinit = bvpinit(t,init_guess);
        sol = bvp4c(@sys_gesamt_fix_tf, @bcfcn_fix_tf, solinit, bvpoptions, p);
        % optimal states
        sol_mesh = sol.x;
        sopt = sol.y(1,:);
        vopt = sol.y(2,:);
        aopt = sol.y(3,:);
        l1opt = sol.y(4,:);
        l2opt = sol.y(5,:);
        l3opt = sol.y(6,:);

        parameter_eqns = [subs(s_t,tsym,0) - s0;...
            subs(v_t,tsym,0) - v0;...
            subs(a_t,tsym,0) - a0;...
            subs(s_t,tsym,tf) - sf;...
            subs(l2_t,tsym,tf);...
            subs(l3_t,tsym,tf)];
        sym_vars = [C1 C2 l1sym c2sym cvsym cssym];
        sym_params = solve(parameter_eqns==0,sym_vars);
        C1_num = double(subs(sym_params.C1,[fasym, fjsym],[fa, fj]));
        C2_num = double(subs(sym_params.C2,[fasym, fjsym],[fa, fj]));
        l1_num = double(subs(sym_params.l1sym,[fasym, fjsym],[fa, fj]));
        c2_num = double(subs(sym_params.c2sym,[fasym, fjsym],[fa, fj]));
        cv_num = double(subs(sym_params.cvsym,[fasym, fjsym],[fa, fj]));
        cs_num = double(subs(sym_params.cssym,[fasym, fjsym],[fa, fj]));
        num_vars = [C1_num, C2_num, l1_num, c2_num, cv_num, cs_num];
        a_t_num = subs(a_t,[sym_vars fasym fjsym],[num_vars fa fj]);
        v_t_num = subs(v_t,[sym_vars fasym fjsym],[num_vars fa fj]);
        s_t_num = subs(s_t,[sym_vars fasym fjsym],[num_vars fa fj]);
        j_t_num = subs(j_t,[sym_vars fasym fjsym],[num_vars fa fj]);
        l3_t_num = subs(l3_t,[sym_vars fasym fjsym],[num_vars fa fj]);
        l2_t_num = subs(l2_t,[sym_vars fasym fjsym],[num_vars fa fj]);
        l1_t_num = subs(l1_t,[sym_vars fasym fjsym],[num_vars fa fj]);
        a_put_t = double(subs(a_t_num,tsym,sol_mesh));
        v_put_t = double(subs(v_t_num,tsym,sol_mesh));
        s_put_t = double(subs(s_t_num,tsym,sol_mesh));
        j_put_t = double(subs(j_t_num,tsym,sol_mesh));
        l1_put_t = double(subs(l1_t_num,tsym,sol_mesh));
        l2_put_t = double(subs(l2_t_num,tsym,sol_mesh));
        l3_put_t = double(subs(l3_t_num,tsym,sol_mesh));
        
        % Nullstellen von j_t berechnen
        zero_j = double(solve(j_t_num == 0,tsym))
        % a(t) an den Nullstellen berechnen
        a_zero_j = double(subs(a_t_num,tsym,zero_j))

        % mit C1 = 0 folgt 
        tmax = -sqrt(fj/fa)*log(l1_num/C2_num*sqrt(fj/fa^3));
        amax = l1_num*sqrt(fj/fa^3)-l1_num*sqrt(fj/fa^3)*log(l1_num/C2_num*sqrt(fj/fa^3))-c2_num/fa;
        l1_approx = (sf+a0*fj/fa-s0-a0*tf*sqrt(fj/fa)-v0*tf)/(tf^3/(6*fa)-tf^3/(2*fa)+tf^2*sqrt(fj/fa^3)-tf*fj/fa^2);
        c2_approx = l1_approx*tf;
        cs_approx = s0 - a0*fj/fa - l1_approx*tf*fj/fa^2;
        cv_approx = v0 + a0*sqrt(fj/fa) + l1_approx*tf*sqrt(fj/fa^3);
        C2_approx = l1_approx*tf/fa + a0;
        amax_of_l1 = l1_num*sqrt(fj/fa^3)-l1_num*sqrt(fj/fa^3)*log(l1_num/(l1_num*tf+a0*fa)*sqrt(fj/fa))-l1_num*tf/fa;

        % Analyse von amax in Abhängigkeit von tf und fj bzw. fa, wobei
        % fj==fa (fj = fa = F) gilt
        F = 0.1:0.1:10;
        T = 100:1:1000;
        [F_grid,T_grid] = meshgrid(F,T);
        % fj und fa sind schon gegeneinander gekürzt
        lam1 = (sf+a0-s0-a0*T_grid-v0*T_grid)./(T_grid.^3./(6*F_grid)-T_grid.^3./(2*F_grid)+T_grid.^2./F_grid-T_grid./F_grid);
        AMAX = lam1./F_grid - lam1./F_grid.*log(lam1./(lam1.*T_grid + a0*F_grid)) - lam1.*T_grid./F_grid;
        figure
        mesh(T_grid,F_grid,AMAX);
        title('maximum acc, for f_a/f_j==1')
        xlabel('t_f [s]')
        ylabel('f_a==f_j')
        zlabel('a_{max} [m/s^2]')

    case 1
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
        delta_opt = sol.parameters;
        tf_opt = delta_opt*p.tf;
        sol_mesh = sol_mesh*delta_opt;

%         parameter_eqns = [subs(s_t,tsym,0) - s0;...
%             subs(v_t,tsym,0) - v0;...
%             subs(a_t,tsym,0) - a0;...
%             subs(s_t,tsym,tf) - sf;...
%             subs(l2_t,tsym,tf);...
%             subs(l3_t,tsym,tf);...
%             subs(H_tf,tsym,tf) + 1];
%         sol = solve(parameter_eqns==0,[C1 C2 tfsym l1sym c2sym cvsym cssym]);

end

% optimal control inputs
for i=1:length(sol_mesh)
   jopt(:,i) = uopt(sol.y(:,i),p); % Steuerung
end

%%
figure
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

figure
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

figure 
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