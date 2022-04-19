% Kreisfahrt Ruhelage Zeitlösung    
% fx = 1; fy = 1; fr = 1; kapparef_curve = 0.1; p.fx = fx; p.fy = fy; p.fr = fr; p.kapparef_curve = kapparef_curve;
syms sr_sym v_sym l1_sym l2_sym ax_sym t_sym tau_sym
X_sym = [sr_sym; v_sym; l1_sym; l2_sym];
U_sym = [ax_sym];

fx = 1; fy = 1; fr = 1; kapparef_straight = 0.0; kapparef_curve = 0.1; 
p.fx = fx; p.fy = fy; p.fr = fr; p.kapparef_straight = kapparef_straight; p.kapparef_curve = kapparef_curve; 

f = [v_sym;...
    -l2_sym/p.fx;...
    0;...
    -(2*p.fy*v_sym^3*p.kapparef_curve^2 + l1_sym)];

X_init_fsolve = [3; 0.1; 0; 0];
fsolve_options = optimoptions('fsolve','MaxIterations',1e5,'MaxFunctionEvaluations',1e5,'Algorithm','levenberg-marquardt');
[RL,fval] = fsolve(@(X) trim_function(X,p),X_init_fsolve,fsolve_options);
X_RL = round(RL(1:3),9); U_RL = round(RL(4),9);
A_sym = jacobian(f,X_sym);
B_sym = jacobian(f,U_sym);
A_RL = double(subs(A_sym,[X_sym; U_sym],[0; X_RL; U_RL]));
B_RL = double(subs(B_sym,[X_sym; U_sym],[0; X_RL; U_RL]));
X0_t = [0; X_RL(1); X_RL(2); X_RL(3)]; U0_t = U_RL;
X_t = expm(A_RL*t_sym)*X0_t;
X_t_simplified = vpa(X_t,8);
% int_anteil = expm(A_RL_transf*(t_sym-tau_sym))*B_RL_transf*U_RL_transf;


clearvars sr_sym v_sym dr_sym psir_sym l1_sym l2_sym l3_sym l4_sym ax_sym kappa_sym tau_sym X_sym U_sym ...
    f X_init_fsolve fsolve_options RL_transf fval X_RL_transf U_RL_transf A_sym B_sym -except t_sym X_t X_t_simplified X0_t U0_t A_RL_transf B_RL_transf
%---------------------------------------------------------------------------
function f = trim_function(X,p)
v=X(1); l1=X(2); l2=X(3); ax=X(4);
f = [ax;...
    0;...
    -(2*p.fy*v^3*p.kapparef_curve^2 + l1);...
    p.fx*ax + l2;...
    1/2*p.fx*ax^2 + 1/2*p.fy*p.kapparef_curve^2*v^4 + l1*v + l2*ax + 1];
end
