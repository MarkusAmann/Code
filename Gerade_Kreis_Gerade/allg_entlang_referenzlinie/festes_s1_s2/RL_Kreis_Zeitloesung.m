% Kreisfahrt Ruhelage Zeitlösung    
% fx = 1; fy = 1; fr = 1; kapparef_curve = 0.1; p.fx = fx; p.fy = fy; p.fr = fr; p.kapparef_curve = kapparef_curve;
syms sr_sym v_sym dr_sym psir_sym l1_sym l2_sym l3_sym l4_sym ax_sym kappa_sym t_sym tau_sym
X_sym = [sr_sym; v_sym; dr_sym; psir_sym; l1_sym; l2_sym; l3_sym; l4_sym];
U_sym = [ax_sym; kappa_sym];

f = [v_sym*cos(psir_sym)/(1-dr_sym*p.kapparef_curve);...
    ax_sym;...
    v_sym*sin(psir_sym);...
    kappa_sym*v_sym - p.kapparef_curve*v_sym*cos(psir_sym)/(1-dr_sym*p.kapparef_curve);...
    0;...
    -(2*p.fy*v_sym^3*kappa_sym^2 + l1_sym*cos(psir_sym)/(1-dr_sym*p.kapparef_curve) + l3_sym*sin(psir_sym) + l4_sym*kappa_sym - l4_sym*p.kapparef_curve*cos(psir_sym)/(1-dr_sym*p.kapparef_curve));...
    -(p.fr*dr_sym + p.kapparef_curve*l1_sym*v_sym*cos(psir_sym)/(1-dr_sym*p.kapparef_curve)^2 - p.kapparef_curve^2*l4_sym*v_sym*cos(psir_sym)/(1-dr_sym*p.kapparef_curve)^2);...
    -(l3_sym*v_sym*cos(psir_sym) - l1_sym*v_sym*sin(psir_sym)/(1-dr_sym*p.kapparef_curve) + p.kapparef_curve*l4_sym*v_sym*sin(psir_sym)/(1-dr_sym*p.kapparef_curve))];


X_init_fsolve = [3; 0.1; 0; -0.1; 0; 0; -2; 0; p.kapparef_curve];
fsolve_options = optimoptions('fsolve','MaxIterations',1e5,'MaxFunctionEvaluations',1e5,'Algorithm','levenberg-marquardt');
[RL_transf,fval] = fsolve(@(X) trim_function(X,p),X_init_fsolve,fsolve_options);
X_RL_transf = round(RL_transf(1:7),9); U_RL_transf = round(RL_transf(8:9),9);
A_sym = jacobian(f,X_sym);
B_sym = jacobian(f,U_sym);
A_RL_transf = double(subs(A_sym,[X_sym; U_sym],[0; X_RL_transf; U_RL_transf]));
B_RL_transf = double(subs(B_sym,[X_sym; U_sym],[0; X_RL_transf; U_RL_transf]));
X0_t = [0; X_RL_transf]; U0_t = U_RL_transf;
X_t = expm(A_RL_transf*t_sym)*X0_t;
X_t_simplified = vpa(X_t,8);
% int_anteil = expm(A_RL_transf*(t_sym-tau_sym))*B_RL_transf*U_RL_transf;


clearvars sr_sym v_sym dr_sym psir_sym l1_sym l2_sym l3_sym l4_sym ax_sym kappa_sym tau_sym X_sym U_sym ...
    f X_init_fsolve fsolve_options RL_transf fval X_RL_transf U_RL_transf A_sym B_sym -except t_sym X_t X_t_simplified X0_t U0_t A_RL_transf B_RL_transf
%---------------------------------------------------------------------------
function f = trim_function(X,p)
v=X(1); dr=X(2); psir=X(3); l1=X(4); l2=X(5); l3=X(6); l4=X(7); ax=X(8); kappa=X(9);
kappa_transf = kappa + p.kapparef_curve*cos(psir)/(1-dr*p.kapparef_curve);
f = [ax;...
    v*sin(psir);...
    kappa_transf*v - p.kapparef_curve*v*cos(psir)/(1-dr*p.kapparef_curve);...
    -(2*p.fy*v^3*kappa_transf^2 + l1*cos(psir)/(1-dr*p.kapparef_curve) + l3*sin(psir) + l4*kappa_transf - l4*p.kapparef_curve*cos(psir)/(1-dr*p.kapparef_curve));...
    -(p.fr*dr + p.kapparef_curve*l1*v*cos(psir)/(1-dr*p.kapparef_curve)^2 - p.kapparef_curve^2*l4*v*cos(psir)/(1-dr*p.kapparef_curve)^2);...
    -(l3*v*cos(psir) - l1*v*sin(psir)/(1-dr*p.kapparef_curve) + p.kapparef_curve*l4*v*sin(psir)/(1-dr*p.kapparef_curve));...
    p.fx*ax + l2;...
    kappa_transf*p.fy*v^4 + l4*v;...
    1/2*p.fr*dr^2 + 1/2*p.fx*ax^2 + 1/2*p.fy*kappa_transf^2*v^4 + l1*v*cos(psir)/(1-dr*p.kapparef_curve) + + l2*ax + l3*v*sin(psir) + l4*v*(kappa_transf - p.kapparef_curve*cos(psir)/(1-dr*p.kapparef_curve)) + 1];
end
