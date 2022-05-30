function [Xopt,fval,exitflag,output] = indirect_discretization(p)
t = linspace(p.t0, p.tf, p.N+1);
deltat = mean(diff(t));
p.deltat = deltat;
% X0 = [x0.'; l0.'];
% X_all_init = X0;
% opt = optimset('Display','iter'); % Optionen
% for i = 1:N
%     X_i_init = fsolve(@X_init,X0,opt,p);
%     X_all_init = [X_all_init; X_i_init];
%     X0 = X_i_init;
%     p.x0 = X0(1:2).';
%     p.l0 = X0(3:4).';
% end
% p.x0 = x0; p.l0 = l0;
X_all_init = zeros(4*(p.N+1),1);
opt = optimoptions('fsolve','Display','iter','MaxFunctionEvaluations',1e06,'MaxIterations',1e05); % Optionen
[Xopt,fval,exitflag,output] = fsolve(@eqns,X_all_init,opt,p); % Numerische Losung mit fsolve
end
