function [Xopt,fval,exitflag,output] = indirect_discretization(p)
X0 = [0;0.001;0;0;0;0;0;0];
X_mat = repmat(X0,p.N+1);
X_all_init = X_mat(:,1);
% X_all_init = zeros(8*(p.N+1),1);
opt = optimoptions('fsolve','Display','iter','MaxFunctionEvaluations',p.maxFunEval,'MaxIterations',p.maxIter); % Optionen
[Xopt,fval,exitflag,output] = fsolve(@eqns,X_all_init,opt,p); % Numerische Losung mit fsolve
end
