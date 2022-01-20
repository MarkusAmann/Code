function [Xopt,fval,exitflag,output] = indirect_discretization(p,varargin)

if nargin == 1
    X0 = [0;0.001;0;0;0;0;0;0];
    X_mat = repmat(X0,p.N+1);
    X_all_init = X_mat(:,1);
    % X_all_init = zeros(8*(p.N+1),1);
    opt = optimoptions('fsolve','Display','iter','MaxFunctionEvaluations',p.maxFunEval,'MaxIterations',p.maxIter); % Optionen
    [Xopt,fval,exitflag,output] = fsolve(@eqns_fix_tf,X_all_init,opt,p); % Numerische Losung mit fsolve
else
    tf_free = varargin;
    X0 = [0;1;0;0;0;0;0;0];
    X_mat = repmat(X0,p.N+1);
    X_all_init = X_mat(:,1);
    X_all_init = [X_all_init;0.5]; % init guess for delta*tf0 = tf
    % X_all_init = zeros(8*(p.N+1),1);
    opt = optimoptions('fsolve','Display','iter','MaxFunctionEvaluations',p.maxFunEval,'MaxIterations',p.maxIter); % Optionen
    [Xopt,fval,exitflag,output] = fsolve(@eqns_free_tf,X_all_init,opt,p); % Numerische Losung mit fsolve
end
