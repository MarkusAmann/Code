function f = sys_gesamt_free_tf(t,X,region,param,p)
nu_tilde = param(1);
delta_t1 = param(2);

x = X(1:3); 
u = uopt(X,p);
f = delta_t1*[fsys(x,u,p);...% kanonische Gleichungen
kangln(X,u,p)];
end