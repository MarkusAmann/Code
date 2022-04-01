function f = sys_gesamt_free_tf(t,X,region,param,p)
nu_tilde1 = param(1);
nu_tilde2 = param(2);
delta_t1 = param(3);
delta_t2 = param(4);

x = X(1:2); 
u = uopt(X,p);
switch region
    case 1
        f = delta_t1*[fsys(x,u,p);...% kanonische Gleichungen
            kangln(X,u,p,region)];
    case 2 
        f = delta_t2*[fsys(x,u,p);...% kanonische Gleichungen
            kangln(X,u,p,region)];
end
end