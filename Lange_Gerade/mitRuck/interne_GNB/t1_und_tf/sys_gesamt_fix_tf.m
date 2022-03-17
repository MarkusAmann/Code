function f = sys_gesamt_fix_tf(t,X,region,param,p)
nu_tilde = param(1);
delta_t1 = param(2);
x = X(1:3); 
u = uopt(X,p);

switch region
    case 1
        f = delta_t1*[fsys(x,u,p);...% kanonische Gleichungen
            kangln(X,u,p)];
    case 2 
        f = [fsys(x,u,p);...% kanonische Gleichungen
            kangln(X,u,p)];
end
end