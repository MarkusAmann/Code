function f = sys_gesamt_free_tf(t,X,param,p)
delta_tf = param;
x = X(1:2); 
u = uopt(X,p);
f = delta_tf*[fsys(x,u,p);...% kanonische Gleichungen
kangln(X,u,p)];
end