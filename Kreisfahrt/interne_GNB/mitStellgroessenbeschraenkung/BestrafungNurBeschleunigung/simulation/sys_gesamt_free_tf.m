function f = sys_gesamt_free_tf(t,X,param,p)
x = X(1:4); 
u = uopt(X,p);
f = param*[fsys(x,u,p);...% kanonische Gleichungen
kangln(X,u,p)];
end