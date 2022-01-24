function f = sys_gesamt_fix_tf(t,X,p)
x = X(1:4); 
u = uopt(X,p);
f = [fsys(x,u,p);...% kanonische Gleichungen
kangln(X,u,p)];
end