function f = fsys(X,u,p)
a = u;
s = X(1);
v = X(2);

f = [v;...
    a];
end