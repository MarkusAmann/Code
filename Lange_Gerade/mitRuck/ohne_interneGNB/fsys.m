function f = fsys(X,u,p)
j = u;
s = X(1);
v = X(2);
a = X(3);

f = [v;...
    a;...
    j];
end