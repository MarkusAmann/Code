function f = fsys(X,u,p)
j = u;
s = X(1);
v = X(2);
ax = X(3);

f = [v;...
    ax;...
    j];
end