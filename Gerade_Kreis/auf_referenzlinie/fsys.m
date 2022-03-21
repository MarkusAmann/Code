function f = fsys(X,u,p)
ax = u;
s = X(1);
v = X(2);

f = [v;...
    ax];
end