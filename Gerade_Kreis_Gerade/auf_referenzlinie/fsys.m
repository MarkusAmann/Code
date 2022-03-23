function f = fsys(X,u,p,region)
ax = u;
s = X(1);
v = X(2);

f = [v;...
    ax];
end