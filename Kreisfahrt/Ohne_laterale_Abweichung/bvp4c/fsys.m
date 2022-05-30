function f = fsys(X,u,p)
ax = u;
sr = X(1);
v = X(2);

f = [v;...
    ax];
end