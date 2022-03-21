function f = fsys_straight(X,u,p)
ax = u;
sr = X(1);
v = X(2);
dr = X(3);
psir = X(4);

f = [v*cos(psir)/(1-dr*p.kapparef_straight);...
    ax;...
    v*sin(psir);...
    p.kapparef_straight*v - p.kapparef_straight*v*cos(psir)/(1-dr*p.kapparef_straight)];

end