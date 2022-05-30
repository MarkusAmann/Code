function f = kangln(X,u,p)

ax = u;
sr = X(1);
v = X(2);
l1 = X(3);
l2 = X(4);


% ohne Zustandsbestrafung im Gütefunktional
f = [0;...
    -(2*p.fy*v^3*p.kapparef^2 + l1)];
end
