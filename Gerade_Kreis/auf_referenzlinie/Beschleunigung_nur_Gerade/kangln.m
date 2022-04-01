function f = kangln(X,u,p,region)

ax = u;
s = X(1);
v = X(2);
l1 = X(3);
l2 = X(4);

switch region
    case 1 % Gerade
        f = [0;...
            -(2*p.fy*v^3*p.kapparef_straight^2 + l1)];
    case 2 % Kurve
        f = [0;...
            -(2*p.fy*v^3*p.kapparef_curve^2 + l1)];
end
end