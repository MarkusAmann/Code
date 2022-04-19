function f = kangln(X,u,p,region)

j = u;
s = X(1);
v = X(2);
ax = X(3);
l1 = X(4);
l2 = X(5);
l3 = X(6);

switch region
    case 1 % Gerade
        f = [0;...
            -(2*p.fy*v^3*p.kapparef_straight^2 + l1);...
            -(p.fx*ax + l2)];
    case 2 % Kurve
        f = [0;...
            -(2*p.fy*v^3*p.kapparef_curve^2 + l1);...
            -(p.fx*ax + l2)];
end
end