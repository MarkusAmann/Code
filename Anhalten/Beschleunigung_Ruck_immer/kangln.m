function f = kangln(X,u,p,region)

j = u;
s = X(1);
v = X(2);
ax = X(3);
l1 = X(4);
l2 = X(5);
l3 = X(6);

f = [0;...
    -l1;...
    -(p.fa*ax + l2)];
end