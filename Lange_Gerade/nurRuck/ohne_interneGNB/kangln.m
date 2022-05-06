function f = kangln(X,u,p)
j = u;
s = X(1);
v = X(2);
a = X(3);
l1 = X(4);
l2 = X(5);
l3 = X(6);

% J = tf+int(1/2*fa*a^2+1/2*fj*j^2)
f = [0;...
    -l1;...
    -l2];
end