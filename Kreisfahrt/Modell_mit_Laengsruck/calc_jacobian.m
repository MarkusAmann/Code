syms X1 X2 X3 X4 X5 X6 X7 X8 tf k fx fy fr
F = tf*[X2*cos(X4)/(1-X3*k);...
    -X6/fx;...
    X2*sin(X4);...
    -X8/(fy*X2^2) - k*X2*cos(X4)/(1-X3*k);...
    0;...
    -2*X8^2/(fy*X2^2) - X5*cos(X4)/(1-X3*k) - X7*sin(X4) + X8^2/(fy*X2^3) + X8*k*cos(X4)/(1-X3*k);...
    -fr*X3 - k*X5*X2*cos(X4)/(1-X3*k)^2 + k^2*X8*X2*cos(X4)/(1-X3*k)^2;...
    X5*X2*sin(X4)/(1-X3*k) - X7*X2*cos(X4) - k*X8*X2*sin(X4)/(1-X3*k)];
dFdX = jacobian(F,[X1,X2,X3,X4,X5,X6,X7,X8]);