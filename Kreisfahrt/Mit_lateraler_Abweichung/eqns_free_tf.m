function F_eqns = eqns_free_tf(X,p) % Gleichungen in Residuenform

F_eqns = 1/2*X(8*(p.N+1)+1)*([X(2)*cos(X(4))/(1-X(3)*p.kapparef);...
    -X(6)/p.fx;...
    X(2)*sin(X(4));...
    -X(8)/(p.fy*X(2)^2) - p.kapparef*X(2)*cos(X(4))/(1-X(3)*p.kapparef);...
    0;...
    -2*X(8)^2/(p.fy*X(2)^2) - X(5)*cos(X(4))/(1-X(3)*p.kapparef) - X(7)*sin(X(4)) + X(8)^2/(p.fy*X(2)^3) + X(8)*p.kapparef*cos(X(4))/(1-X(3)*p.kapparef);...
    -p.kapparef*X(5)*X(2)*cos(X(4))/(1-X(3)*p.kapparef)^2 + p.kapparef^2*X(8)*X(2)*cos(X(4))/(1-X(3)*p.kapparef)^2;...
    X(5)*X(2)*sin(X(4))/(1-X(3)*p.kapparef) - X(7)*X(2)*cos(X(4)) - p.kapparef*X(8)*X(2)*sin(X(4))/(1-X(3)*p.kapparef)] + ...
    [X(10)*cos(X(12))/(1-X(11)*p.kapparef);...
    -X(14)/p.fx;...
    X(10)*sin(X(12));...
    -X(16)/(p.fy*X(10)^2) - p.kapparef*X(10)*cos(X(12))/(1-X(11)*p.kapparef);...
    0;...
    -2*X(16)^2/(p.fy*X(10)^2) - X(13)*cos(X(12))/(1-X(11)*p.kapparef) - X(15)*sin(X(12)) + X(16)^2/(p.fy*X(10)^3) + X(16)*p.kapparef*cos(X(12))/(1-X(11)*p.kapparef);...
    -p.kapparef*X(13)*X(10)*cos(X(12))/(1-X(11)*p.kapparef)^2 + p.kapparef^2*X(16)*X(10)*cos(X(12))/(1-X(11)*p.kapparef)^2;...
    X(13)*X(10)*sin(X(12))/(1-X(11)*p.kapparef) - X(15)*X(10)*cos(X(12)) - p.kapparef*X(16)*X(10)*sin(X(12))/(1-X(11)*p.kapparef)])*p.deltat + ...
    [X(1);X(2);X(3);X(4);X(5);X(6);X(7);X(8)] - [X(9);X(10);X(11);X(12);X(13);X(14);X(15);X(16)];
for k = 1:p.N-1
    F_eqns_k = 1/2*X(8*(p.N+1)+1)*([X(8*k+2)*cos(X(8*k+4))/(1-X(8*k+3)*p.kapparef);...
    -X(8*k+6)/p.fx;...
    X(8*k+2)*sin(X(8*k+4));...
    -X(8*k+8)/(p.fy*X(8*k+2)^2) - p.kapparef*X(8*k+2)*cos(X(8*k+4))/(1-X(8*k+3)*p.kapparef);...
    0;...
    -2*X(8*k+8)^2/(p.fy*X(8*k+2)^2) - X(8*k+5)*cos(X(8*k+4))/(1-X(8*k+3)*p.kapparef) - X(8*k+7)*sin(X(8*k+4)) + X(8*k+8)^2/(p.fy*X(8*k+2)^3) + X(8*k+8)*p.kapparef*cos(X(8*k+4))/(1-X(8*k+3)*p.kapparef);...
    -p.kapparef*X(8*k+5)*X(8*k+2)*cos(X(8*k+4))/(1-X(8*k+3)*p.kapparef)^2 + p.kapparef^2*X(8*k+8)*X(8*k+2)*cos(X(8*k+4))/(1-X(8*k+3)*p.kapparef)^2;...
    X(8*k+5)*X(8*k+2)*sin(X(8*k+4))/(1-X(8*k+3)*p.kapparef) - X(8*k+7)*X(8*k+2)*cos(X(8*k+4)) - p.kapparef*X(8*k+8)*X(8*k+2)*sin(X(8*k+4))/(1-X(8*k+3)*p.kapparef)] + ...
    [X(8*(k+1)+2)*cos(X(8*(k+1)+4))/(1-X(8*(k+1)+3)*p.kapparef);...
    -X(8*(k+1)+6)/p.fx;...
    X(8*(k+1)+2)*sin(X(8*(k+1)+4));...
    -X(8*(k+1)+8)/(p.fy*X(8*(k+1)+2)^2) - p.kapparef*X(8*(k+1)+2)*cos(X(8*(k+1)+4))/(1-X(8*(k+1)+3)*p.kapparef);...
    0;...
    -2*X(8*(k+1)+8)^2/(p.fy*X(8*(k+1)+2)^2) - X(8*(k+1)+5)*cos(X(8*(k+1)+4))/(1-X(8*(k+1)+3)*p.kapparef) - X(8*(k+1)+7)*sin(X(8*(k+1)+4)) + X(8*(k+1)+8)^2/(p.fy*X(8*(k+1)+2)^3) + X(8*(k+1)+8)*p.kapparef*cos(X(8*(k+1)+4))/(1-X(8*(k+1)+3)*p.kapparef);...
    -p.kapparef*X(8*(k+1)+5)*X(8*(k+1)+2)*cos(X(8*(k+1)+4))/(1-X(8*(k+1)+3)*p.kapparef)^2 + p.kapparef^2*X(8*(k+1)+8)*X(8*(k+1)+2)*cos(X(8*(k+1)+4))/(1-X(8*(k+1)+3)*p.kapparef)^2;...
    X(8*(k+1)+5)*X(8*(k+1)+2)*sin(X(8*(k+1)+4))/(1-X(8*(k+1)+3)*p.kapparef) - X(8*(k+1)+7)*X(8*(k+1)+2)*cos(X(8*(k+1)+4)) - p.kapparef*X(8*(k+1)+8)*X(8*(k+1)+2)*sin(X(8*(k+1)+4))/(1-X(8*(k+1)+3)*p.kapparef)])*p.deltat + ...
    [X(8*k+1);X(8*k+2);X(8*k+3);X(8*k+4);X(8*k+5);X(8*k+6);X(8*k+7);X(8*k+8)] - [X(8*(k+1)+1);X(8*(k+1)+2);X(8*(k+1)+3);X(8*(k+1)+4);X(8*(k+1)+5);X(8*(k+1)+6);X(8*(k+1)+7);X(8*(k+1)+8)];

    F_eqns = [F_eqns; F_eqns_k];
end
X_eqns_t0 = [X(1) - p.x0(1);...
    X(2) - p.x0(2);...
    X(3) - p.x0(3);...
    X(4) - p.x0(4)];

% wenn dr und psir bei tf fest sind
G_eqns_tf = [X(8*p.N+1) - p.sf;...
    X(8*p.N+3) - p.drf;...
    X(8*p.N+4) - p.psirf;...
    X(8*p.N+6);...
    X(8*p.N+8)^2/(2*p.fy*X(8*p.N+2)^2) + X(8*p.N+5)*X(8*p.N+2) - X(8*p.N+8)^2/(p.fy*X(8*p.N+2)^2) - X(8*p.N+8)*X(8*p.N+2)*p.kapparef + 1];

% wenn dr und psir bei tf frei sind
% G_eqns_tf = [X(8*p.N+1) - p.sf;...
%     X(8*p.N+6);...
%     X(8*p.N+7);...
%     X(8*p.N+8)];
F_eqns = [F_eqns; X_eqns_t0; G_eqns_tf];
end
