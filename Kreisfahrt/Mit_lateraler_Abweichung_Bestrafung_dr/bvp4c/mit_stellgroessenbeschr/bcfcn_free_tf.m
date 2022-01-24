function res = bcfcn_free_tf(X0,Xf,param,p)
res = [X0(1) - p.x0(1);...
    X0(2) - p.x0(2);...
    X0(3) - p.x0(3);...
    X0(4) - p.x0(4);...
    Xf(1) - p.sf;...
    Xf(3) - p.drf;...
    Xf(4) - p.psirf;...
    Xf(6);...
    Xf(8)^2/(2*p.fy*Xf(2)^2) + Xf(5)*Xf(2) - Xf(8)^2/(p.fy*Xf(2)^2) - Xf(8)*Xf(2)*p.kapparef + 1];
end