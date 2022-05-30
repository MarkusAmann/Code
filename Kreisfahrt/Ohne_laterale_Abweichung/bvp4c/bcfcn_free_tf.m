function res = bcfcn_free_tf(X0,Xf,param,p)
delta_t = param; 
srf = Xf(1);
vf = Xf(2);
l1f = Xf(3);
l2f = Xf(4);
uf = uopt(Xf,p);
axf = uf;

% H(tf)+1=0 ohne Zustandsbestrafung im Gütefunktional
% J=tf+int(1/2*fx*ax^2+1/2*fy*kappa^2*v^4)
H_tf = 1/2*p.fx*axf^2 + 1/2*p.fy*p.kapparef^2*vf^4 + l1f*vf + l2f*axf;

% nur sf ist festgelegt, vf, drf und psirf sind frei
res = [X0(1) - p.x0(1);...
    X0(2) - p.x0(2);...
    srf - p.sf;...
    l2f;...
    H_tf + 1];

end
