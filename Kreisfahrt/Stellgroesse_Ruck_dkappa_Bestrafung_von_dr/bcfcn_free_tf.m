function res = bcfcn_free_tf(X0,Xf,param,p)
delta_t = param; 
srf = Xf(1);
vf = Xf(2);
af = Xf(3);
drf = Xf(4);
psirf = Xf(5);
kappaf = Xf(6);
l1f = Xf(7);
l2f = Xf(8);
l3f = Xf(9);
l4f = Xf(10);
l5f = Xf(11);
l6f = Xf(12);
uf = uopt(Xf,p);
jf = uf(1);
dkappaf = uf(2);

jyf = dkappaf*vf^2 + 2*kappaf*af*vf;
switch p.use_dr
    case 0
        % H(tf)+1=0 ohne Zustandsbestrafung im Gütefunktional
        % J=tf+int(1/2*fjx*j^2+1/2*fax*a^2+1/2*fay*kappa^2*v^4)
        H_tf = 1/2*p.fjx*jf^2 + 1/2*p.fax*af^2 + 1/2*p.fay*kappaf^2*vf^4 + 1/2*p.fjy*jyf^2 + l1f*vf*cos(psirf)/(1-drf*p.kapparef) + l2f*af + l3f*jf + l4f*vf*sin(psirf) + ...
            l5f*vf*(kappaf - p.kapparef*cos(psirf)/(1-drf*p.kapparef)) + l6f*dkappaf;
    case 1
        % H(tf)+1=0 Bestrafung von dr im Gütefunktional
        % J=tf+int(1/2*fjx*j^2+1/2*fax*a^2+1/2*fay*kappa^2*v^4+1/2*fr*dr^2)
        H_tf = 1/2*p.fr*drf^2 + 1/2*p.fjx*jf^2 + 1/2*p.fax*af^2 + 1/2*p.fay*kappaf^2*vf^4 + 1/2*p.fjy*jyf^2 + l1f*vf*cos(psirf)/(1-drf*p.kapparef) + l2f*af + l3f*jf + l4f*vf*sin(psirf) + ...
            l5f*vf*(kappaf - p.kapparef*cos(psirf)/(1-drf*p.kapparef)) + l6f*dkappaf;
end

% sf und psirf sind festgelegt, vf, af, drf und kappaf sind frei
res = [X0(1) - p.x0(1);...
    X0(2) - p.x0(2);...
    X0(3) - p.x0(3);...
    X0(4) - p.x0(4);...
    X0(5) - p.x0(5);...
    X0(6) - p.x0(6);...
    srf - p.sf;...
    psirf - p.psirf;...
    l2f;...
    l3f;...
    l4f;...
    l6f;...
    H_tf + 1];

end
