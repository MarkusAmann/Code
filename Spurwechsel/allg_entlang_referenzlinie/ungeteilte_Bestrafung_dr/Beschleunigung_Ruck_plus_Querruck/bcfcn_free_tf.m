function res = bcfcn_free_tf(XL,XR,param,p)
nu_tilde1 = param(1);
nu_tilde2 = param(2);
delta_t1 = param(3);
delta_t2 = param(4);

X0 = XL(:,1);
Xf = XR(:,end);
X_internalGNB_L = [XL(1:6,end); XR(7:12,1)];
X_internalGNB_R = [XR(1:6,1); XL(7:12,end)];
% X_internalGNB_L = XL(:,end);
% X_internalGNB_R = XR(:,1);

% Endzustand
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

% Stellgrößen tf, Kurve
uf = uopt(Xf,p);
jf = uf(1);
dkappaf = uf(2);

% linker Randzustand von t1
sL = X_internalGNB_L(1);
vL = X_internalGNB_L(2);
aL = X_internalGNB_L(3);
drL = X_internalGNB_L(4);
psirL = X_internalGNB_L(5);
kappaL = X_internalGNB_L(6);
l1L = X_internalGNB_L(7);
l2L = X_internalGNB_L(8);
l3L = X_internalGNB_L(9);
l4L = X_internalGNB_L(10);
l5L = X_internalGNB_L(11);
l6L = X_internalGNB_L(12);

% rechter Randzustand von t1
sR = X_internalGNB_R(1);
vR = X_internalGNB_R(2);
aR = X_internalGNB_R(3);
drR = X_internalGNB_R(4);
psirR = X_internalGNB_R(5);
kappaR = X_internalGNB_R(6);
l1R = X_internalGNB_R(7);
l2R = X_internalGNB_R(8);
l3R = X_internalGNB_R(9);
l4R = X_internalGNB_R(10);
l5R = X_internalGNB_R(11);
l6R = X_internalGNB_R(12);

% Stellgrößen linker und rechter Rand von t1
uL = uopt(X_internalGNB_L,p);
jL = uL(1);
dkappaL = uL(2); 

uR = uopt(X_internalGNB_R,p);
jR = uR(1);
dkappaR = uR(2); 

jyf = dkappaf*vf^2 + 2*kappaf*af*vf;
jyL = dkappaL*vL^2 + 2*kappaL*aL*vL;
jyR = dkappaR*vR^2 + 2*kappaR*aR*vR;

switch p.use_dr
    case 1
        % H(tf)+1=0 mit Zustandsbestrafung im Gütefunktional
        H_tf = 1/2*p.fr*drf^2 + 1/2*p.fjx*jf^2 + 1/2*p.fax*af^2 + 1/2*p.fjy*jyf^2 + 1/2*p.fay*kappaf^2*vf^4 + ...
            l1f*vf*cos(psirf)/(1-drf*p.kapparef) + l2f*af + l3f*jf + l4f*vf*sin(psirf) + l5f*vf*(kappaf - p.kapparef*cos(psirf)/(1-drf*p.kapparef)) + l6f*dkappaf;
        H_L = 1/2*p.fr*drL^2 + 1/2*p.fjx*jL^2 + 1/2*p.fax*aL^2 + 1/2*p.fjy*jyL^2 + 1/2*p.fay*kappaL^2*vL^4 + ...
            l1L*vL*cos(psirL)/(1-drL*p.kapparef) + l2L*aL + l3L*jL + l4L*vL*sin(psirL) + l5L*vL*(kappaL - p.kapparef*cos(psirL)/(1-drL*p.kapparef)) + l6L*dkappaL;
        H_R = 1/2*p.fr*drR^2 + 1/2*p.fjx*jR^2 + 1/2*p.fax*aR^2 + 1/2*p.fjy*jyR^2 + 1/2*p.fay*kappaR^2*vR^4 + ...
            l1R*vR*cos(psirR)/(1-drR*p.kapparef) + l2R*aR + l3R*jR + l4R*vR*sin(psirR) + l5R*vR*(kappaR - p.kapparef*cos(psirR)/(1-drR*p.kapparef)) + l6R*dkappaR;
    case 0
        % H(tf)+1=0 ohne Zustandsbestrafung im Gütefunktional
        H_tf = 1/2*p.fjx*jf^2 + 1/2*p.fax*af^2 + 1/2*p.fjy*jyf^2 + 1/2*p.fay*kappaf^2*vf^4 + ...
            l1f*vf*cos(psirf)/(1-drf*p.kapparef) + l2f*af + l3f*jf + l4f*vf*sin(psirf) + l5f*vf*(kappaf - p.kapparef*cos(psirf)/(1-drf*p.kapparef)) + l6f*dkappaf;
        H_L = 1/2*p.fjx*jL^2 + 1/2*p.fax*aL^2 + 1/2*p.fjy*jyL^2 + 1/2*p.fay*kappaL^2*vL^4 + ...
            l1L*vL*cos(psirL)/(1-drL*p.kapparef) + l2L*aL + l3L*jL + l4L*vL*sin(psirL) + l5L*vL*(kappaL - p.kapparef*cos(psirL)/(1-drL*p.kapparef)) + l6L*dkappaL;
        H_R = 1/2*p.fjx*jR^2 + 1/2*p.fax*aR^2 + 1/2*p.fjy*jyR^2 + 1/2*p.fay*kappaR^2*vR^4 + ...
            l1R*vR*cos(psirR)/(1-drR*p.kapparef) + l2R*aR + l3R*jR + l4R*vR*sin(psirR) + l5R*vR*(kappaR - p.kapparef*cos(psirR)/(1-drR*p.kapparef)) + l6R*dkappaR;
end

% sf, psirf sind festgelegt, vf, af, drf, kappaf sind frei
% res = [X0(1) - p.x0(1);...
%     X0(2) - p.x0(2);...
%     X0(3) - p.x0(3);...
%     X0(4) - p.x0(4);...
%     X0(5) - p.x0(5);...
%     X0(6) - p.x0(6);...
%     srf - p.sf;...
%     psirf - p.psirf;...
%     l2f;...
%     l3f;
%     l4f;...
%     l6f;...
%     H_tf + 1;...
%     sL - sR;... % Stetigkeit am internen Randwert
%     vL - vR;...
%     drL - drR;...
%     aL - aR;...
%     psirL - psirR;...
%     kappaL - kappaR;...
%     l1L - l1R - 2*nu_tilde1;... % Transversalitätsbedingung für l1 (l1(t1-) = l1(t1+) + dg/dx|t1*nu_tilde)
%     l2L - l2R;...
%     l3L - l3R;...
%     l4L - l4R - 2*nu_tilde2;...
%     l5L - l5R;...
%     l6L - l6R;...
%     sL - p.s1;... % zusätzliche interne Randbedingung
%     drL - p.dr1;...
%     H_L - H_R + 0]; % Transversalitätsbedingung für Umschaltpunkt

% sf, drf, psirf sind festgelegt, vf, af, kappaf sind frei
res = [X0(1) - p.x0(1);...
    X0(2) - p.x0(2);...
    X0(3) - p.x0(3);...
    X0(4) - p.x0(4);...
    X0(5) - p.x0(5);...
    X0(6) - p.x0(6);...
    srf - p.sf;...
    drf - p.drf;...
    psirf - p.psirf;...
    l2f;...
    l3f;
    l6f;...
    H_tf + 1;...
    sL - sR;... % Stetigkeit am internen Randwert
    vL - vR;...
    drL - drR;...
    aL - aR;...
    psirL - psirR;...
    kappaL - kappaR;...
    l1L - l1R - 2*nu_tilde1;... % Transversalitätsbedingung für l1 (l1(t1-) = l1(t1+) + dg/dx|t1*nu_tilde)
    l2L - l2R;...
    l3L - l3R;...
    l4L - l4R - 2*nu_tilde2;...
    l5L - l5R;...
    l6L - l6R;...
    sL - p.s1;... % zusätzliche interne Randbedingung
    drL - p.dr1;...
    H_L - H_R + 0]; % Transversalitätsbedingung für Umschaltpunkt

end
