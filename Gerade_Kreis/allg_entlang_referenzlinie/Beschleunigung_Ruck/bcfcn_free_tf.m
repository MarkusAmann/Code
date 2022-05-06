function res = bcfcn_free_tf(XL,XR,param,p)
nu_tilde = param(1);
delta_t1 = param(2);
delta_t2 = param(3);

X0 = XL(:,1);
Xf = XR(:,end);
X_internalGNB_L = [XL(1:5,end); XR(6:10,1)];
X_internalGNB_R = [XR(1:5,1); XL(6:10,end)];
% X_internalGNB_L = XL(:,end);
% X_internalGNB_R = XR(:,1);

% Endzustand
srf = Xf(1);
vf = Xf(2);
af = Xf(3);
drf = Xf(4);
psirf = Xf(5);
l1f = Xf(6);
l2f = Xf(7);
l3f = Xf(8);
l4f = Xf(9);
l5f = Xf(10);

% Stellgrößen tf, Kurve
uf = uopt(Xf,p);
jf = uf(1);
kappaf = uf(2);

% linker Randzustand von t1
sL = X_internalGNB_L(1);
vL = X_internalGNB_L(2);
aL = X_internalGNB_L(3);
drL = X_internalGNB_L(4);
psirL = X_internalGNB_L(5);
l1L = X_internalGNB_L(6);
l2L = X_internalGNB_L(7);
l3L = X_internalGNB_L(8);
l4L = X_internalGNB_L(9);
l5L = X_internalGNB_L(10);

% rechter Randzustand von t1
sR = X_internalGNB_R(1);
vR = X_internalGNB_R(2);
aR = X_internalGNB_R(3);
drR = X_internalGNB_R(4);
psirR = X_internalGNB_R(5);
l1R = X_internalGNB_R(6);
l2R = X_internalGNB_R(7);
l3R = X_internalGNB_R(8);
l4R = X_internalGNB_R(9);
l5R = X_internalGNB_R(10);

% Stellgrößen linker und rechter Rand von t1
uL = uopt(X_internalGNB_L,p);
jL = uL(1);
kappaL = uL(2); 

uR = uopt(X_internalGNB_R,p);
jR = uR(1);
kappaR = uR(2); 

switch p.use_dr
    case 1
        % H(tf)+1=0 mit Zustandsbestrafung im Gütefunktional
        % J=tf+int(1/2*fa*ax^2+1/2*fy*kappa^2*v^4+1/2*fr*dr^2)
        H_tf = 1/2*p.fr*drf^2 + 1/2*p.fj*jf^2 + 1/2*p.fa*af^2 + 1/2*p.fy*kappaf^2*vf^4 + ...
            l1f*vf*cos(psirf)/(1-drf*p.kapparef_curve) + l2f*af + l3f*jf + l4f*vf*sin(psirf) + l5f*vf*(kappaf - p.kapparef_curve*cos(psirf)/(1-drf*p.kapparef_curve));
        H_L = 1/2*p.fr*drL^2 + 1/2*p.fj*jL^2 + 1/2*p.fa*aL^2 + 1/2*p.fy*kappaL^2*vL^4 + ...
            l1L*vL*cos(psirL)/(1-drL*p.kapparef_straight) + l2L*aL + l3L*jL + l4L*vL*sin(psirL) + l5L*vL*(kappaL - p.kapparef_straight*cos(psirL)/(1-drL*p.kapparef_straight));
        H_R = 1/2*p.fr*drR^2 + 1/2*p.fj*jR^2 + 1/2*p.fa*aR^2 + 1/2*p.fy*kappaR^2*vR^4 + ...
            l1R*vR*cos(psirR)/(1-drR*p.kapparef_curve) + l2R*aR + l3R*jR + l4R*vR*sin(psirR) + l5R*vR*(kappaR - p.kapparef_curve*cos(psirR)/(1-drR*p.kapparef_curve));
    case 0
        % H(tf)+1=0 ohne Zustandsbestrafung im Gütefunktional
        % J=tf+int(1/2*fa*ax^2+1/2*fy*kappa^2*v^4)
        H_tf = 1/2*p.fj*jf^2 + 1/2*p.fa*af^2 + 1/2*p.fy*kappaf^2*vf^4 + ...
            l1f*vf*cos(psirf)/(1-drf*p.kapparef_curve) + l2f*af + l3f*jf + l4f*vf*sin(psirf) + l5f*vf*(kappaf - p.kapparef_curve*cos(psirf)/(1-drf*p.kapparef_curve));
        H_L = 1/2*p.fj*jL^2 + 1/2*p.fa*aL^2 + 1/2*p.fy*kappaL^2*vL^4 + ...
            l1L*vL*cos(psirL)/(1-drL*p.kapparef_straight) + l2L*aL + l3L*jL + l4L*vL*sin(psirL) + l5L*vL*(kappaL - p.kapparef_straight*cos(psirL)/(1-drL*p.kapparef_straight));
        H_R = 1/2*p.fj*jR^2 + 1/2*p.fa*aR^2 + 1/2*p.fy*kappaR^2*vR^4 + ...
            l1R*vR*cos(psirR)/(1-drR*p.kapparef_curve) + l2R*aR + l3R*jR + l4R*vR*sin(psirR) + l5R*vR*(kappaR - p.kapparef_curve*cos(psirR)/(1-drR*p.kapparef_curve));
end

% sf ist festgelegt, vf, af, drf, psirf sind frei
% res = [X0(1) - p.x0(1);...
%     X0(2) - p.x0(2);...
%     X0(3) - p.x0(3);...
%     X0(4) - p.x0(4);...
%     X0(5) - p.x0(5);...
%     srf - p.sf;...
%     l2f;...
%     l3f;
%     l4f;...
%     l5f;...
%     H_tf + 1;...
%     sL - sR;... % Stetigkeit am internen Randwert
%     vL - vR;...
%     drL - drR;...
%     aL - aR;...
%     psirL - psirR;...
%     l1L - l1R - 2*nu_tilde;... % Transversalitätsbedingung für l1 (l1(t1-) = l1(t1+) + dg/dx|t1*nu_tilde)
%     l2L - l2R;...
%     l3L - l3R;...
%     l4L - l4R;...
%     l5L - l5R;...
%     sL - p.s1;... % zusätzliche interne Randbedingung
%     H_L - H_R + 0]; % Transversalitätsbedingung für Umschaltpunkt

% sf, drf, psirf sind festgelegt, vf, af sind frei
% res = [X0(1) - p.x0(1);...
%     X0(2) - p.x0(2);...
%     X0(3) - p.x0(3);...
%     X0(4) - p.x0(4);...
%     X0(5) - p.x0(5);...
%     srf - p.sf;...
%     drf - p.drf;...
%     psirf - p.psirf;...
%     l2f;...
%     l3f;
%     H_tf + 1;...
%     sL - sR;... % Stetigkeit am internen Randwert
%     vL - vR;...
%     drL - drR;...
%     aL - aR;...
%     psirL - psirR;...
%     l1L - l1R - 2*nu_tilde;... % Transversalitätsbedingung für l1 (l1(t1-) = l1(t1+) + dg/dx|t1*nu_tilde)
%     l2L - l2R;...
%     l3L - l3R;...
%     l4L - l4R;...
%     l5L - l5R;...
%     sL - p.s1;... % zusätzliche interne Randbedingung
%     H_L - H_R + 0]; % Transversalitätsbedingung für Umschaltpunkt

% sf, psirf sind festgelegt, vf, af, drf sind frei
res = [X0(1) - p.x0(1);...
    X0(2) - p.x0(2);...
    X0(3) - p.x0(3);...
    X0(4) - p.x0(4);...
    X0(5) - p.x0(5);...
    srf - p.sf;...
    psirf - p.psirf;...
    l2f;...
    l3f;
    l4f;...
    H_tf + 1;...
    sL - sR;... % Stetigkeit am internen Randwert
    vL - vR;...
    drL - drR;...
    aL - aR;...
    psirL - psirR;...
    l1L - l1R - 2*nu_tilde;... % Transversalitätsbedingung für l1 (l1(t1-) = l1(t1+) + dg/dx|t1*nu_tilde)
    l2L - l2R;...
    l3L - l3R;...
    l4L - l4R;...
    l5L - l5R;...
    sL - p.s1;... % zusätzliche interne Randbedingung
    H_L - H_R + 0]; % Transversalitätsbedingung für Umschaltpunkt

% sf, af, psirf sind festgelegt, vf, drf sind frei
% res = [X0(1) - p.x0(1);...
%     X0(2) - p.x0(2);...
%     X0(3) - p.x0(3);...
%     X0(4) - p.x0(4);...
%     X0(5) - p.x0(5);...
%     srf - p.sf;...
%     af - 0;...
%     psirf - p.psirf;...
%     l2f;...
%     l4f;...
%     H_tf + 1;...
%     sL - sR;... % Stetigkeit am internen Randwert
%     vL - vR;...
%     drL - drR;...
%     aL - aR;...
%     psirL - psirR;...
%     l1L - l1R - 2*nu_tilde;... % Transversalitätsbedingung für l1 (l1(t1-) = l1(t1+) + dg/dx|t1*nu_tilde)
%     l2L - l2R;...
%     l3L - l3R;...
%     l4L - l4R;...
%     l5L - l5R;...
%     sL - p.s1;... % zusätzliche interne Randbedingung
%     H_L - H_R + 0]; % Transversalitätsbedingung für Umschaltpunkt
end
