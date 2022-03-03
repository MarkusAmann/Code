function res = bcfcn_free_tf(XL,XR,param,p)
nu_tilde = param(1);
delta_t1 = param(2);
delta_t2 = param(3);

X0 = XL(:,1);
Xf = XR(:,end);
X_internalGNB_L = XL(:,end);
X_internalGNB_R = XR(:,1);

% Endzustand
srf = Xf(1);
vf = Xf(2);
drf = Xf(3);
psirf = Xf(4);
l1f = Xf(5);
l2f = Xf(6);
l3f = Xf(7);
l4f = Xf(8);

% Stellgrößen tf
uf = uopt(Xf,p);
axf = uf(1);
kappaf = uf(2);

% linker Randzustand von t1
srL = X_internalGNB_L(1);
vL = X_internalGNB_L(2);
drL = X_internalGNB_L(3);
psirL = X_internalGNB_L(4);
l1L = X_internalGNB_L(5);
l2L = X_internalGNB_L(6);
l3L = X_internalGNB_L(7);
l4L = X_internalGNB_L(8);

% rechter Randzustand von t1
srR = X_internalGNB_R(1);
vR = X_internalGNB_R(2);
drR = X_internalGNB_R(3);
psirR = X_internalGNB_R(4);
l1R = X_internalGNB_R(5);
l2R = X_internalGNB_R(6);
l3R = X_internalGNB_R(7);
l4R = X_internalGNB_R(8);

% Stellgrößen linker und rechter Rand von t1
uL = uopt(X_internalGNB_L,p);
axL = uL(1);
kappaL = uL(2);

uR = uopt(X_internalGNB_R,p);
axR = uR(1);
kappaR = uR(2);

% sf, drf, psirf sind festgelegt, nur vf ist frei
% res = [X0(1) - p.x0(1);...
%     X0(2) - p.x0(2);...
%     X0(3) - p.x0(3);...
%     X0(4) - p.x0(4);...
%     srf - p.sf;...
%     drf - p.drf;...
%     psirf - p.psirf;...
%     l2f];

% Hamiltonfunktion linker und rechter Rand von t1 für freie t1
H_L = 1/2*p.fx*axL^2 + 1/2*p.fy*kappaL^2*vL^4 + delta_t1*(l1L*vL*cos(psirL)/(1-drL*p.kapparef) + ...
    l2L*axL + l3L*vL*sin(psirL) + l4L*vL*(kappaL - p.kapparef*cos(psirL)/(1-drL*p.kapparef)));
H_R = 1/2*p.fx*axR^2 + 1/2*p.fy*kappaR^2*vR^4 + delta_t2*(l1R*vR*cos(psirR)/(1-drR*p.kapparef) + ...
    l2R*axR + l3R*vR*sin(psirR) + l4R*vR*(kappaR - p.kapparef*cos(psirR)/(1-drR*p.kapparef)));
H_L - H_R;

% Hamiltonfunktion tf
Hf = 1/2*p.fx*axf^2 + 1/2*p.fy*kappaf^2*vf^4 + 1*(l1f*vf*cos(psirf)/(1-drf*p.kapparef) + ...
    l2f*axf + l3f*vf*sin(psirf) + l4f*vf*(kappaf - p.kapparef*cos(psirf)/(1-drf*p.kapparef)));

% nur sf ist festgelegt, vf, drf und psirf sind frei
res = [X0(1) - p.x0(1);... % Anfangswerte
    X0(2) - p.x0(2);...
    X0(3) - p.x0(3);...
    X0(4) - p.x0(4);...
    srf - p.sf;... % Endwerte
    l2f;...
    l3f;...
    l4f;...
    X_internalGNB_L(1:6) - X_internalGNB_R(1:6);... % Stetigkeit am internen Randwert
    X_internalGNB_L(8) - X_internalGNB_R(8);...
    X_internalGNB_L(7) - X_internalGNB_R(7) - nu_tilde;... % Transversalitätsbedingung für l3 (l3(t1-) = l3(t1+) + dg/dx|t1*nu_tilde)
    X_internalGNB_L(3) - p.dr1;... % zusätzliche interne Randbedingung
    Hf + 1;... % Transversalitätsbedingung für freie tf
    H_L - H_R]; % Transversalitätsbedingung für freie t1

% sf und psirf sind festgelegt, vf und drf sind frei
% res = [X0(1) - p.x0(1);...
%     X0(2) - p.x0(2);...
%     X0(3) - p.x0(3);...
%     X0(4) - p.x0(4);...
%     srf - p.sf;...
%     psirf - p.psirf;...
%     l2f;...
%     l3f];

% sf und drf sind festgelegt, vf und psirf sind frei, macht keinen Sinn,
% weil psirf nicht frei sein darf
% res = [X0(1) - p.x0(1);... % Anfangswerte
%     X0(2) - p.x0(2);...
%     X0(3) - p.x0(3);...
%     X0(4) - p.x0(4);...
%     srf - p.sf;... % Endwerte
%     drf - p.drf;...
%     l2f;...
%     l4f;...
%     X_internalGNB_L(1:6) - X_internalGNB_R(1:6);... % Stetigkeit am internen Randwert
%     X_internalGNB_L(8) - X_internalGNB_R(8);...
%     X_internalGNB_L(7) - X_internalGNB_R(7) - nu_tilde;... % Transversalitätsbedingung für l3 (l3(t1-) = l3(t1+) + dg/dx|t1*nu_tilde)
%     X_internalGNB_L(3) - p.dr1];... % zusätzliche interne Randbedingung
%     H_L - H_R]; % Transversalitätsbedingung für freie t1

end