function res = bcfcn_free_tf(XL,XR,param,p)
nu_tilde1 = param(1);
nu_tilde2 = param(2);
delta_t1 = param(3);
delta_t2 = param(4);

X0 = XL(:,1);
Xf = XR(:,end);
X_internalGNB_L = [XL(1:2,end); XR(3:4,1)];
X_internalGNB_R = [XR(1:2,1); XL(3:4,end)];

% Endzustand
srf = Xf(1);
vf = Xf(2);
l1f = Xf(3);
l2f = Xf(4);

% Stellgrößen tf, Kurve
axf = uopt(Xf,p);
kappaf = p.kapparef_curve;

% linker Randzustand von t1
sL = X_internalGNB_L(1);
vL = X_internalGNB_L(2);
l1L = X_internalGNB_L(3);
l2L = X_internalGNB_L(4);

% rechter Randzustand von t1
sR = X_internalGNB_R(1);
vR = X_internalGNB_R(2);
l1R = X_internalGNB_R(3);
l2R = X_internalGNB_R(4);

% Stellgrößen linker und rechter Rand von t1
axL = uopt(X_internalGNB_L,p);
kappaL = p.kapparef_straight; % Gerade

axR = uopt(X_internalGNB_R,p);
kappaR = p.kapparef_curve; % Kurve

% H(tf)+1=0 ohne Zustandsbestrafung im Gütefunktional
% J=tf+int(1/2*fx*ax^2+1/2*fy*kappa^2*v^4)
H_tf = 1/2*p.fx*axf^2 + 1/2*p.fy*kappaf^2*vf^4 + l1f*vf + l2f*axf;
H_L = 1/2*p.fx*axL^2 + 1/2*p.fy*kappaL^2*vL^4 + l1L*vL + l2L*axL;
H_R = 1/2*p.fx*axR^2 + 1/2*p.fy*kappaR^2*vR^4 + l1R*vR + l2R*axR;

% sf ist festgelegt, vf ist frei
res = [X0(1) - p.x0(1);...
    X0(2) - p.x0(2);...
    srf - p.sf;...
    l2f;...
    H_tf + 1;...
    sL - sR;... % Stetigkeit am internen Randwert
    vL - vR;...
    l1L - l1R - 2*nu_tilde1;... % Transversalitätsbedingung für l1 (l1(t1-) = l1(t1+) + dg/dx|t1*nu_tilde)
    l2L - l2R - 2*nu_tilde2;...
    sL - p.s1;... % zusätzliche interne Randbedingung
    vR - vf;...
    H_L - H_R + 0]; % Transversalitätsbedingung für Umschaltpunkt

end
