function res = bcfcn_free_tf(XL,XR,param,p)
nu_tilde = param(1);
delta_tf = param(2);

X0 = XL(:,1);
Xf = XR(:,end);
X_internalGNB_L = [XL(1:3,end); XR(4:6,1)];
X_internalGNB_R = [XR(1:3,1); XL(4:6,end)];
% X_internalGNB_L = XL(:,end);
% X_internalGNB_R = XR(:,1);

% Endzustand
srf = Xf(1);
vf = Xf(2);
axf = Xf(3);
l1f = Xf(4);
l2f = Xf(5);
l3f = Xf(6);

% Stellgrößen tf, Kurve
jf = uopt(Xf,p);

% linker Randzustand von t1
sL = X_internalGNB_L(1);
vL = X_internalGNB_L(2);
axL = X_internalGNB_L(3);
l1L = X_internalGNB_L(4);
l2L = X_internalGNB_L(5);
l3L = X_internalGNB_L(6);

% rechter Randzustand von t1
sR = X_internalGNB_R(1);
vR = X_internalGNB_R(2);
axR = X_internalGNB_R(3);
l1R = X_internalGNB_R(4);
l2R = X_internalGNB_R(5);
l3R = X_internalGNB_R(6);

% Stellgrößen linker und rechter Rand von t1
jL = uopt(X_internalGNB_L,p);
jR = uopt(X_internalGNB_R,p);

% H(tf)+1=0 ohne Zustandsbestrafung im Gütefunktional
% J=tf+int(1/2*fa*ax^2+1/2*fy*kappa^2*v^4)
H_tf = 1/2*p.fj*jf^2 + 1/2*p.fa*axf^2 + l1f*vf + l2f*axf + l3f*jf;
H_L = 1/2*p.fj*jL^2 + 1/2*p.fa*axL^2 + l1L*vL + l2L*axL + l3L*jL;
H_R = 1/2*p.fj*jR^2 + 1/2*p.fa*axR^2 + l1R*vR + l2R*axR + l3R*jR;

% sf ist festgelegt, vf, axf sind frei
res = [X0(1) - p.x0(1);...
    X0(2) - p.x0(2);...
    X0(3) - p.x0(3);...
    srf - p.sf;...
    l2f;...
    l3f;...
    H_tf + 1;...
    sL - sR;... % Stetigkeit am internen Randwert
    vL - vR;...
    axL - axR;...
    l1L - l1R - 2*nu_tilde;... % Transversalitätsbedingung für l1 (l1(t1-) = l1(t1+) + dg/dx|t1*nu_tilde)
    l2L - l2R;...
    l3L - l3R;...
    sL - p.s1;... % zusätzliche interne Randbedingung
    ]; 

end
