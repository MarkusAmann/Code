function res = bcfcn_free_tf(XL,XR,param,p)
nu_tilde_t1 = param(1);
nu_tilde_t2 = param(2);
delta_t1 = param(3);
delta_t2 = param(4);
delta_t3 = param(5);

X0 = XL(:,1);
Xf = XR(:,end);
X_internal_L_1 = XL(:,2);
X_internal_L_2 = XL(:,end);
X_internal_R_1 = XR(:,1);
X_internal_R_2 = XR(:,2);

%% Endpunkt
% Endzustand
srf = Xf(1);
vf = Xf(2);
l1f = Xf(3);
l2f = Xf(4);

% Stellgrößen tf, Gerade
axf = uopt(Xf,p);
kappaf = p.kapparef_straight;

%% interner RB t1
% linker Randzustand von t1
sL_1 = X_internal_L_1(1);
vL_1 = X_internal_L_1(2);
l1L_1 = X_internal_L_1(3);
l2L_1 = X_internal_L_1(4);

% rechter Randzustand von t1
sR_1 = X_internal_R_1(1);
vR_1 = X_internal_R_1(2);
l1R_1 = X_internal_R_1(3);
l2R_1 = X_internal_R_1(4);

% Stellgrößen linker und rechter Rand von t1
axL_1 = uopt(X_internal_L_1,p);
kappaL_1 = p.kapparef_straight; 

axR_1 = uopt(X_internal_R_1,p);
kappaR_1 = p.kapparef_curve;

%% interner RB t2
% linker Randzustand von t2
sL_2 = X_internal_L_2(1);
vL_2 = X_internal_L_2(2);
l1L_2 = X_internal_L_2(3);
l2L_2 = X_internal_L_2(4);

% rechter Randzustand von t2
sR_2 = X_internal_R_2(1);
vR_2 = X_internal_R_2(2);
l1R_2 = X_internal_R_2(3);
l2R_2 = X_internal_R_2(4);

% Stellgrößen linker und rechter Rand von t2
axL_2 = uopt(X_internal_L_2,p);
kappaL_2 = p.kapparef_curve;

axR_2 = uopt(X_internal_R_2,p);
kappaR_2 = p.kapparef_straight;

% H(tf)+1=0 ohne Zustandsbestrafung im Gütefunktional
% J=tf+int(1/2*fx*ax^2+1/2*fy*kappa^2*v^4)
H_tf = 1/2*p.fx*axf^2 + 1/2*p.fy*kappaf^2*vf^4 + l1f*vf + l2f*axf;
H_L_1 = 1/2*p.fx*axL_1^2 + 1/2*p.fy*kappaL_1^2*vL_1^4 + l1L_1*vL_1 + l2L_1*axL_1;
H_R_1 = 1/2*p.fx*axR_1^2 + 1/2*p.fy*kappaR_1^2*vR_1^4 + l1R_1*vR_1 + l2R_1*axR_1;
H_L_2 = 1/2*p.fx*axL_2^2 + 1/2*p.fy*kappaL_2^2*vL_2^4 + l1L_2*vL_2 + l2L_2*axL_2;
H_R_2 = 1/2*p.fx*axR_2^2 + 1/2*p.fy*kappaR_2^2*vR_2^4 + l1R_2*vR_2 + l2R_2*axR_2;

% sf ist festgelegt, vf ist frei
res = [X0(1) - p.x0(1);...
    X0(2) - p.x0(2);...
    srf - p.sf;...
    l2f;...
    H_tf + 1;...
    sL_1 - sR_1;... % Stetigkeit am internen Randpunkt t1
    vL_1 - vR_1;...
    l1L_1 - l1R_1 - 2*nu_tilde_t1;... % Transversalitätsbedingung für l1 (l1(t1-) = l1(t1+) + dg/dx|t1*nu_tilde_t1)
    l2L_1 - l2R_1;...
    sL_1 - p.s1;... % zusätzliche interne Randbedingung für t1
    H_L_1 - H_R_1 - 0;... % Transversalitätsbedingung für Umschaltpunkt t1
    sL_2 - sR_2;... % Stetigkeit am internen Randpunkt t2
    vL_2 - vR_2;...
    l1L_2 - l1R_2 - 2*nu_tilde_t2;... % Transversalitätsbedingung für l1 (l1(t2-) = l1(t2+) + dg/dx|t1*nu_tilde_t2)
    l2L_2 - l2R_2;...
    sL_2 - p.s2;... % zusätzliche interne Randbedingung für t2
    H_L_2 - H_R_2 - 0]; % Transversalitätsbedingung für Umschaltpunkt t2
end
