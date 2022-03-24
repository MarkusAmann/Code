function res = bcfcn_free_tf(XL,XR,param,p)
nu_tilde1_t1 = param(1);
nu_tilde2_t1 = param(2);
nu_tilde1_t2 = param(3);
nu_tilde2_t2 = param(4);
delta_t1 = param(5);
delta_t2 = param(6);
delta_t3 = param(7);

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
drf = Xf(3);
psirf = Xf(4);
l1f = Xf(5);
l2f = Xf(6);
l3f = Xf(7);
l4f = Xf(8);

% Stellgrößen tf, Gerade
uf = uopt(Xf,p);
axf = uf(1);
kappaf = uf(2);

%% interner RB t1
% linker Randzustand von t1
sL_1 = X_internal_L_1(1);
vL_1 = X_internal_L_1(2);
drL_1 = X_internal_L_1(3);
psirL_1 = X_internal_L_1(4);
l1L_1 = X_internal_L_1(5);
l2L_1 = X_internal_L_1(6);
l3L_1 = X_internal_L_1(7);
l4L_1 = X_internal_L_1(8);

% rechter Randzustand von t1
sR_1 = X_internal_R_1(1);
vR_1 = X_internal_R_1(2);
drR_1 = X_internal_R_1(3);
psirR_1 = X_internal_R_1(4);
l1R_1 = X_internal_R_1(5);
l2R_1 = X_internal_R_1(6);
l3R_1 = X_internal_R_1(7);
l4R_1 = X_internal_R_1(8);

% Stellgrößen linker und rechter Rand von t1
uL_1 = uopt(X_internal_L_1,p);
axL_1 = uL_1(1);
kappaL_1 = uL_1(2); 

uR_1 = uopt(X_internal_R_1,p);
axR_1 = uR_1(1);
kappaR_1 = uR_1(2); 

%% interner RB t2
% linker Randzustand von t2
sL_2 = X_internal_L_2(1);
vL_2 = X_internal_L_2(2);
drL_2 = X_internal_L_2(3);
psirL_2 = X_internal_L_2(4);
l1L_2 = X_internal_L_2(5);
l2L_2 = X_internal_L_2(6);
l3L_2 = X_internal_L_2(7);
l4L_2 = X_internal_L_2(8);

% rechter Randzustand von t1
sR_2 = X_internal_R_2(1);
vR_2 = X_internal_R_2(2);
drR_2 = X_internal_R_2(3);
psirR_2 = X_internal_R_2(4);
l1R_2 = X_internal_R_2(5);
l2R_2 = X_internal_R_2(6);
l3R_2 = X_internal_R_2(7);
l4R_2 = X_internal_R_2(8);

% Stellgrößen linker und rechter Rand von t1
uL_2 = uopt(X_internal_L_2,p);
axL_2 = uL_2(1);
kappaL_2 = uL_2(2); 

uR_2 = uopt(X_internal_R_2,p);
axR_2 = uR_2(1);
kappaR_2 = uR_2(2); 

switch p.use_dr
    case 1
        % H(tf)+1=0 mit Zustandsbestrafung im Gütefunktional
        % J=tf+int(1/2*fx*ax^2+1/2*fy*kappa^2*v^4+1/2*fr*dr^2)
        H_tf = 1/2*p.fr*drf^2 + 1/2*p.fx*axf^2 + 1/2*p.fy*kappaf^2*vf^4 + ...
            l1f*vf*cos(psirf)/(1-drf*p.kapparef_straight) + l2f*axf + l3f*vf*sin(psirf) + l4f*vf*(kappaf - p.kapparef_straight*cos(psirf)/(1-drf*p.kapparef_straight));
        H_L_1 = 1/2*p.fr*drL_1^2 + 1/2*p.fx*axL_1^2 + 1/2*p.fy*kappaL_1^2*vL_1^4 + ...
            l1L_1*vL_1*cos(psirL_1)/(1-drL_1*p.kapparef_straight) + l2L_1*axL_1 + l3L_1*vL_1*sin(psirL_1) + l4L_1*vL_1*(kappaL_1 - p.kapparef_straight*cos(psirL_1)/(1-drL_1*p.kapparef_straight));
        H_R_1 = 1/2*p.fr*drR_1^2 + 1/2*p.fx*axR_1^2 + 1/2*p.fy*kappaR_1^2*vR_1^4 + ...
            l1R_1*vR_1*cos(psirR_1)/(1-drR_1*p.kapparef_curve) + l2R_1*axR_1 + l3R_1*vR_1*sin(psirR_1) + l4R_1*vR_1*(kappaR_1 - p.kapparef_curve*cos(psirR_1)/(1-drR_1*p.kapparef_curve));
        H_L_2 = 1/2*p.fr*drL_2^2 + 1/2*p.fx*axL_2^2 + 1/2*p.fy*kappaL_2^2*vL_2^4 + ...
            l1L_2*vL_2*cos(psirL_2)/(1-drL_2*p.kapparef_curve) + l2L_2*axL_2 + l3L_2*vL_2*sin(psirL_2) + l4L_2*vL_2*(kappaL_2 - p.kapparef_curve*cos(psirL_2)/(1-drL_2*p.kapparef_curve));
        H_R_2 = 1/2*p.fr*drR_2^2 + 1/2*p.fx*axR_2^2 + 1/2*p.fy*kappaR_2^2*vR_2^4 + ...
            l1R_2*vR_2*cos(psirR_2)/(1-drR_2*p.kapparef_straight) + l2R_2*axR_2 + l3R_2*vR_2*sin(psirR_2) + l4R_2*vR_2*(kappaR_2 - p.kapparef_straight*cos(psirR_2)/(1-drR_2*p.kapparef_straight));
    case 0
        % H(tf)+1=0 ohne Zustandsbestrafung im Gütefunktional
        % J=tf+int(1/2*fx*ax^2+1/2*fy*kappa^2*v^4)
        H_tf = 1/2*p.fx*axf^2 + 1/2*p.fy*kappaf^2*vf^4 + ...
            l1f*vf*cos(psirf)/(1-drf*p.kapparef_straight) + l2f*axf + l3f*vf*sin(psirf) + l4f*vf*(kappaf - p.kapparef_straight*cos(psirf)/(1-drf*p.kapparef_straight));
        H_L_1 = 1/2*p.fx*axL_1^2 + 1/2*p.fy*kappaL_1^2*vL_1^4 + ...
            l1L_1*vL_1*cos(psirL_1)/(1-drL_1*p.kapparef_straight) + l2L_1*axL_1 + l3L_1*vL_1*sin(psirL_1) + l4L_1*vL_1*(kappaL_1 - p.kapparef_straight*cos(psirL_1)/(1-drL_1*p.kapparef_straight));
        H_R_1 = 1/2*p.fx*axR_1^2 + 1/2*p.fy*kappaR_1^2*vR_1^4 + ...
            l1R_1*vR_1*cos(psirR_1)/(1-drR_1*p.kapparef_curve) + l2R_1*axR_1 + l3R_1*vR_1*sin(psirR_1) + l4R_1*vR_1*(kappaR_1 - p.kapparef_curve*cos(psirR_1)/(1-drR_1*p.kapparef_curve));
        H_L_2 = 1/2*p.fx*axL_2^2 + 1/2*p.fy*kappaL_2^2*vL_2^4 + ...
            l1L_2*vL_2*cos(psirL_2)/(1-drL_2*p.kapparef_curve) + l2L_2*axL_2 + l3L_2*vL_2*sin(psirL_2) + l4L_2*vL_2*(kappaL_2 - p.kapparef_curve*cos(psirL_2)/(1-drL_2*p.kapparef_curve));
        H_R_2 = 1/2*p.fx*axR_2^2 + 1/2*p.fy*kappaR_2^2*vR_2^4 + ...
            l1R_2*vR_2*cos(psirR_2)/(1-drR_2*p.kapparef_straight) + l2R_2*axR_2 + l3R_2*vR_2*sin(psirR_2) + l4R_2*vR_2*(kappaR_2 - p.kapparef_straight*cos(psirR_2)/(1-drR_2*p.kapparef_straight));
end

%%
% sf ist festgelegt, vf, drf, psirf sind frei
% res = [X0(1) - p.x0(1);...
%     X0(2) - p.x0(2);...
%     X0(3) - p.x0(3);...
%     X0(4) - p.x0(4);...
%     srf - p.sf;...
%     l2f;...
%     l3f;...
%     l4f;...
%     H_tf + 1;...
%     sL_1 - sR_1;... % Stetigkeit am internen Randpunkt t1
%     vL_1 - vR_1;...
%     drL_1 - drR_1;...
%     psirL_1 - psirR_1;...
%     l1L_1 - l1R_1 - 2*nu_tilde1_t1;... % Transversalitätsbedingung für l1 (l1(t1-) = l1(t1+) + dg/dx|t1*nu_tilde_t1)
%     l2L_1 - l2R_1;...
%     l3L_1 - l3R_1 - 2*nu_tilde2_t1;...
%     l4L_1 - l4R_1;...
%     sL_1 - p.s1;... % zusätzliche interne Randbedingung für t1
%     drL_1 - p.dr1;...
%     H_L_1 - H_R_1 - 0;... % Transversalitätsbedingung für Umschaltpunkt t1
%     sL_2 - sR_2;... % Stetigkeit am internen Randpunkt t2
%     vL_2 - vR_2;...
%     drL_2 - drR_2;...
%     psirL_2 - psirR_2;...
%     l1L_2 - l1R_2 - 2*nu_tilde1_t2;... % Transversalitätsbedingung für l1 (l1(t2-) = l1(t2+) + dg/dx|t2*nu_tilde_t2)
%     l2L_2 - l2R_2;...
%     l3L_2 - l3R_2 - 2*nu_tilde2_t2;...
%     l4L_2 - l4R_2;...
%     sL_2 - p.s2;... % zusätzliche interne Randbedingung für t2
%     drL_2 - p.dr2;...
%     H_L_2 - H_R_2 - 0]; % Transversalitätsbedingung für Umschaltpunkt t2

%%
% sf, drf, psirf sind festgelegt, vf ist frei
res = [X0(1) - p.x0(1);...
    X0(2) - p.x0(2);...
    X0(3) - p.x0(3);...
    X0(4) - p.x0(4);...
    srf - p.sf;...
    drf - p.drf;...
    psirf - p.psirf;...
    l2f;...
    H_tf + 1;...
    sL_1 - sR_1;... % Stetigkeit am internen Randpunkt t1
    vL_1 - vR_1;...
    drL_1 - drR_1;...
    psirL_1 - psirR_1;...
    l1L_1 - l1R_1 - 2*nu_tilde1_t1;... % Transversalitätsbedingung für l1 (l1(t1-) = l1(t1+) + dg/dx|t1*nu_tilde_t1)
    l2L_1 - l2R_1;...
    l3L_1 - l3R_1 - 2*nu_tilde2_t1;...
    l4L_1 - l4R_1;...
    sL_1 - p.s1;... % zusätzliche interne Randbedingung für t1
    drL_1 - p.dr1;...
    H_L_1 - H_R_1 - 0;... % Transversalitätsbedingung für Umschaltpunkt t1
    sL_2 - sR_2;... % Stetigkeit am internen Randpunkt t2
    vL_2 - vR_2;...
    drL_2 - drR_2;...
    psirL_2 - psirR_2;...
    l1L_2 - l1R_2 - 2*nu_tilde1_t2;... % Transversalitätsbedingung für l1 (l1(t2-) = l1(t2+) + dg/dx|t2*nu_tilde_t2)
    l2L_2 - l2R_2;...
    l3L_2 - l3R_2 - 2*nu_tilde2_t2;...
    l4L_2 - l4R_2;...
    sL_2 - p.s2;... % zusätzliche interne Randbedingung für t2
    drL_2 - p.dr2;...
    H_L_2 - H_R_2 - 0]; % Transversalitätsbedingung für Umschaltpunkt t2

%%
% sf, psirf sind festgelegt, vf, drf sind frei
% res = [X0(1) - p.x0(1);...
%     X0(2) - p.x0(2);...
%     X0(3) - p.x0(3);...
%     X0(4) - p.x0(4);...
%     srf - p.sf;...
%     psirf - p.psirf;...
%     l2f;...
%     l3f;
%     H_tf + 1;...
%     sL_1 - sR_1;... % Stetigkeit am internen Randpunkt t1
%     vL_1 - vR_1;...
%     drL_1 - drR_1;...
%     psirL_1 - psirR_1;...
%     l1L_1 - l1R_1 - 2*nu_tilde1_t1;... % Transversalitätsbedingung für l1 (l1(t1-) = l1(t1+) + dg/dx|t1*nu_tilde_t1)
%     l2L_1 - l2R_1;...
%     l3L_1 - l3R_1 - 2*nu_tilde2_t1;...
%     l4L_1 - l4R_1;...
%     sL_1 - p.s1;... % zusätzliche interne Randbedingung für t1
%     drL_1 - p.dr1;...
%     H_L_1 - H_R_1 - 0;... % Transversalitätsbedingung für Umschaltpunkt t1
%     sL_2 - sR_2;... % Stetigkeit am internen Randpunkt t2
%     vL_2 - vR_2;...
%     drL_2 - drR_2;...
%     psirL_2 - psirR_2;...
%     l1L_2 - l1R_2 - 2*nu_tilde1_t2;... % Transversalitätsbedingung für l1 (l1(t2-) = l1(t2+) + dg/dx|t1*nu_tilde_t2)
%     l2L_2 - l2R_2;...
%     l3L_2 - l3R_2 - 2*nu_tilde2_t2;...
%     l4L_2 - l4R_2;...
%     sL_2 - p.s2;... % zusätzliche interne Randbedingung für t2
%     drL_2 - p.dr2;...
%     H_L_2 - H_R_2 - 0]; % Transversalitätsbedingung für Umschaltpunkt t2
end
