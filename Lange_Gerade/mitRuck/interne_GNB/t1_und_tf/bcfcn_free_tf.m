function res = bcfcn_free_tf(XL,XR,param,p)
nu_tilde = param(1);
delta_t1 = param(2);
delta_t2 = param(3);

X0 = XL(:,1);
Xf = XR(:,end);
X_internalGNB_L = XL(:,end);
X_internalGNB_R = XR(:,1);

% Endzustand
sf = Xf(1);
vf = Xf(2);
af = Xf(3);
l1f = Xf(4);
l2f = Xf(5);
l3f = Xf(6);

% Stellgrößen tf
uf = uopt(Xf,p);
jf = uf;

% linker Randzustand von t1
sL = X_internalGNB_L(1);
vL = X_internalGNB_L(2);
aL = X_internalGNB_L(3);
l1L = X_internalGNB_L(4);
l2L = X_internalGNB_L(5);
l3L = X_internalGNB_L(6);

% rechter Randzustand von t1
sR = X_internalGNB_R(1);
vR = X_internalGNB_R(2);
aR = X_internalGNB_R(3);
l1R = X_internalGNB_R(4);
l2R = X_internalGNB_R(5);
l3R = X_internalGNB_R(6);

% Stellgrößen linker und rechter Rand von t1
uL = uopt(X_internalGNB_L,p);
jL = uL;

uR = uopt(X_internalGNB_R,p);
jR = uR;

% Hamiltonfunktion linker und rechter Rand von t1 für freie t1 
H_L = 1/2*p.fa*aL^2 + 1/2*p.fj*jL^2 + l1L*vL + l2L*aL + l3L*jL;
H_R = 1/2*p.fa*aR^2 + 1/2*p.fj*jR^2 + l1R*vR + l2R*aR + l3R*jR;
% H(tf)+1=0 ohne Zustandsbestrafung im Gütefunktional
% J=tf+int(1/2*fa*a^2+1/2*fj*j^2)
H_tf = 1/2*p.fa*af^2 + 1/2*p.fj*jf^2 + l1f*vf + l2f*af + l3f*jf;

% nur sf ist festgelegt
res = [X0(1) - p.x0(1);...
    X0(2) - p.x0(2);...
    X0(3) - p.x0(3);...
    sf - p.sf;...
    l2f;...
    l3f;...
    X_internalGNB_L(1:3) - X_internalGNB_R(1:3);... % Stetigkeit am internen Randwert
    X_internalGNB_L(4) - X_internalGNB_R(4) - 2*nu_tilde;... % Transversalitätsbedingung für l1 (l1(t1-) = l1(t1+) + dg/dx|t1*nu_tilde)
    X_internalGNB_L(5:6) - X_internalGNB_R(5:6);...
    X_internalGNB_L(1) - p.s1;... % zusätzliche interne Randbedingung
    H_L - H_R - 1; % Transversalitätsbedingung für Umschaltpunkt
    H_tf + 1]; % Transversalitätsbedingung für Endpunkt
end