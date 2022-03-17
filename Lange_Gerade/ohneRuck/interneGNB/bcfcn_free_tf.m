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
l1f = Xf(3);
l2f = Xf(4);

% Stellgrößen tf
uf = uopt(Xf,p);
af = uf;

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
uL = uopt(X_internalGNB_L,p);
aL = uL;

uR = uopt(X_internalGNB_R,p);
aR = uR;

% Hamiltonfunktion linker und rechter Rand von t1 für freie t1 
H_L = 1*(1/2*p.fa*aL^2 + l1L*vL + l2L*aL);
H_R = 1*(1/2*p.fa*aR^2 + l1R*vR + l2R*aR);
% H(tf)+1=0 ohne Zustandsbestrafung im Gütefunktional
% J=tf+int(1/2*fa*a^2+1/2*fj*j^2)
H_tf = 1*(1/2*p.fa*af^2 + l1f*vf + l2f*af);


res = [X0(1) - p.x0(1);...
    X0(2) - p.x0(2);...
    sf - p.sf;...
    l2f;...
    H_tf + 1;...
    sL - sR;... % Stetigkeit am internen Randwert
    vL - vR;
    l1L - l1R - 2*nu_tilde;... % Transversalitätsbedingung für l1 (l1(t1-) = l1(t1+) + dg/dx|t1*nu_tilde)
    l2L - l2R;...
    sL - p.s1;... % zusätzliche interne Randbedingung
    H_L - H_R - 1    % Transversalitätsbedingung für Umschaltpunkt
    ];
end