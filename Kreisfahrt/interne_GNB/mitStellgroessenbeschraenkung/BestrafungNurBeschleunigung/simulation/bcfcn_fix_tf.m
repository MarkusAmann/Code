function res = bcfcn_fix_tf(XL,XR,p)
X0 = XL(:,1);
Xf = XR(:,end);
X_internalGNB_L = XL(:,end);
X_internalGNB_R = XR(:,1);

srf = Xf(1);
vf = Xf(2);
drf = Xf(3);
psirf = Xf(4);
l1f = Xf(5);
l2f = Xf(6);
l3f = Xf(7);
l4f = Xf(8);
uf = uopt(Xf,p);
axf = uf(1);
kappaf = uf(2);

% sf, drf, psirf sind festgelegt, nur vf ist frei
% res = [X0(1) - p.x0(1);...
%     X0(2) - p.x0(2);...
%     X0(3) - p.x0(3);...
%     X0(4) - p.x0(4);...
%     srf - p.sf;...
%     drf - p.drf;...
%     psirf - p.psirf;...
%     l2f];

% nur sf ist festgelegt, vf, drf und psirf sind frei
res = [X0(1) - p.x0(1);...
    X0(2) - p.x0(2);...
    X0(3) - p.x0(3);...
    X0(4) - p.x0(4);...
    srf - p.sf;...
    l2f;...
    l3f;...
    l4f;...
    X_internalGNB_L(1:2) - X_internalGNB_R(1:2);...
    X_internalGNB_L(3) - p.dr1;...
    X_internalGNB_L(3) - X_internalGNB_R(3);...
    X_internalGNB_L(4:6) - X_internalGNB_R(4:6);...
    X_internalGNB_L(8) - X_internalGNB_R(8)];

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
% res = [X0(1) - p.x0(1);...
%     X0(2) - p.x0(2);...
%     X0(3) - p.x0(3);...
%     X0(4) - p.x0(4);...
%     srf - p.sf;...
%     drf - p.drf;...
%     l2f;...
%     l4f];

end