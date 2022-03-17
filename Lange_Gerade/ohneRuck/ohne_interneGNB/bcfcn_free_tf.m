function res = bcfcn_free_tf(X0,Xf,param,p)
delta_tf = param;
sf = Xf(1);
vf = Xf(2);
l1f = Xf(3);
l2f = Xf(4);
af = uopt(Xf,p);

% H(tf)+1=0 ohne Zustandsbestrafung im Gütefunktional
% J=tf+int(1/2*fa*a^2+1/2*fj*j^2)
H_tf = delta_tf*(1/2*p.fa*af^2 + l1f*vf + l2f*af);

res = [X0(1) - p.x0(1);...
    X0(2) - p.x0(2);...
    sf - p.sf;...
    delta_tf*l2f;...
%     vf - p.vf;
    H_tf + delta_tf];
end