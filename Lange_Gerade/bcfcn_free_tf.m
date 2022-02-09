function res = bcfcn_free_tf(X0,Xf,param,p)

sf = Xf(1);
vf = Xf(2);
af = Xf(3);
l1f = Xf(4);
l2f = Xf(5);
l3f = Xf(6);
jf = uopt(Xf,p);
if jf ~= -l3f/p.fj
    aaa=1;
end
% H(tf)+1=0 ohne Zustandsbestrafung im Gütefunktional
% J=tf+int(1/2*fa*a^2+1/2*fj*j^2)
H_tf = 1/2*p.fa*af^2 + 1/2*p.fj*jf^2 + l1f*vf + l2f*af + l3f*jf;

res = [X0(1) - p.x0(1);...
    X0(2) - p.x0(2);...
    X0(3) - p.x0(3);...
    sf - p.sf;...
    l2f;...
    l3f;...
    H_tf + 1];
end