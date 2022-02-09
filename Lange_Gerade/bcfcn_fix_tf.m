function res = bcfcn_fix_tf(X0,Xf,p)
sf = Xf(1);
vf = Xf(2);
af = Xf(3);
l1f = Xf(4);
l2f = Xf(5);
l3f = Xf(6);
jf = uopt(Xf,p);

% nur sf ist festgelegt
res = [X0(1) - p.x0(1);...
    X0(2) - p.x0(2);...
    X0(3) - p.x0(3);...
    sf - p.sf;...
    l2f;...
    l3f];
end