function res = bcfcn_fix_tf(X0,Xf,p)
sf = Xf(1);
vf = Xf(2);
l1f = Xf(3);
l2f = Xf(4);
af = uopt(Xf,p);

% nur sf ist festgelegt
res = [X0(1) - p.x0(1);...
    X0(2) - p.x0(2);...
    sf - p.sf;...
    l2f];
end