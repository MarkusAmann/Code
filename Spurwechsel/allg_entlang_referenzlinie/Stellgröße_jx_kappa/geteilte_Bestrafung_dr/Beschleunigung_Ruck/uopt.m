function u = uopt(X,p)
sr = X(1);
v = X(2);
a = X(3);
dr = X(4);
psir = X(5);
l1 = X(6);
l2 = X(7);
l3 = X(8);
l4 = X(9);
l5 = X(10);
umax = p.umax; umin = p.umin;
u0 = [-l3/p.fjx; -l5/(p.fay*v^3)]; % unbeschraenkte Stellgroesse
if p.use_umax 
    if u0(1)>umin(1) & u0(1)<umax(1)
        u(1) = u0(1);
    elseif u0(1)<=umin(1)
        u(1) = umin(1);
    else 
        u(1) = umax(1);
    end
    
    if u0(2)>umin(2) & u0(2)<umax(2)
        u(2) = u0(2);
    elseif u0(2)<=umin(2)
        u(2) = umin(2);
    else 
        u(2) = umax(2);
    end
    u = u.';
else 
    u = u0;
end

end