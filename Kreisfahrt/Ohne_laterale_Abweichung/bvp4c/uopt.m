function u = uopt(X,p)
lambda = X(3:4);
umax = p.umax; umin = p.umin;
u0 = -lambda(2)/p.fx; % unbeschraenkte Stellgroesse
if p.use_umax 
    if u0(1)>umin(1) & u0(1)<umax(1)
        u(1) = u0(1);
    elseif u0(1)<=umin(1)
        u(1) = umin(1);
    else 
        u(1) = umax(1);
    end
else 
    u = u0;
end

end