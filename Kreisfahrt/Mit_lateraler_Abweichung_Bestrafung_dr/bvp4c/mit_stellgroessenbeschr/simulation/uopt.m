function u = uopt(X,p)
lambda = X(5:8);
umax = p.umax; umin = p.umin;
u0 = [-lambda(2)/p.fx; -lambda(4)/(p.fy*X(2)^3)]; % unbeschraenkte Stellgroesse
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