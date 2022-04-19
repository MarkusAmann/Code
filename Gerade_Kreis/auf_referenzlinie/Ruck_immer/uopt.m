function u = uopt(X,p)
lambda = X(4:6);
umax = p.umax; umin = p.umin;
u0 = -lambda(3)/p.fx; % unbeschraenkte Stellgroesse
if p.use_umax 
    if u0>umin & u0<umax
        u = u0;
    elseif u0<=umin
        u = umin;
    else 
        u = umax;
    end
else 
    u = u0;
end

end