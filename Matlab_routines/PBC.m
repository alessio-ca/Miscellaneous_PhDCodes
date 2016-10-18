    function vec = PBC(vec,Lx,Ly)
    
    
    % Vector should be in the range 0 -> Lx and 0 -> Ly
    % Therefore, we need to apply the following changes if it's not in this range: 
    if (vec(1) > Lx)
        vec(1) = vec(1)-Lx;
    elseif (vec(1) < 0)
        vec(1) = vec(1)+Lx;
    end
    
    if (vec(2) > Ly)
        vec(2) = vec(2)-Ly;
    elseif (vec(2) < 0)
        vec(2) = vec(2)+Ly;
    end
    
end