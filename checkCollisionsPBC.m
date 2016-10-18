    function check = checkCollisionsPBC(coords,rTrial,iMove,nPart,Lx,Ly)
    
    check = true;
    
    for i=1:nPart
       % Compare only to particle that are not the one being moved
       if (i ~= iMove) 
    
           % Get the displacement vector between particle i and the suggested
           % trial position
           dr = coords(:,i)-rTrial;
           
           % Apply periodic boundary conditions
           drt = distPBC2D(dr,Lx,Ly);
           
           % The magnitude of this vector can be found by a dot product
           drmag2 = drt*drt';
           
           if (drmag2 < 1)
               % The particles collide, no point in checking any further
               check=false;
               break;
           end
       end
    end
    
    end
    

