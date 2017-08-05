function [ur,uth,r_V]=focus_fast_E_Vector_Sph(R)
    ONES = ones(size(R));
    ZEROS = zeros(size(R));

    [theta,phi,r] = Transform.Car2Sph(R.X,R.Y,R.Z);

    [urVx,urVy,urVz] = Transform.Sph2CarVector(theta,phi,ZEROS,ZEROS,ONES);

    [uthVx,uthVy,uthVz] = Transform.Sph2CarVector(theta,phi,ONES,ZEROS,ZEROS);
    
    r_V=[theta,r];
    ur=[urVx,urVy,urVz];
    uth=[uthVx,uthVy,uthVz];
    
end