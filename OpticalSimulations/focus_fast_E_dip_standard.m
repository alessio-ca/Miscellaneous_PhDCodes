function E = focus_fast_E_dip_standard(id,R)
% ESTANDARD Electric field for a dipole with p = 1 at the origin and orineted along z [V/m]
%
% E = ESTANDARD(ID,R) calculates the electric field at positions R (Point)
%   for the induced dipole with p = 1 at the origin and orineted along z.
%   E is a ComplexVector.
%
% See also InducedDipole, Point, ComplexVector.


ONES = ones(size(R));
ZEROS = zeros(size(R));

[theta,phi,r] = Transform.Car2Sph(R.X,R.Y,R.Z);
[urVx,urVy,urVz] = Transform.Sph2CarVector(theta,phi,ZEROS,ZEROS,ONES);
[uthVx,uthVy,uthVz] = Transform.Sph2CarVector(theta,phi,ONES,ZEROS,ZEROS);


kr = (2*pi*sqrt(id.er*id.mr)/id.lambda0)*r;
E = (2*pi*sqrt(id.er*id.mr)/id.lambda0)^3/(4*pi*PhysConst.e0*id.er) * exp(1i*kr)./kr ...
    .* ( ...
    2*cos(theta).*(kr.^-2-1i*kr.^-1) * [urVx,urVy,urVz] ...
    + sin(theta).*(kr.^-2-1i*kr.^-1-1) * [uthVx,uthVy,uthVz] ...
    );
end