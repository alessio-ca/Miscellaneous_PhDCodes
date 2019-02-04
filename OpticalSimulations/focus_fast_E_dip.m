function E = focus_fast_E_dip(id,R,Ei)
% E Electric field [V/m]
%
% E = E(ID,R,Ei) calculates the electric field at positions R (Point)
%   for the induced dipole ID usbject to the electric field Ei.
%   E is a ComplexVector.
%
% See also InducedDipole, InducedDipole.Estandard, Point, ComplexVector.


p = id.alpha*Ei;
R = R-id.rd;

X = [R.X,R.X-1e-12,R.X+1e-12,R.X,R.X,R.X,R.X];
Y = [R.Y,R.Y,R.Y,R.Y-1e-12,R.Y+1e-12,R.Y,R.Y];
Z = [R.Z,R.Z,R.Z,R.Z,R.Z,R.Z-1e-12,R.Z+1e-12];

E = zeros(7,3);

vroty = [cos(pi/2) 0 -sin(pi/2);
    0         1 0;
    sin(pi/2) 0 cos(pi/2)];
vrotx = [1 0 0;
    0 cos(-pi/2) sin(-pi/2);
    0 -sin(-pi/2) cos(-pi/2)];

for n=1:7
    E(n,:) = p(1) * focus_fast_E_dip_standard(id,Point(X(n),Y(n),Z(n)).yrotation(-pi/2))*vroty ...
        + p(2) * focus_fast_E_dip_standard(id,Point(X(n),Y(n),Z(n)).xrotation(pi/2))*vrotx ...
        + p(3) * focus_fast_E_dip_standard(id,Point(X(n),Y(n),Z(n)));
end
end