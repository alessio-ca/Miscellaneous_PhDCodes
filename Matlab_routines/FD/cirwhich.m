function outp=cirwhich(x,y,z,r,xin,yin,zin,zasp)
% usage: outp=cirwhich(x,y,z,r,xin,yin,zin,zasp)
% Function to find which cirlce a particular point belongs to
% x,y,z,r specify a set of spheres.
% xin,yin,zin specify a point
% zasp is an aspect ratio one z pixel=zasp*one (x,y) pixel.
% IH 7/3/02
j=1;
outp(j)=0;
% This is a quick botch to handle simple x,y data
if size(z,1)<size(x,1)
    z=zeros(size(x,1),1);
end

% This is another quick botch for the circumstance where just a single r value is provided
if size(r,2)==1
    tmp=r;
    r=zeros(size(x,1),1);
    r(:)=tmp;
    
end

%Now go find the points
for i=1:size(x)
    dist2=(x(i)-xin)^2+(y(i)-yin)^2+zasp^2*(z(i)-zin)^2;
    if (dist2<r(i)^2);
        outp(j)=i;
        j=j+1;
    end
end
