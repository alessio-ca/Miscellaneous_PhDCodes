function outp=cirpaint2(x,y,z,r,orgim,zorg,zsp,flag);
% Usage outp=cirpaint2(x,y,z,r,orgim,zorg,zsp,flag);
% Function to plot images overlaid with circles or spheres specified by x,y,z,r
% If no image, orgim, is supplied then a blank image is created
% zorg if a 3D set of spheres is specified then this sets the number of the first slice 
% in the stack for which spheres are drawn.
% Could do with putting in a paint style flag 1 for fill, 0 for outline
% IH 7/3/02

% 19/11/02 - Put in some mods to work with a particular dataset Exper 12 p42 IH
%
% 23/11/02 - Further mod puts in an aspect ratio in the call and code to measure the size of the
%            input file rather than hardcode it.
%          - Implemented the flag option, 1 plots solid spheres, 0 outlines
% 4/3/03   - Extending the flag option, Jasna would like black outlines!
%            Now handles white outlines(0),white circles(1), black outlines(2), black circles(3), no overlay(4) 

%xsize=256; %512;
%ysize=256; %512; Now in call
%zsize=59; %1;
%zsp=1.333; %2.67; %aspect ratio

xsize=size(orgim,1);
ysize=size(orgim,2);
zsize=size(orgim,3);

outp=zeros(xsize,ysize,zsize);
% This is a quick botch to handle simple x,y data
if size(z,1)<size(x,1)
    z=ones(size(x,1),1);
end

for m=1:size(x,1) % Loop over whole list of shapes.
        xs=x(m);
        ys=y(m);
        zs=z(m);
        for p=zorg:zorg+zsize-1 %loop over slices
            if abs(zsp*(zs-p))<r(m)
                %Draw a circle in a plane
                ra=sqrt(r(m)^2-zsp^2*(zs-p)^2);
                for n=-pi/2:(1/ra):pi/2 %Can put a step size in here to make lines outlined circles fill a little better.
                    xs1=round(xs+ra*cos(n));
                    xs2=round(xs-ra*cos(n));
                    ys1=round(ys+ra*sin(n));
                    if xs1<1 xs1=1; end;
                    if xs1>xsize xs1=xsize; end;
                    if xs2<1 xs2=1; end;
                    if xs2>xsize xs2=xsize; end;
                    if ys1>ysize ys1=ysize; end;
                    if ys1<1 ys1=1; end;
          %This is where any flag options are enacted
                    if flag==1
                        outp(xs2:xs1,ys1,p-zorg+1)=255; %Plot as solid circles
                    elseif flag==0
                        outp(xs1,ys1,p-zorg+1)=255; %Plot as white outlines
                        outp(xs2,ys1,p-zorg+1)=255; %Plot as white outlines
                    elseif flag==2
                        outp(xs1,ys1,p-zorg+1)=-255; %Plot as black outlines
                        outp(xs2,ys1,p-zorg+1)=-255; %Plot as black outlines
                    elseif flag==3
                        outp(xs2:xs1,ys1,p-zorg+1)=-255; %Plot as solid black circles
                    elseif flag==4
                        outp(xs2:xs1,ys1,p-zorg+1)=0; %No overlay
                    end
                end
            end
        end
end
outp=uint8(double(outp)+double(orgim));