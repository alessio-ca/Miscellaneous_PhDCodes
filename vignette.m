function masked_image = vignette(image,radius,value)
% out=cntrd(im,mx,sz,interactive)
% 
% PURPOSE:  
%           Mask a vignetted image by setting the locations outside of a
%           radius equal to value. The default value is 0.
%           
% 
% INPUT:
% image:    The original image
% radius:   The radius of the desired vignette in pixels. If the radius is
%           -1, just return the original image.
% value:    Set to value. OPTIONAL. Default to 0.
%
% NOTES:
%
% OUTPUT:  A masked image
%
% CREATED: Eric M. Furst, University of Delaware, July 23, 2013
%  Modifications:


if nargin==2
   value=0; 
end

if radius<0

    masked_image = image;

else    
    % Make mask
    % Get the size of image
    % 
    [ix, iy] = size(image);

    % cx and cy are the center coordinates of the circle
    cx=ix/2;
    cy=iy/2;

    % 
    [x,y]=meshgrid(-(cx-1):(ix-cx),-(cy-1):(iy-cy));
    c_mask=((x.^2+y.^2)<=radius^2)';
    mask_outside = (c_mask-1)*(-value) + c_mask;

    masked_image = image.*c_mask + mask_outside;
end