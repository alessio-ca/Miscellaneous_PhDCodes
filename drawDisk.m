function circlemask = drawDisk(mask,coordCenter,radCenter)
iY = size(mask,1);
iX = size(mask,2);
[x,y]=meshgrid(1:iX,1:iY);

circlemask = ((x-coordCenter(1)).^2 + (y-coordCenter(2)).^2 <= radCenter.^2);
end
