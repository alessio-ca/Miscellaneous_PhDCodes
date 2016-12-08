function ShowFilter(image,filter,cnt,colour)
%
% USAGE:    ShowFilter(image,filter,cnt)
% PURPOSE:  Overlays a filter on an image to show the tracked particle location.
%           N.B! Linked to the usage of AnalyseTrack particle tracking
%
%
% INPUT:
% image:    current image of interest
% filter:   filter which has been used to cross-correlate the particle(s)
%           position(s)
% cnt:      position of the tracked centroid (in pixel)
% colour:   colour of the overlayed filter (must be inserted in the form of
%           an rgb vector)
%
% NOTES:
%
% OUTPUT:
%
% CREATED: Alessio Caciagli, University of Cambridge, 06/04/2016

filter_TOT=filter(1,1)*ones(size(image));
if size(cnt,1)==1
    if round(cnt(1,1))+floor(size(filter,2)/2)<length(filter) || round(cnt(1,2))+floor(size(filter,1)/2)<length(filter)
        square1=[(round(cnt(1,2))-length(filter)+1:round(cnt(1,2)))+floor(size(filter,1)/2);(round(cnt(1,1))-length(filter)+1:round(cnt(1,1)))+floor(size(filter,2)/2)];
        offset=2-min(round(cnt(1,1))-length(filter)+1+floor(size(filter,2)/2),round(cnt(1,2))-length(filter)+1+floor(size(filter,1)/2));
        square1=square1(:,offset:end);
        filter_TOT(square1(1,:),square1(2,:)) = filter(offset:end,offset:end);
    else
        square1=[(round(cnt(1,2))-length(filter)+1:round(cnt(1,2)))+floor(size(filter,1)/2);(round(cnt(1,1))-length(filter)+1:round(cnt(1,1)))+floor(size(filter,2)/2)];
        filter_TOT(square1(1,:),square1(2,:)) = filter;
    end
elseif size(cnt,1)==2
    if round(cnt(1,1))<length(filter) || round(cnt(1,2))<length(filter)
        filter_TOT=zeros(size(filter));
    else
        square1=[(round(cnt(1,2))-length(filter)+1:round(cnt(1,2)))+floor(size(filter,1)/2);(round(cnt(1,1))-length(filter)+1:round(cnt(1,1)))+floor(size(filter,2)/2)];
        filter_TOT(square1(1,:),square1(2,:)) = filter;
    end
    if round(cnt(2,1))>=length(filter) && round(cnt(2,2))>=length(filter)
        square2=[(round(cnt(2,2))-length(filter)+1:round(cnt(2,2)))+floor(size(filter,1)/2);(round(cnt(2,1))-length(filter)+1:round(cnt(2,1)))+floor(size(filter,2)/2)];
        filter_TOT(square2(1,:),square2(2,:)) = max(filter,filter_TOT(square2(1,:),square2(2,:)));
    end
else
    error('Error in the dimensions of cnt.');
end
filter_TOT(filter_TOT<=filter(1,1))=0;
r = colour(1);
g = colour(2);
b = colour(3);
rgbImage = cat(3,filter_TOT, filter_TOT, filter_TOT); % Convert the gray image into RGB so we can get it into hsv space.
hsv = rgb2hsv(rgbImage);
hsvOfThisPixel = rgb2hsv([r g b]); % Make the hue plane of the hsv image this same hue
hsv(:, :, 1) = hsvOfThisPixel(1);  % Make the saturation 1 so we'll see color (originally it was zero which meant that the image would be gray).
hsv(:, :, 2) = 1;
filter_TOT = uint8(255*hsv2rgb(hsv)); % Convert back to rgb space.

%Show an overlay between image and filter mask
X_Lim=get(gca,'XLim'); %Commands to keep the zoom in imshow during filter creation
Y_Lim=get(gca,'YLim');
figure(1)
imshow(imadjust(image))
hold on
h = imshow(filter_TOT);
set(h, 'AlphaData', 0.5)
hold off
zoom reset
set(gca,'XLim',X_Lim,'YLim',Y_Lim);

end

    