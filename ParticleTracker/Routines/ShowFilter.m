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
% MODIFIED: Alessio Caciagli, University of Cambridge, 21/12/2016

%Create a representatopm of the filter in a full-frame image
filter_TOT=zeros(size(image));
if round(cnt(1,1))+floor(size(filter,2)/2)<length(filter) || round(cnt(1,2))+floor(size(filter,1)/2)<length(filter)
    %Centroid is next to edges
    
    %Create square of size(filter) centered in cnt, with pixel coordinates
    squareX=(round(cnt(1,2))-length(filter)+1:round(cnt(1,2)))+floor(size(filter,1)/2);
    squareY=(round(cnt(1,1))-length(filter)+1:round(cnt(1,1)))+floor(size(filter,2)/2);
    
    %Calculate the offset from the edge (minimum between X & Y)
    offsetX=2-(round(cnt(1,2))-length(filter)+1+floor(size(filter,1)/2));
    offsetY=2-(round(cnt(1,1))-length(filter)+1+floor(size(filter,2)/2));

    if offsetX>offsetY
        offsetY=1;
    else
        offsetX=1;
    end
    
    %Cut the filter to the offset value
    squareX = squareX(offsetX:end);
    squareY = squareY(offsetY:end);
   
    filter_TOT(squareX,squareY) = filter(offsetX:end,offsetY:end);
else
    square=[(round(cnt(1,2))-length(filter)+1:round(cnt(1,2)))+floor(size(filter,1)/2);(round(cnt(1,1))-length(filter)+1:round(cnt(1,1)))+floor(size(filter,2)/2)];
    filter_TOT(square(1,:),square(2,:)) = filter;
end

line_x = [round(cnt(1,1))-length(filter),round(cnt(1,1))+length(filter)];
line_y = [round(cnt(1,2))-length(filter),round(cnt(1,2))+length(filter)];

line_x(line_x < 0) = 1;
line_x(line_x > size(image,2))=size(image,2);

line_y(line_y < 0) = 1;
line_y(line_y > size(image,1))=size(image,1);



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
clf
imshow(image)
hold on
h = imshow(filter_TOT);
set(h, 'AlphaData',1-mean(image(:)))
line([line_x(1),line_x(2)],[line_y(1),line_y(1)],'Color','w')
line([line_x(2),line_x(2)],[line_y(1),line_y(2)],'Color','w')
line([line_x(2),line_x(1)],[line_y(2),line_y(2)],'Color','w')
line([line_x(1),line_x(1)],[line_y(2),line_y(1)],'Color','w')
hold off
zoom reset
set(gca,'XLim',X_Lim,'YLim',Y_Lim);

end

    