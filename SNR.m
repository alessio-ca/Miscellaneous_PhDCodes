function [signal,noise]=SNR(image,feature,maskradius)
% [signal,noise]=SNR(image,feature,maskradius)
% 
% PURPOSE: Calculate the SNR for a single image. We calculate the
% distribution of pixels that consist of our signal and the noise. The
% signal are the pixels associated with tracked particles. The noise
% consists of all pixels that are not signal. We do this by masking the
% tracked particles (assigning their region with a negative value). Note
% that the signal area is two pixels smaller than the feature size and the
% noise area is two TIMES feature size larger. This is to prevent the edges
% of particles from contributing to the signal or noise.
%
% Load the image with command: I = double(imread('test.tif'));
%   
%           
% 
% INPUT:
%
% image: The image to analyze
% 
% feature:  Expected size of the particle (diameter) that constitutes the
% signal region. A value 5 seems to work well.
%
% OPTIONAL INPUT: 
%
% maskradius: If an image is vignetted, this is the radius
% of the image to use in calculations. Otherwise, this defaults to a value
% 450. A maskradius equal to -1 will ignore the vignette.
%
% NOTES:
%
% OUTPUT:  
% signal: A vector of pixel values of the tracked particles
% noise: A vector of pixel values of the background image excluding tacked
% particles
%
% CREATED: Eric M. Furst, University of Delaware, July 25, 2013
%  Modifications:

if mod(feature,2) == 0 || feature < 3
    warning('Feature size must be an odd value >= 3.');
    signal=[];
    noise=[];
    return;
end    

% Set the maskradius to a default value
if nargin==2
   maskradius=450; 
end


% read in file
% image = double(imread(filename));

% Bandpass filter
imagebp = bpass(image,1,feature);

% Find locations of the brightest pixels Need to change the 'th'
% (threshold) argument to something that is specified or determined.
% A rough guide is 0.6*max(imagebp(:))
th = 0.5*max(imagebp(:));
pk = pkfnd(imagebp, th, feature);
disp(['Maximum imagebp value is ', num2str(max(imagebp(:)))]);
disp(['Using values > ', num2str(th)]);
disp(['located ' num2str(size(pk,1)) ' particles.']);

% Now we have particle positions. Need to mask the image at these
% locations. Make the particle mask. This is an array with -1 for
% values within a radius of half the feature size.
sz = 2*feature+1;
r = (sz-1)/2;
[x, y]=meshgrid(-(r):(sz-r-1),-(r):(sz-r-1));
particlemask=((x.^2+y.^2)<=r^2);
particlemask=((x.^2+y.^2)>r^2)-particlemask;

% go through list of particle locations
% make a circle of radius feature at each location with value -1
imagenoise = image;
for i=1:size(pk,1);
    x = pk(i,1);
    y = pk(i,2);
    imagenoise((y-r):(y+r),(x-r):(x+r))=...
        image((y-r):(y+r),(x-r):(x+r)).*particlemask;
end 

% Make the mask for the signal
sz = feature-2;
r = (sz-1)/2;
[x, y]=meshgrid(-(r):(sz-r-1),-(r):(sz-r-1));
particlemask=((x.^2+y.^2)<=r^2);
particlemask=((x.^2+y.^2)>r^2)-particlemask;

% go through list of particle locations
% make a circle of radius feature at each location with value -1
imagesignal = image;
for i=1:size(pk,1);
    x = pk(i,1);
    y = pk(i,2);
    imagesignal((y-r):(y+r),(x-r):(x+r))=...
        image((y-r):(y+r),(x-r):(x+r)).*particlemask;
end 

imagesignal = -imagesignal;

% Finally, mask the image to exclude vignette regions
imagenoise = vignette(imagenoise,maskradius,-1);
imagesignal = vignette(imagesignal,maskradius,-1);

% Now that we have maskimage, let's calculate the noise
% Statistics of all regions with positive values > 0
index = find(imagenoise>0);
noise = image(index);
index = find(imagesignal>0);
signal = image(index);

imagesc(imagenoise./max(imagenoise(:)));
disp(['Mean signal: ', num2str(mean(signal))]);
disp(['Mean noise: ', num2str(mean(noise))]);
disp(['Standard deviation noise: ', num2str(std(noise))]);
disp(['Signal-to-noise: ', num2str((mean(signal)-mean(noise))/std(noise)),...
    '(',num2str(10*log10((mean(signal)-mean(noise))/std(noise))),' dB)' ]);

end
