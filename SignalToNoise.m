%Routine to perform a pre-track SNR analysis on a sample frame. Useful to
%individuate the correct EXCESS size of the particle to have a good SNR.
%The EXCESS size is a size larger than the particle's radius but smaller
%than the spacing between particles (if doing multiparticle tracking)

%N.B! This is a "trial-and-error" process. Perform the first analysis by
%guessing an initial size for the particles. The final histogram should
%show two separated bell-shaped curves (representing signal and noise). If
%not, change the value for the particle size and reiterate the process.
%Moreover, the signal-to-noise ratio should at least be higher than 5dB.

close all
%Input: enter mother directory name (it must end with a "\" sign)
folder1='F:\2016-10-07 PS_stuck\';
%Input: enter child directory name
folder2='PS_PFSon_MicroscopeON_1000fps_LED0_LampON.07Oct2016_18.16.29';
%Input: enter filename stamp
stamp='PS_PFSon_MicroscopeON_1000fps_LED0_LampON.07Oct2016_18.16.29_';

pic=[folder1,folder2,'\',stamp,'000000.tiff']
%Show original image
if ndims(imread(pic))==2
    I = double(imread(pic));
else
    I = double(rgb2gray(imread(pic)));
end
I=double(imread(pic));
figure
%imshow(imadjust(I./max(I(:))))
imshow(mat2gray(I))

%INPUT: set the EXCESS size of the particle (in pixel)
featdiam = 6; 

%Perform SNR
featsize = floor(featdiam/2);
if mod(featsize,2) == 0
    featsize = featsize + 1;
end
figure
[signal,noise]=SNR(I,featsize); 

[signalhist,signalbin] = hist(signal ,(max(signal)-min(signal)));
[noisehist,noisebin] = hist(noise,(max(noise)-min(noise)));
sigpix = sum(signalhist); 
noisepix = sum(noisehist);
signalhist = signalhist/sigpix;
noisehist = noisehist/noisepix;

%Show histograms of signal and noise
figure
bar(signalbin,signalhist)
hold on
bar(noisebin,noisehist)
hold off