%Routine to perform a pre-track ImageFiltering analysis on a sample frame.
%Useful to individuate the correct parameters for a good filtering,
%necessary for the correct particle tracking in Fluo.

%This routine serves to see if the image filtering procedure is
%appropriate for your serie of image. It shows the location of the detected
%particles in red circles (blue circles are unrefined positions which do
%not appear in the final result).

%Carefully check the filtered image before doing the real tracking on the
%overall serie.

%The refinement routine is a centroid location algorithm which provides
%sub-pixel resolution for particle tracking. See the docuentation of
%cntrd.m

close all
%Input: enter mother directory name (it must end with a "\" sign)
folder1='F:\2016-10-07 PS_stuck\';
%Input: enter child directory name
folder2='PS_PFSon_MicroscopeON_1000fps_LED0_LampON.07Oct2016_18.16.29';
%Input: enter filename stamp
stamp='PS_PFSon_MicroscopeON_1000fps_LED0_LampON.07Oct2016_18.16.29_';

pic=[folder1,folder2,'\',stamp,'000000.tiff']

%Show original image
figure
imshow(pic)

featdiam=5; %should be the same size got in the SignalToNoise pre-tracking, if you tried it.

if ndims(imread(pic))==2
    I = double(imread(pic));
else
    I = double(rgb2gray(imread(pic)));
end
featsize = floor(featdiam/2);
if mod(featsize,2) == 0
    featsize = featsize + 1;
end
Ib = bpass(I,1,featsize); %Apply a high and low-pass filters. See bpass documentation for info about the parameters
Ibscaled=Ib./max(max(Ib));
%Show contrasted image
figure
imshow(mat2gray(Ibscaled))

threshold=0.5*max(Ib(:)); %Set threshold for identifying the brightest pixel (60% of max intensity of the contrasted image)
pk = pkfnd(Ib,threshold,floor(featdiam/2));
%Show contrasted image with superimposed positions of the particles
figure
imshow(Ibscaled)
hold on
%Blue circles: particle pos before image refinement
plot(pk(:,1),pk(:,2),'bo','MarkerSize',featdiam,'LineWidth',2)
cnt=cntrd(Ib,pk,featdiam+3-mod(featdiam,2));
%Red circles: particle pos after image refinement
plot(cnt(:,1),cnt(:,2),'ro','MarkerSize',featdiam,'LineWidth',2)
hold off