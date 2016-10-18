%Before running make sure the number of pixels in x and y is 512, otherwise
%probably some major changes would have to be made. 
%
%Aspect ratio is length in microns in z direction divided by the number of 
%slices in z (number of pixels) and this all divided by the length in microns 
%in x (or y, they are usualy the same) divided by the number of pixels
%after downsizing, e.i. 256. So it is (z [micr/pix])/(x [micr/pix]).
% Also keep your eyes on handles.zsp in gui, it's the aspect ratio and
% needs to be changed for different series. Do not forget Snell
% correction!!!

close all, clear all;
nSlice  = 261;                     % Pixels in z. Should be even because of midY=nSlice/2 and nFrames/2 (p 47, 123, 127). I changed to round(..)
nPixOut = 512;                  % Pixels in xy. Should be even since there is nPix/2. I changed to round(...)    
nTarget = 133;                           % Histogram equil. In imageJ go to Image->Stack->Plot Z-axis profile and chose the highest intensity slice.
titl    = '20um-creamed-02_gui_synthNoise_constr';  % Name of input tif file
fname=[titl '.tif'];      % Original image stack presented in a single tif. I converted it to greyscale in imageJ.
outName = [titl '.mat'];       % Cleaned and equalized output.
Grayscl = 1;    % If you set 1 then it assumes you input greyscale, otherwise it will convert RGB into greyscale.

    bigstackprocJB
close all
clearvars -except titl
 disp('Hit any key to proceed to FD')
 pause
fname = [titl '.mat']      % Output of bigstackprocJB.
out_file = [titl '.txt']   % Output of FD.
out_test = [titl '_test.txt']  % Outputs all drops with garbage, but can check in gui if they at least cover mosto of the spheres.
FDres    = [titl '_FDres.mat'];
nSn      = 1.4150/1.4475 ; %1.405/1.45; Snells law correction
boxX=144.72; % Box in xy
boxY=boxX;
boxZ=nSn*261.84;            % Each step is boxZ/(nslices - 1). I here divided since the original nSlices was 213.
rFirst=10.0;                  % Shouldn't be too small, cause try to make it at least 10 pixels. Smallest sphere radius.
maxjj=20;                   % Number of steps, i.e. maxjj * stepR = the biggest sphere that will be analyzed. Cause usually r is less than 6 mum.
stepR=0.1;                 % Following Jasnas Thesis (p 108), increment size is limited by the accuracy of this image analisys method, 1/3 of a pixel. 
RfConst=1.0e8;           % Wiener signal to noise level. Assuming that noise is white and indipendent of frequency. Better image -> higher number (1e7 to 1e9).
maxNoise = 0.1;          % Noise needed for zeros in fip!!! IJ: changed to 200, was 100.
maxnumObjects=700;           % Initial filtration by threshold so that you have no more than this number of objects per radius.
maxRegionsPerRadius = 500;   % How many spheres per radius we store after the second step of filtration.
minRadiiNum=4;               % How many times the sphere should be found in order to be treated as real one. So far worked 4 as best.
ThreshMaxI          = 0.25;  % 0.5  %0.6         % filtering: only regions with max intensity > ThreshMaxI are considered. For bad images decrease this number (0.1 to 0.5). 
Thresholding        = 0.75 * ThreshMaxI; %0.8    % CM of the region is computed using voxels with intensity > Thresholding
maxnumPixels        = 30000;                    % Get rid of it. Don't change
    FD  % Enter=accept current sphere; 1+enter=reject current sphere; 2+enter=reject all that follow; 4+enter=enter number of spheres to accept in a row. 