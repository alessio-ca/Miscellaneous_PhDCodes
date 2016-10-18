close all, clear all;
nSlice  = 204;
nPixOut = 512;
titl    = 'Series004';
fname=[titl '.tif'];
outName = [titl '.mat'];
Grayscl = 1;

    bigstackprocJB2a

close all
clearvars -except titl

fname = [titl '.mat']
out_file = [titl '.txt']
out_test = [titl '_test.txt']
FDres    = [titl '_FDres.mat'];
boxX=50;
boxY=boxX;
rFirst=.75;
maxjj=70;
stepR=0.05;
RfConst=1.0e8;
maxNoise = 0.1;
maxnumObjects=700;
maxRegionsPerRadius = 500;
minRadiiNum=4;
ThreshMaxI          = 0.25;
Thresholding        = 0.75 * ThreshMaxI;
maxnumPixels        = 30000;

nSn_start      = .94;
nSn_i = .005;
for i=1:50
    nSn = nSn_start+i*nSn_i;
    boxZ=nSn*50.14;
    FD_nSn_opt
    load('temp_spheres')
    out(i,1) = nSn;
    out(i,2) = mean(temp_spheres(:,5));
end

save('nSn_trials', 'out', '-ascii')