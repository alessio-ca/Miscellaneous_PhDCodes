nFirstNice = 1;
nLastNice = nSlice;
nStartSlice = 1;   % slice 26 in 40xlens_185slices_4xaveraged is corrupted

A = uint8(zeros(nPixOut, nPixOut, nSlice));   %185

if Grayscl                  % IJ: Reading the image.
    for i=1:nSlice
        [A(:,:,i)]=imread(fname,i); 
    end
else
    for i=1:nSlice 
        [A(:,:,i)] = rgb2gray(imread(fname,i));
    end
end

targethist = imhist(A(:,:,nTarget));
imshow(A(:,:,1))
figure, imshow(A(:,:,nTarget))
figure, imshow(A(:,:,nSlice))
for i=1:nSlice     
    A(:,:,i)=histeq(A(:,:,i),targethist);
end 
figure, imshow(A(:,:,1))
figure, imshow(A(:,:,nTarget))
figure, imshow(A(:,:,nSlice))

% Correction for background intensity variation
% 4nov02ave.tif is the average intensity of 4nov02norm.tif
% as calculated using the z-project function in ImageJ
% Found that the best bet was the average over a limited range
bkg=zeros(nPixOut, nPixOut);
for i=nFirstNice:nLastNice
    bkg=bkg+double(A(:,:,i));
end
bkg=bkg/(nLastNice-nFirstNice+1);
%I think I should now fit a background function
%Use that suggested by Russ on p185
%B(x,y)=a1+a2.x+a3.y+a4.x^2+a5.y^2+a6.xy
%
[x,y]=meshgrid(1:nPixOut,1:nPixOut);

TolX = 1.e-3; TolFun = 1.e-3; MaxFunEvals = 1.0e+9; %MaxFunEvals = 100; MaxFunEvals = 1.e+8
Options = optimset('TolX',TolX,'TolFun',TolFun,'MaxFunEvals',MaxFunEvals);
astart=[180 0 0 0 0 0];                         % starting point for minimization
afinish=fminsearch('bkgfunc',astart,Options,double(bkg),x,y)
%These are parameter values obtained using the code above
%afinish=[191.5518    0.0355   -0.0603   -0.0001    0.0000   -0.0000]; 
bkgfit=afinish(1)+afinish(2)*x+afinish(3)*y+afinish(4)*(x.*x)+afinish(5)*(y.*y)+afinish(6)*(x.*y);

for i=1:nSlice %Process code
     bkgcorr=imdivide(double(A(:,:,i)),bkgfit);
     bkgcorr=uint8(bkgcorr*max(max(bkgfit)));
     % This code only required to create a tif stack
     A(:,:,i)=bkgcorr;
end 
figure, imshow(A(:,:,1))
figure, imshow(A(:,:,nTarget))
figure, imshow(A(:,:,nSlice))
save(outName,'A')