clear all;
img=imread('F:\AC_2018_01_11\ME_F108_SiOil_Frames_IF\10x_0.45NA_ODs_ME_SiOil_F108_10Vst_0.5mlmin.11Jan2018_19.48.31_000000.tiff');


%Detector parameters (see sections for typical values)
lpDiam = 10;
CannyParam = [0.5,5.];
dilDiam = 6;
radCorr = 1.15;
mindiam = 50;
maxdiam = 400;


%Contrast enhancement
imgJ=adapthisteq(img);

figure(1)
imshowpair(img,imgJ,'montage')

figure(2)
subplot(1,2,1)
imhist(img,64)
subplot(1,2,2)
imhist(imgJ,64)

%Low-pass filtering and high-pass filtering
%Typical values for lpDiam:
% - 10
h = fspecial('disk',lpDiam);
imgJ2 = imfilter(imgJ,h);
imgJ3 = imsubtract(imgJ2,imgJ);
imgJ4 = adapthisteq(imgJ3);

figure(3)
imshowpair(imgJ2,imgJ4,'montage')

%%
%Canny-edge detection & dilation
%Typical values for Canny:
% - 0.3,5.
% - 0.4,1.
% - 0.4,5.
% - 0.6,1.
% - 0.6,7.
%Typical values for dilDiam:
% - 1
% - 3
% - 5
% - 10
bw=edge(imgJ4,'canny',CannyParam(1),CannyParam(2));
se=strel('disk',dilDiam);
bw2=imdilate(bw,se);

figure(4)
imshowpair(bw,bw2,'montage')
bw2=not(bw2);

%Connecting components & filling holes
obj=bwconncomp(bw2);
LM = labelmatrix(obj);
stats = regionprops(LM, 'MajorAxisLength','MinorAxisLength','Area');
[~,idx]=find((mindiam<=[stats.MinorAxisLength])& ([stats.MajorAxisLength]<=maxdiam));

bw3 = ismember(LM,idx);
bw4=imfill(bw3,'holes');

figure(5)
imshowpair(bw3,bw4,'montage')

%Fit with circles
obj=bwconncomp(bw4);
props = regionprops(obj,'ConvexHull');
ND_temp = obj.NumObjects;
diami=zeros(ND_temp);
ND=0;
mask=zeros(size(bw4));
CircleFitFail=zeros(size(bw4));

dropPos = zeros(ND_temp,2);
for j=1:ND_temp
  
    if(size(props(j).ConvexHull,1)<=24) %Check minimum numbers of points for a circle fit
        continue
    end
    x=props(j).ConvexHull(:,1);
    y=props(j).ConvexHull(:,2);
    k=convhull(x,y);
    XY=[x(k),y(k)];
    [CC,RC]=fitcircle(XY);
    RC=radCorr*RC;
    if(2*RC < mindiam || 2*RC > maxdiam || ... %Check on size
            any(XY(:) < 10) || ... %Check on low border
            any(XY(:,1) > size(bw4,2)-10) || ... %Check on high border, coord1
            any(XY(:,2) > size(bw4,1)-10))      %Check on high border, coord2
        
        tempMask = drawDisk(mask,CC,RC);
        CircleFitFail=CircleFitFail | tempMask;

        continue
    end
    ND=ND+1;
    dropPos(ND,:)=CC;
    tempMask = drawDisk(mask,CC,RC);
    mask=mask | tempMask;
    diami(j)=2*RC;    
end
diami = diami(diami~=0);
%%
figure(6)
imshow(img);
green = cat(3, zeros(size(img)), ones(size(img)), zeros(size(img)));
red = cat(3, ones(size(img)), zeros(size(img)), zeros(size(img)));
hold on
h = imshow(green);
g = imshow(red);
hold off
set(h, 'AlphaData', 0.5*mask)
set(g, 'AlphaData', 0.5*CircleFitFail)

diami*(100/169.33)



    
        

            










