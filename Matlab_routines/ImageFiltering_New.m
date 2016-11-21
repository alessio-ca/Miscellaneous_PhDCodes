%Routine to perform a pre-track ImageFiltering analysis on a sample frame. Useful to
%individuate the correct parameters for a good filtering, necessary for a
%correct particle tracking

%N.B! This routine serves to see if the image filtering procedure is
%appropriate for your serie of image. Carefully check the filtered image
%before doing the real tracking on the overall serie. You are warned!

%The refinement routine is a centroid location algorithm which provides
%sub-pixel resolution for particle tracking. It does correct the position of the brightest pixel and excludes particles which are too close to the pic edges. See the documentation of
%cntrd.m

close all
ext=input('Indicate folder of your files ','s');
%cd(ext)
ext=input('Indicate extension of your files ','s');
pic='\\np-nobelium\OE_Personal\ac2014\Desktop\Oil_Colloids_DJ\Diffusivity_2014_07_21_331fps\solution_2.21Jul2014_22.42.27\solution_2.21Jul2014_22.42.27_000382.tiff';
ext='\\np-nobelium\OE_Personal\ac2014\Desktop\Oil_Colloids_DJ\hig_exp_low_fps_50fps\Solution_Colloid_3_2014_07_18\Solution_Colloid.18Jul2014_14.01.45_001081.tiff';


%Show original image
figure
imshow(pic)
hold on
M = size(imread(pic),1);
N = size(imread(pic),2);
for k = 1:20:M
    x = [1 N];
    y = [k k];
    plot(x,y,'Color','w','LineStyle','-');
    plot(x,y,'Color','k','LineStyle',':');
end

for k = 1:20:N
    x = [k k];
    y = [1 M];
    plot(x,y,'Color','w','LineStyle','-');
    plot(x,y,'Color','k','LineStyle',':');
end

hold off

%Select appropriate conversion based on pic format
if ismatrix(imread(pic))
    I = double(imread(pic));
else
    I = double(rgb2gray(imread(pic)));
end

Ibscaled=I./max(max(I));
%Show scaled image
figure
imshow(Ibscaled)

%Apply convo filter (smears out noise)& contrast enhance the picture
threshold=0.2*max(Ibscaled(:));
im=imadjust(bpassA(1-Ibscaled,1,1,threshold));
figure
imshow(im)
figure
im_enh=imadjust(im-imopen(im,strel('square',15)));
imshow(im_enh)

%Binarize im
level = graythresh(im_enh)*3;
while level > 1
    level = level*0.9;
end
bw = im2bw(im_enh,level);
bw = bwareaopen(bw, 50);
bw = imfill(bw,'holes');


%Find region & centriud
im_prop=regionprops(bw,im_enh,'Area','PixelIdxList','WeightedCentroid');
while size(im_prop,1) > 1 && level > 0.2
    level = level*0.9;
    bw = im2bw(im_enh,level);
    bw = bwareaopen(bw, 50);
   im_prop=regionprops(bw,im_enh,'Area','PixelIdxList','WeightedCentroid');
end
bw = imfill(bw,'holes');
[circles,radii,metric]=imfindcircles(bw,[5 15]);
figure
imshow(bw)
figure
imshow(im_enh)
hold on
viscircles(circles,radii,'EdgeColor','b')
plot(im_prop.WeightedCentroid(1),im_prop.WeightedCentroid(2),'ro','MarkerSize',1,'LineWidth',2)
 % cnt=cntrd(Ib,pk,featdiam+3-mod(featdiam,2));
 % Blue circles: particle pos after image refinement 
 % plot(cnt(:,1),cnt(:,2),'bo','MarkerSize',featdiam,'LineWidth',2)
 im_prop.WeightedCentroid
hold off