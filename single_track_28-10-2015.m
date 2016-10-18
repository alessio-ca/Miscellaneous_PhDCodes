function pos_list=single_track(filehead,ext,ext1,first,last)

% USEAGE:   pos_list=single_track(filehead,ext,first,last,feature);
% PURPOSE:  
%           For each image in a directory
%              Read image file
%               Filter the image (bpass.m)
%               Convolve and enhance contrast
%               Binarize and finds weighted centroid
%               Concatenate positions to pos_list
%           
% 
% INPUT:
% filehead: A string common to the image file names
%           i.e. we assume that the filename is of the form 'framexxxx'
%           n.b! It also assumes a name with 6 digits 
%
% ext:      The extension of the files in the serie (e.g., png)
%
% first:    First frame to read
% last:      Final frame to read
% 
%
% NOTES:
%
% OUTPUT:  
%
% CREATED: Alessio Caciagli, University of Cambridge, October 26, 2015
cd(ext1)
pos_list=[];
for frame=first:last
    
    if mod(frame,10)==0
        disp(['Frame number: ' num2str(frame)]);
    end
    
    % read in file
    pic=[filehead, num2str(frame,'%06u'),'.',ext];
  
    % Convert accordingly to extension
    if ismatrix(imread(pic))
        image = double(imread(pic));
    else
        image = double(rgb2gray(imread(pic)));
    end
    
    image_sc=image./max(max(image));
    
    %Apply convo filter (smears out noise)& contrast enhance the picture
    threshold=0.2*max(image_sc(:));
    im_impr=imadjust(bpassA(1-image_sc,1,1,threshold));
    im_enha=imadjust(im_impr-imopen(im_impr,strel('square',15))); %Might need to adjust the size of the opening feature


    %Binarize im
    level = graythresh(im_enha); %auto threshold
    bw = im2bw(im_enha,level*2.5); %binarize with higher threshold (Might need to adjust the multiplication factor)
    bw = bwareaopen(bw, 100); %remove features smaller than 50 (noise or artifact) % DJ changed 50 to 100
    bw = imfill(bw,'holes'); %fill holes in shapes (i.e., donuts)

    %Find region & centroid
    im_reg=regionprops(bw,im_enha,'Area','PixelIdxList','WeightedCentroid');
    if size(im_reg,1)>1 %If it finds another region bigger than 50, considers the biggest 
        [~,idx]=max([im_reg.Area]);
        im_reg=im_reg(idx);
    end
    
    cnt=[im_reg.WeightedCentroid(1),im_reg.WeightedCentroid(2)];
    if isempty(pos_list)==false
        dist=norm(cnt-pos_list(end,2:3));
        if dist > 5
            
            cd('\\np-nobelium\OE_Personal\ac2014\Documents\MATLAB');
            error('Problem with centroid at frame %d',frame)
        end
    end
    pos_list = [pos_list;[frame,cnt]];
end
cd('\\np-nobelium\OE_Personal\ac2014\Documents\MATLAB');
