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
r_big = 15;
r_small = 5;
for frame=first:last
    
    if mod(frame,10)==0
        disp(['Frame number: ' num2str(frame)]);
    end
    
    % read in file
    pic=[filehead, num2str(frame,'%05u'),'.',ext];
  
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
    level = graythresh(im_enha)*2.5; %auto threshold
    while level >= 1
        level = level*0.9;
    end
    bw = im2bw(im_enha,level); %binarize with higher threshold (Might need to adjust the multiplication factor)
    bw = bwareaopen(bw, 50); %remove features smaller than 50 (noise or artifact) % DJ changed 50 to 100
    bw = imfill(bw,'holes'); %fill holes in shapes (i.e., donuts)

    %Find region & centroid
    im_reg=regionprops(bw,im_enha,'Area','PixelIdxList','WeightedCentroid');
    [~,idx]=max([im_reg.Area]);
    im_reg_max=im_reg(idx);    
    
    while isempty(im_reg)
        level = level*0.9; %auto threshold
        bw = im2bw(im_enha,level); %binarize with higher threshold (Might need to adjust the multiplication factor)
        bw = bwareaopen(bw, 50); %remove features smaller than 50 (noise or artifact) 
        bw = imfill(bw,'holes'); %fill holes in shapes (i.e., donuts)

        %Find region & centroid
        im_reg=regionprops(bw,im_enha,'Area','PixelIdxList','WeightedCentroid');
        [~,idx]=max([im_reg.Area]);
        im_reg_max=im_reg(idx);
    end
    cnt=[im_reg_max.WeightedCentroid(1),im_reg_max.WeightedCentroid(2)];
    if im_reg_max.Area > 250
       [circles,~,metric]=imfindcircles(bw,[r_small r_big]);
       cnt = circles;
       if size(circles,1) > 1
            [~,idx]=max(metric);
            cnt = circles(idx,:);
       end
    end
    if isempty(pos_list)==false
        dist=norm(cnt-pos_list(end,2:3));
        index=1;
        while dist > 5 && index <= size(im_reg,1)
            cnt = [im_reg(index).WeightedCentroid(1),im_reg(index).WeightedCentroid(2)];
            dist = norm(cnt-pos_list(end,2:3));
            index = index + 1;
        end
        if dist > 5
            level = level*1.1;
            if level > 1
                level = 0.95;
            end
            bw = im2bw(im_enha,level); %binarize with higher threshold (Might need to adjust the multiplication factor)
            bw = bwareaopen(bw, 50); %remove features smaller than 50 (noise or artifact) 
            bw = imfill(bw,'holes'); %fill holes in shapes (i.e., donuts)

            %Find region & centroid
            im_reg=regionprops(bw,im_enha,'Area','PixelIdxList','WeightedCentroid');
            while isempty(im_reg)
                level = level*0.9;
                bw = im2bw(im_enha,level); %binarize with higher threshold (Might need to adjust the multiplication factor)
                bw = bwareaopen(bw, 50); %remove features smaller than 50 (noise or artifact) 
                bw = imfill(bw,'holes'); %fill holes in shapes (i.e., donuts)

                %Find region & centroid
                im_reg=regionprops(bw,im_enha,'Area','PixelIdxList','WeightedCentroid');
               
                [~,idx]=max([im_reg.Area]);
                im_reg_max=im_reg(idx);
            end
            cnt=[im_reg_max.WeightedCentroid(1),im_reg_max.WeightedCentroid(2)];
            if im_reg_max.Area > 250
                [circles,~,metric]=imfindcircles(bw,[r_small r_big]);
                cnt = circles;
                if size(circles,1) > 1
                    [~,idx]=max(metric);
                    cnt = circles(idx,:);
                end
            end
            dist=norm(cnt-pos_list(end,2:3));
            index=1;
            while dist > 5 && index <= size(im_reg,1)
                cnt = [im_reg(index).WeightedCentroid(1),im_reg(index).WeightedCentroid(2)];
                dist = norm(cnt-pos_list(end,2:3));
                index = index + 1;
            end
            if dist > 5
                level = graythresh(im_impr)*3; %auto threshold
                while level >= 1
                    level = level*0.9;
                end
                bw = im2bw(im_impr,level); %binarize with higher threshold (Might need to adjust the multiplication factor)
                bw = bwareaopen(bw, 100); %remove features smaller than 50 (noise or artifact) % DJ changed 50 to 100
                bw = imfill(bw,'holes'); %fill holes in shapes (i.e., donuts)

                %Find region & centroid
                im_reg=regionprops(bw,im_impr,'Area','PixelIdxList','WeightedCentroid');
                [~,idx]=max([im_reg.Area]);
                im_reg_max=im_reg(idx);    
    
                while isempty(im_reg)
                    level = level*0.9; %auto threshold
                    bw = im2bw(im_impr,level); %binarize with higher threshold (Might need to adjust the multiplication factor)
                    bw = bwareaopen(bw, 50); %remove features smaller than 50 (noise or artifact)
                    bw = imfill(bw,'holes'); %fill holes in shapes (i.e., donuts)

                    %Find region & centroid
                    im_reg=regionprops(bw,im_impr,'Area','PixelIdxList','WeightedCentroid');
                    [~,idx]=max([im_reg.Area]);
                    im_reg_max=im_reg(idx);
                end
                cnt=[im_reg_max.WeightedCentroid(1),im_reg_max.WeightedCentroid(2)];
                if im_reg_max.Area > 250
                    [circles,~,metric]=imfindcircles(bw,[r_small r_big]);
                    cnt = circles;
                    if size(circles,1) > 1
                        [~,idx]=max(metric);
                        cnt = circles(idx,:);
                    end
                end
                dist=norm(cnt-pos_list(end,2:3));
                if dist > 5
                    cnt
                    pos_list(end,2:3)
                    imshow(im_enha)
                    hold on
                    %viscircles(circles,radii,'EdgeColor','b')
                    plot(cnt(1),cnt(2),'ro','MarkerSize',1,'LineWidth',2)
                    plot(pos_list(end,2),pos_list(end,3),'go','MarkerSize',1,'LineWidth',2)
                    hold off
                    cd('\\np-nobelium\OE_Personal\ac2014\Documents\MATLAB');
                    error('Problem with centroid at frame %d',frame)
                end
            end
        end
    end
    pos_list = [pos_list;[frame,cnt]];
end
cd('\\np-nobelium\OE_Personal\ac2014\Documents\MATLAB');
