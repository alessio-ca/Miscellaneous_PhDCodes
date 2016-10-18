function pos_list=FtrackA(filehead,ext,first,last,feature,threshold,singletrack)
% 
% USEAGE:   pos_list=Ftrack(filehead,first,last,feature);
% PURPOSE:  
%           For each image in a directory
%              Read image file
%               Filter the image (bpass.m)
%               Find brightest pixels (pkfnd.m)
%               Correct position with centroid (cntrd.m)
%               Concatenate positions to pos_list
%           
% 
% INPUT:
% filehead: A string common to the image file names
%           i.e. we assume that the filename is of the form 'framexxxx'
%
% ext:      The extension of the files in the serie (e.g., png)
%
% first:    First frame to read
% last:      Final frame to read
%
% feature:  Expected size of the particle (diameter) 
% threshold: Threshold for locating the brightest pixels (from 0 to 1)
% singletrack: specify if you are going to do single or multiple particle
% tracking
%
% NOTES:
%
% OUTPUT:  
%
% CREATED: Eric M. Furst, University of Delaware, July 23, 2013
% MODIFIED: Alessio Caciagli, University of Cambridge, October 15, 2015

if mod(feature,2) == 0
    warning('feature size must be an odd value');
    out=[];
    return;
end

pos_list=[];
for frame=first:last
    
    if mod(frame,10)==0
        disp(['Frame number: ' num2str(frame)]);
    end
    
    % read in file
    pic=[filehead, num2str(frame,'%05u'),'.',ext];
  
    
    if ndims(imread(pic))==2
        image = double(imread(pic));
    else
        image = double(rgb2gray(imread(pic)));
    end
    
    % Bandpass filter
    imagebp = bpass(image,1,feature);
   
    %DEBUGGING STUFF
%     Ibscaled=imagebp./max(max(imagebp));
%     figure
%     imshow(Ibscaled)
    
    % Find locations of the brightest pixels Need to change the 'th'
    % (threshold) argument to something that is specified or determined.
    % Current value is 8000. A rough guide is to accept 60-70% of the
    % brightest pixels
    % pk = pkfnd(imagebp, 8000, feature);
    pk = pkfnd(imagebp, threshold*max(imagebp(:)), feature);
     if str2double(singletrack) == 0
             im_prop=regionprops(im2bw(imagebp./max(max(imagebp)),0.2),'Area','PixelIdxList');
             %im_prop([im_prop.Area] < 10)=[];
             [area_obj,idx]=min([im_prop.Area]);
             imagebp2=ones(size(imagebp));
             imagebp2(im_prop(idx).PixelIdxList)=3;
             imagebp3=imagebp.*imagebp2;
             pk = pkfnd(imagebp3, threshold*max(imagebp3(:)), feature);
             n=0;
             while size(pk,1) > 1 && n < 4
                  n=n+1;
                  im_prop=regionprops(im2bw(imagebp./max(max(imagebp)),0.2+0.1*n),'Area','PixelIdxList');
%                   im_prop([im_prop.Area] < 10)=[];
                  [area_obj,idx]=min([im_prop.Area]);
                  imagebp2=ones(size(imagebp));
                  imagebp2(im_prop(idx).PixelIdxList)=3;
                  imagebp3=imagebp.*imagebp2;
                  pk = pkfnd(imagebp3, threshold*max(imagebp3(:)), feature);
             end
%              if isempty(pos_list)==false &&(size(pk,1) ~= 1 ||  norm(pk'-round(pos_list(1:2,end)))>5)
%                     image = 255-image;
%                     imagebp = bpassA(image,1,feature); %Apply a high and low-pass filters. See bpass documentation for info about the parameters
%                     im_prop=regionprops(im2bw(imagebp./max(imagebp(:)),0.9),imagebp./max(imagebp(:)),'Area','PixelIdxList','WeightedCentroid');
%                     pk=double(int8(im_prop.WeightedCentroid));
%              end
             imagebp=imagebp3;
             clear imagebp3;
    end
    % Refine location estimates using centroid
    cnt = cntrdA(imagebp, pk, feature+2);
    % cntrd can also accept "interactive" mode
    % cnt = cntrd(imagebp, pk, feature+2, 1);
    
    % Add frame number to tracking data
    cnt(:,5) = frame;
    
    % You can also add pk to cnt to check if
    % corrections are greater than one pixel.
    %     cnt = [pk(:,2), cnt];
    %     cnt = [pk(:,1), cnt];
    
    % Concatenate the new frame to the existing data
    pos_list = [pos_list, cnt'];
    
    % Code that can be used to make images to monitor tracking
    % Uncomment to use
%     colormap(gray);
%     imagesc(image,[8000 48000]);
%     axis image
%     hold on;
%     plot(pk(:,1),pk(:,2),'go','MarkerSize',10,'LineWidth',1);
%     plot(cnt(:,1),cnt(:,2),'y.'); 
%     hold off    
%     pause(0.05); 
end

% Format the position list so that we have four columns: x, y, m_0, m_2, frame
pos_list = pos_list';

% It's better to separate the particle tracking from the trajectory
% analysis. Once the particles are located in each frame, we can re-run the
% trajectory anslysis using different parameters.
%
% After generating the unsorted position lists, run track.m to generate the
% particle trajectories. The routine takes the following input:
%      For the input data structure (positionlist):
%            (x)          (y)          (t)
%      pos = 3.60000      5.00000      0.00000
%            15.1000      22.6000      0.00000
%            4.10000      5.50000      1.00000 
%            15.9000      20.7000      2.00000
%            6.20000      4.30000      2.00000
% 
% Use the command:
%   trajectories = track([pos_list(:,1:2),pos_list(:,5)],10);
