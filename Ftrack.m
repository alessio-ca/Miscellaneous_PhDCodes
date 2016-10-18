function pos_list=Ftrack(filehead,first,last,feature)
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
%           i.e. we assume that the filename is of the form 'framexxxx.tif'
%
% start:    First frame to read
% end:      Final frame to read
%
% feature:  Expected size of the particle (diameter)     
%
% NOTES:
%
% OUTPUT:  
%
% CREATED: Eric M. Furst, University of Delaware, July 23, 2013
%  Modifications:

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
    image = double(imread([filehead, num2str(frame,'%04u'),'.tif']));
    
    % Bandpass filter
    imagebp = bpass(image,1,feature);
    
    % Find locations of the brightest pixels Need to change the 'th'
    % (threshold) argument to something that is specified or determined.
    % Current value is 8000. A rough guide is to accept 60-70% of the
    % brightest pixels 0.4*max(imagebp(:))
%    pk = pkfnd(imagebp, 8000, feature);
    pk = pkfnd(imagebp, 0.4*max(imagebp(:)), feature);
    
    % Refine location estimates using centroid
    cnt = cntrd(imagebp, pk, feature+2);
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
