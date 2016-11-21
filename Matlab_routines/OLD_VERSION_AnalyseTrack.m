function pos_list=AnalyseTrack(filehead,first,last,filter_list,feature)
% 
% USAGE:   AnalyseTrack(filehead,first,last,filter)
% PURPOSE:  
%           For each image in a directory:
%               Read image file
%               Substract background
%               Cross correlate with filter
%               Find max position in correlation matrix
%               Refine max position with 2nd order polynomial fitting
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
% filter:   Mask used for cross correlation
% feature:  Maximum diameter of the particles
%
% NOTES:
%
% OUTPUT:  list of particle positions
%
% CREATED: Alessio Caciagli, University of Cambridge, 11/12/2015

pos_list=[];
filter=filter_list{1};
count=1;
frame = first;
loop=0;
looprep='A';

while(frame<=last)
    
    if mod(frame,10)==0
        disp(['Frame number: ' num2str(frame)]);
    end
    
    % read in file
    image = double(imread([filehead, num2str(frame,'%06u'),'.tiff']));
    
    % Background substraction (comment one of the two)
    
    %Non-uniform 
    %background = imopen(image,strel('disk',feature)); 
    
    %Uniform
    background = mean2(image);
    image_ref=mat2gray(image - background); %Background substraction
    
    % Cross correlation
    norm=filter./(sum(filter(:)));
    c = normxcorr2(norm,image_ref); %Cross correlation between mask and background substracted image
    [ypeak, xpeak] = find(c==max(c(:))); %Maximum of cross correlation
    peakint=max(c(:));
    
    % Refine location estimates using polynomial fitting
    subc = c((ypeak-2:ypeak+2),(xpeak-2:xpeak+2));
    fit=ParabFitting(subc);
    xpeak=xpeak+fit(1)-3-floor(size(filter,2)/2);
    ypeak=ypeak+fit(2)-3-floor(size(filter,1)/2);
    cnt = [frame,xpeak,ypeak];
    
    
    if (~isempty(pos_list))
        % Check if the step is too large (due to probable error in particle tracker) 
        stepdist=sqrt(sum((cnt(2:3)-pos_list(end,2:3)).^2));
        if stepdist > 5
            count=count-1;
            if count==0
                ShowFilter(image_ref,filter,cnt(2:3),[1 0 0])
                disp(['Error in tracking (Frame ',num2str(frame),'). Redefine the mask.']);
                count = length(filter_list)+1;
                %Why analyzeTrack give different visualization with imshow
                %than filterPractice?
                filter_list{count}=CreateFilter(image_ref);
                filter = filter_list{count};
                if any(cellfun(@(X) isequal(diag(filter_list{end}),diag(X)),filter_list(1:end-1)))
                    filter_list(count)=[];
                    count=count-1;
                    filter=filter_list{count};
                end
            else
                filter=filter_list{count};
            end
            frame = frame - 1;
        end
    end
    
    if isempty(pos_list)
        stepdist=0;
    end
    % Concatenate the new frame to the existing data
    pos_list = [pos_list;cnt];
    
    if stepdist > 5
        pos_list(end,:)=[];
    end
    frame=frame+1;
    
end

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
%   trajectories = track(pos_list,10);
disp('Tracking complete!')
end



