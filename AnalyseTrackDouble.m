function pos_list=AnalyseTrackDouble(filehead,first,last,filter_list,searchSize,min_feature)
% 
% USAGE:   AnalyseTrackDouble(filehead,first,last,filter)
% PURPOSE:
%           For each image in a directory:
%               Read image file
%               Substract background
%               Cross correlate with filter
%               Find max position in correlation matrix
%               Refine max position with 2nd order polynomial fitting
%               Concatenate positions to pos_list
%               (version for 2-particle tracking)
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
% OUTPUT:
%
% CREATED: Alessio Caciagli, University of Cambridge, 04/04/2016
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
    [peakint,ix]=sort(c(:),'descend');
    [ypeak,xpeak]=ind2sub(size(c),ix(1:searchSize*searchSize)); % position in the matrix
    peakint=peakint(1:searchSize*searchSize);
    cond1=ypeak>length(filter)+1 & (ypeak<length(c)-length(filter)-1);
    cond2=xpeak>length(filter)+1 & (xpeak<length(c)-length(filter)-1);
    ypeak=ypeak(cond1&cond2);
    xpeak=xpeak(cond1&cond2);
    peakint=peakint(cond1&cond2);
    peakdist=sqrt(sum(abs([ypeak,xpeak] - [ones(length(xpeak),1)*ypeak(1),ones(length(xpeak),1)*xpeak(1)]).^2,2));
    
    ypeak = [ypeak(1);ypeak(find(peakdist>searchSize,1))];
    xpeak = [xpeak(1);xpeak(find(peakdist>searchSize,1))]; 
    peakint=[peakint(1);peakint(find(peakdist>searchSize,1))];%Find maximum of cross correlation
    
    % Refine location estimates using polynomial fitting
    subc(:,:,1) = c((ypeak(1)-2:ypeak(1)+2),(xpeak(1)-2:xpeak(1)+2));
    subc(:,:,2) = c((ypeak(2)-2:ypeak(2)+2),(xpeak(2)-2:xpeak(2)+2));
    
    fit=ParabFitting(subc);
    xpeak(1)=xpeak(1)+fit(1,1)-3-floor(size(filter,2)/2);
    ypeak(1)=ypeak(1)+fit(1,2)-3-floor(size(filter,1)/2);
    xpeak(2)=xpeak(2)+fit(2,1)-3-floor(size(filter,2)/2);
    ypeak(2)=ypeak(2)+fit(2,2)-3-floor(size(filter,1)/2);
    cnt = [frame*ones(2,1),xpeak,ypeak];
    
    %This section activates when tracking is attempted to be reprised.
    if loop==1
        stepdist=sqrt(sum((cnt(:,2:3) - cnttemp(:,2:3)).^2,2));
        % Check if the step is too large (due to probable error in particle
        % tracker) or if the peak intensity is too faint (non-correct
        % particle localization).
        if any(stepdist>5)||any(peakint(1:2)<0.3)
            %Since I have two detection, order might have been shifted: try
            %by circularly shifting the cnt pair and check again.
            cnt=circshift(cnt,[1,0]);
            stepdist=sqrt(sum((cnt(:,2:3) - cnttemp(:,2:3)).^2,2));
            if any(stepdist>5) || any(peakint(1:2)<0.3)
                %If check is negative, loop thorugh all the available
                %filters.
                count=count-1;
                if count==0
                    %If all filters were checked, skip the frame. Update
                    %cnntemp, reset count and filter and upcount frame.
                    cnttemp=cnt;
                    count = length(filter_list);
                    filter = filter_list{count};
                    frame=frame+1;
                else
                    filter=filter_list{count};
                end
                frame=frame-1;
            end
        end
        %If check has been passed, show proposed point for reprise.
        if ~any(stepdist>5) && ~any(peakint(1:2)<0.3)
            ShowFilter(image_ref,filter,cnt(:,2:3),[0 1 0])
            looprep='A';
            while (looprep ~= 'Y') && (looprep ~= 'N')
                looprep=input('Reprise? [Y/N] ','s');
                if looprep=='Y'
                    loop=0;
                    if cnttemp ~= pos_list(end-1:end,:)
                        pos_list=[pos_list;cnttemp];
                    end
                elseif looprep=='N'
                    %If answer is negative, try if better result is
                    %achieved with other filters (if available).
                    %If all filters were checked, skip the frame. Update
                    %cnntemp, reset count and filter and upcount frame.
                    count=count-1;
                    if count==0
                        cnttemp=cnt;
                        count = length(filter_list);
                        filter = filter_list{count};
                        frame=frame+1;
                    end
                    frame=frame-1;
                else
                    disp('Incorrect input from user.')
                    looprep = 'A';
                end
            end
        end
    end
    
    %This section activates if tracking is active.
    if (~isempty(pos_list)) && loop==0
        % Check if the step is too large (due to probable error in particle
        % tracker) or if the peak intensity is too faint (non-correct
        % particle localization).
        stepdist=sqrt(sum((cnt(:,2:3) - pos_list(end-1:end,2:3)).^2,2));
        %As before
        if any(stepdist>5) || any(peakint(1:2)<0.3)
            cnt=circshift(cnt,[1,0]);
            stepdist=sqrt(sum((cnt(:,2:3) - pos_list(end-1:end,2:3)).^2,2));
            %As before
            if any(stepdist>5) || any(peakint(1:2)<0.3)
                count=count-1;
                %As before
                if count==0
                    loop=2;
                    %Decide whether to redefine filter or skip until things
                    %are better
                    while (loop ~= 0) && (loop ~= 1)
                        ShowFilter(image_ref,filter,cnt(:,2:3),[1 0 0])
                        loop=input(['Error in tracking (Frame ',num2str(frame),'). Would you like to redefine the mask or skip until tracking reprises? [0/1] ']);
                        switch(loop)
                            case 0
                            case 1
                                count = length(filter_list);
                                filter = filter_list{count};
                                frame = frame + 1;
                                cnttemp=pos_list(end-1:end,:);
                            otherwise
                                disp('Incorrect input from user.')
                                loop = 2;
                        end
                    end
                    %This section activates when you want to redefine the
                    %mask
                    if loop==0
                        disp('Redefine the mask.');
                        count = length(filter_list)+1;
                        %Why analyzeTrack give different visualization with imshow
                        %than filterPractice?
                        filter_list{count}=CreateFilterDouble(image_ref);
                        filter = filter_list{count};
                        %If filter is a duplicate or a preexisting one,
                        %delete it.
                        if any(cellfun(@(X) isequal(diag(filter_list{end}),diag(X)),filter_list(1:end-1)))
                            filter_list(count)=[];
                            count=count-1;
                            filter=filter_list{count};
                        end
                    end
                else
                    filter=filter_list{count};
                end
                frame = frame - 1;
            end
        end
    end
    
    
    %Only for initial frame! Record first entries in pos_list
    if isempty(pos_list)
        stepdist=0;
    end
    % Concatenate the new frame to the existing data
    pos_list = [pos_list;cnt];
    
    %If check was negative or you didn't authorize tracking reprise, delete
    %last (two) entries.
    if any(stepdist>5) || any(peakint(1:2)<0.3) || looprep=='N'
        pos_list(end-1:end,:)=[];
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



