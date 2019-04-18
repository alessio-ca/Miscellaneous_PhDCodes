function [pos_list,filter_list]=AnalyseTrack_Movie(file,filter_list,feature,step_max,pos_list)
% 
% USAGE:   AnalyseTrack_Movie(file,filter,feature,step_max,pos_list)
% PURPOSE:  
%           For each image in a .movie file:
%               Read image file
%               Substract background
%               Cross correlate with filter
%               Find max position in correlation matrix
%               Refine max position with 2nd order polynomial fitting
%               Concatenate positions to pos_list
%
%           If concatenation is unsuccessful, prompts the user to redefine
%           the filter or to ignore the frame and proceed until tracking
%           reprises.
%           
% 
% INPUT:
% file:     A 3D array containing the fresh part of the .movie file
% filter:   Mask used for cross correlation
% feature:  Maximum diameter of the particles (for uneven background substraction) (optional)
% step_max: Maximum allowed difference between consecutive frames in pixels (optional, default=5)
% pos_list: List of particle positions (optional,reprised from previous execution on
%           older parts of the .movie file if specified, default = null)
%
% NOTES:
%
% OUTPUT:  
% pos_list: list of particle positions
% filter_list: list of used filters
%
% CREATED: Alessio Caciagli, University of Cambridge, 27/09/2016

%Variable number of inputs
if nargin < 5
    pos_list=[];
    frameoffset=0;
    firstframe = 1;
    squareX = 1;
    squareY = 1;
    edge_pad = 1;
else
    frameoffset=pos_list(end,1);
    firstframe = 0;
end

if nargin < 4
    step_max=5; %Default step max
end

if nargin < 3
    feature=0;
end 

%Definition stuff
filter=filter_list{1,1};
filter_pres=filter_list{1,2};
count=1;
loop=0;
looprep='A';
filter_indic=1;
lastframe=0;
frame=1;
disp(' ')
disp(['Tracking active with step max ',num2str(step_max)]);
disp(' ')
%Tracking loop
while(frame<=size(file,1))
    
    if mod(frame,1000)==0
        if frame ~= lastframe
            disp(['Frame number: ' num2str(frame+frameoffset)]);
        end
        lastframe=frame;
    end
    
    % Read in image
    image = squeeze(file(frame,:,:));
    
    % Background substraction (comment one of the two)
    
    %Non-uniform
    %background = imopen(image,strel('disk',feature)); 
    
    %Uniform
    background = mean2(double(image));
    
    
    % Cross correlation and corr max calculation 
    image_ref=(double(image) - background)./double(max(image(:))); %Background substraction
    if ~(firstframe && frame == 1)
        %Window size is max between 2*stepmax & 2*filter dimension
        winsize = max(step_max,length(filter));
        
        %If previous pos is too close to edge, cut the winsize
        wincenter = [round(pos_list(size(pos_list,1),2)),round(pos_list(size(pos_list,1),3))];
        squareX = (wincenter(1)-2*winsize+1:wincenter(1))+winsize;
        squareY = (wincenter(2)-2*winsize+1:wincenter(2))+winsize;
        squareX = squareX(squareX>0);
        squareX = squareX(squareX<=size(image,2));
        squareY = squareY(squareY>0);
        squareY = squareY(squareY<=size(image,1));
        
        if any(wincenter<ceil(length(filter)/2)) || any(fliplr(size(image)) - wincenter < ceil(length(filter)/2))
            edge_pad = 1;
        else
            edge_pad = 0;
        end
            
        %Cut the window from the full image
        image_ref = image_ref(squareY,squareX);
    end
    
    % Cross correlation and corr max calculation    
    norm=filter./sum(sum(abs(filter)));
    %c = normxcorr2_general(norm,image_ref,size(norm,1)*size(norm,2)*0.2); %Consider an overlap of minimum 20% of norm (zero-pad the boundaries)
    
    if edge_pad == 0
        c = normxcorr2_general(norm,image_ref,numel(norm)); %Consider complete overlap
    else
        c = normxcorr2_general(norm,image_ref,floor(numel(norm)/2)); %Consider half overlap
    end
    
    [ypeak, xpeak] = find(c==max(c(:))); 
    peakint=max(c(:));
        
    % Refine location estimates using polynomial fitting
    subc = c((ypeak-2:ypeak+2),(xpeak-2:xpeak+2));
    fit=ParabFitting(subc);
    xpeak=(squareX(1)-1) + xpeak + fit(1) - 3 - floor(size(filter,2)/2);
    ypeak=(squareY(1)-1) + ypeak + fit(2) - 3 - floor(size(filter,1)/2);
    cnt = [frame+frameoffset,xpeak,ypeak];
    
    %This section activates when tracking is attempted to be reprised.
    if loop==1
        stepdist=sqrt(sum((cnt(2:3) - cnttemp(2:3)).^2,2));
        % Check if the step is too large (due to probable error in particle
        % tracker) or if the peak intensity is too faint (non-correct
        % particle localization).
        if stepdist>step_max || peakint<0.6
            %If check is negative, loop thorugh all the available
            %filters.
            count=count-1;
            if count==0
                %If all filters were checked, skip the frame. Update
                %cnntemp, reset count and filter and upcount frame.
                cnttemp=cnt;
                count = size(filter_list,1);
                filter = filter_list{count,1};
                filter_pres = filter_list{count,2};
                frame=frame+1;
            else
                filter=filter_list{count,1};
                filter_pres = filter_list{count,2};                   
            end
            frame=frame-1;
        end
        %If check has been passed, show proposed point for reprise.
        if ~(stepdist>step_max) && ~(peakint<0.6)
            ShowFilter(mat2gray(image),filter_pres,cnt(2:3),[1 0 0])
            looprep='A';
            while (looprep ~= 'Y') && (looprep ~= 'N') && (looprep ~= 'Q')
                looprep=input('Reprise or interrupt? [Y/N/Q] ','s');
                if looprep=='Y'
                    loop=0;
                    if cnttemp ~= pos_list(end,:)
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
                        count = size(filter_list,1);
                        filter = filter_list{count,1};
                        filter_pres = filter_list{count,2};
                        frame=frame+1;
                    end
                    frame=frame-1;
                elseif looprep=='Q'
                    frame = size(file,1) + 1;
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
        stepdist=sqrt(sum((cnt(2:3) - pos_list(end,2:3)).^2,2));
        %As before
        if stepdist>step_max || peakint<0.6
            if (filter_indic<size(filter_list,1) && filter_indic>0)
                count = size(filter_list,1);
            else
                count=count-1;
            end
            %As before
            if count==0
                loop=3;
                %Decide whether to redefine filter or skip until things
                %are better
                while (loop ~= 0) && (loop ~= 1) && (loop ~= 2)
                    ShowFilter(mat2gray(image),filter_pres,cnt(2:3),[1 0 0])
                    loop=input(['Error in tracking (Frame ',num2str(frame),'). Would you like to redefine the mask (0), skip until tracking reprises (1) or interrupt (2)? ']);
                    if ~isscalar(loop)
                        loop = 3;
                    end
                    switch(loop)
                        case 0
                        case 1
                            count = size(filter_list,1);
                            filter = filter_list{count,1};
                            filter_pres = filter_list{count,2};
                            frame = frame + 1;
                            cnttemp=pos_list(end,:);
                        case 2
                            frame = size(file,1) + 1;
                        otherwise
                            disp('Incorrect input from user.')
                            loop = 3;
                    end
                end
                %This section activates when you want to redefine the
                %mask
                if loop==0
                    disp('Redefine the mask.');
                    count = size(filter_list,1)+1;
                    %Create filter on unprocessed image
                    [filter_list{count,1},filter_list{count,2}]=CreateFilter(image);
                    filter = filter_list{count,1};
                    filter_pres = filter_list{count,2};

                    %If filter is a duplicate or a preexisting one,
                    %delete it.
                    
                    if any(cellfun(@(X) isequal(diag(filter_list{end,1}),diag(X)),filter_list(1:end-1,1)))
                        filter_list(count,:)=[];
                        count=count-1;
                        filter=filter_list{count,1};
                        filter_pres=filter_list{count,1};
                    end
                end
                    
            else
                filter=filter_list{count,1};
                filter_pres=filter_list{count,1};
            end
            frame = frame - 1;
        end
    end
    
    if isempty(pos_list)
        stepdist=0;
    end
    % Concatenate the new frame to the existing data
    pos_list = [pos_list;cnt];
    filter_indic=count;
    
    % Check that new entry meets all the requirements
    if stepdist > step_max || peakint<0.6 || looprep=='N'
        pos_list(end,:)=[];
        filter_indic=0;
    end
    frame=frame+1;
end
%To have full compatibility with other data analysis routines, frames that
%were unsucessful in tracking are padded with zero
steps=diff(pos_list(:,1));
idx=find(steps~=1);
while(~isempty(idx))
    pos_list=[pos_list(1:idx(1),:);[pos_list(idx(1),1)+1,0,0];pos_list(idx(1)+1:end,:)];
    idx=idx+1;
    idx(1)=[];
end

% It's better to separate the particle tracking from the trajectory
% analysis. Once the particles are located in each frame, we can re-run the
% trajectory anslysis using different parameters.

disp('Partial tracking complete.')
end



