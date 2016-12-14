function [pos_list,filter_list]=AnalyseTrack(filehead,first,last,filter_list,feature,step_max,opfilter)
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
%           If concatenation is unsuccessful, prompts the user to redefine
%           the filter or to ignore the frame and proceed until tracking
%           reprises.
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
% feature:  Maximum diameter of the particles (for uneven background substraction) (optional, default = 0)
% step_max: Maximum allowed difference between consecutive frames (optional, default = 5)
% opfilter: If active, loads the last filter in filter_list (optional,
%           default = 0)
%
% NOTES:
%
% OUTPUT:   list of particle positions
%           list of filters
%
% CREATED: Alessio Caciagli, University of Cambridge

%Variable number of inputs
if nargin < 7
    opfilter=0;
end
if nargin < 6
    step_max=5;
end
if nargin < 5
    feature=0;
end 
pos_list=[];

if opfilter == 1
    filter=filter_list{end};
else
    filter=filter_list{1};
end

%Definition stuff
count=1;
frame = first;
loop=0;
looprep='A';
filter_indic=1;
lastframe=0;

disp(' ')
disp(['Tracking active with step max ',num2str(step_max)]);
disp(' ')

%Tracking loop
while(frame<=last)
    
    if mod(frame,1000)==0
        if frame ~= lastframe
            disp(['Frame number: ' num2str(frame)]);
        end
        lastframe=frame;
    end
    
    % Read in image
    image = double(imread([filehead, num2str(frame,'%06u'),'.tiff']));
    
    
    
    % Background substraction (comment one of the two)
    
    %Non-uniform
    %background = imopen(image,strel('disk',feature)); 
    
    %Uniform
    background = mean2(image);
    
    image_ref=mat2gray(image - background); %Background substraction
    
    % Cross correlation and corr max calculation
    norm=filter./(sum(filter(:)));
    c = normxcorr2(norm,image_ref); 
    [ypeak, xpeak] = find(c==max(c(:))); 
    peakint=max(c(:));
    
    % If peak location is too close to edges (less than filter/2 length), assign peak pos randomly
    if all([ypeak,xpeak]-[floor(size(norm,2)/2),floor(size(norm,1)/2)]>0)==0 || all([size(c,1),size(c,2)]-[ypeak,xpeak]>[floor(size(norm,2)/2),floor(size(norm,1)/2)])==0
        fake=floor(floor(max(size(norm,1)/2,size(norm,2)/2)) + floor(min(size(c,1),size(c,2))-max(size(norm,1)/2,size(norm,2)/2))*rand(1,2));
        ypeak=fake(1);
        xpeak=fake(2);
    end
        
    % Refine location estimates using polynomial fitting
    subc = c((ypeak-2:ypeak+2),(xpeak-2:xpeak+2));
    fit=ParabFitting(subc);
    xpeak=xpeak+fit(1)-3-floor(size(filter,2)/2);
    ypeak=ypeak+fit(2)-3-floor(size(filter,1)/2);
    cnt = [frame,xpeak,ypeak];
    cornercheck=[sum((cnt(2:3)-[1,1]).^2),sum((cnt(2:3)-[1,size(image,1)]).^2),sum((cnt(2:3)-[size(image,2),1]).^2),sum((cnt(2:3)-[size(image,2),size(image,1)]).^2)];
    
    % If centroid is detected in corners, assign it randomly
    if any(cornercheck<(max(size(filter,1),size(filter,2))/2).^2)
        cnt(2:3)= min(size(image,1),size(image,2))*rand(1,2);
    end
    
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
                count = length(filter_list);
                filter = filter_list{count};
                frame=frame+1;
            else
                filter=filter_list{count};
            end
            frame=frame-1;
        end
        %If check has been passed, show proposed point for reprise.
        if ~(stepdist>step_max) && ~(peakint<0.6)
            filter_pres=filter;
            if (filter(1,1)==1)
                filter_pres=1-filter;
            end
            ShowFilter(image_ref,filter_pres,cnt(2:3),[0 1 0])
            looprep='A';
            while (looprep ~= 'Y') && (looprep ~= 'N')
                looprep=input('Reprise? [Y/N] ','s');
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
        stepdist=sqrt(sum((cnt(2:3) - pos_list(end,2:3)).^2,2));
        %As before
        if stepdist>step_max || peakint<0.6
            if (filter_indic<length(filter_list) && filter_indic>0)
                count = length(filter_list);
            else
                count=count-1;
            end
            %As before
            if count==0
                loop=2;
                %Decide whether to redefine filter or skip until things
                %are better
                while (loop ~= 0) && (loop ~= 1)
                    filter_pres=filter;
                    if (filter(1,1)==1)
                        filter_pres=1-filter;
                    end
                    ShowFilter(image_ref,filter_pres,cnt(2:3),[1 0 0])
                    loop=input(['Error in tracking (Frame ',num2str(frame),'). Would you like to redefine the mask or skip until tracking reprises? [0/1] ']);
                    if ~isscalar(loop)
                        loop = 2;
                    end
                    switch(loop)
                        case 0
                        case 1
                            count = length(filter_list);
                            filter = filter_list{count};
                            frame = frame + 1;
                            cnttemp=pos_list(end,:);
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
                    
                    filter_list{count}=CreateFilter(image_ref);
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
%were unsucessful in tracking are padded with zero.
steps=diff(pos_list(:,1));
idx=find(steps~=1);
while(~isempty(idx))
    pos_list=[pos_list(1:idx(1),:);[pos_list(idx(1),1)+1,0,0];pos_list(idx(1)+1:end,:)];
    idx=idx+1;
    idx(1)=[];
end


% It's better to separate the particle tracking from the trajectory
% analysis. Once the particles are located in each frame, we can re-run the
% trajectory analysis using different parameters.

disp('Tracking complete!')
end



