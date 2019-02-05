function [pos_list,histostep]=AnalyseTrack_Parallel_Movie(file,filter_list,ppm,progressStep,feature,pos_list)
%
% USAGE:   AnalyseTrack_Parallel_Movie(file,filter,ppm,progressStep,feature,step_max,pos_list)
% PURPOSE:
%           For each image in a directory:
%               Read image file
%               Substract background
%               Cross correlate with filter
%               Find max position in correlation matrix
%               Refine max position with 2nd order polynomial fitting
%
%           N.B! Optimized for parallelization
%
%
% INPUT:
% file:         A 3D array containing the fresh part of the .movie file
%
% filter:       Mask used for cross correlation
% ppm:          Java object that provides a progress bar for the parallel loop
% progressStep: number of steps for updating the progress bar
% feature:      Maximum diameter of the particles (for uneven background substraction) (optional)
% pos_list:     List of particle positions (optional,reprised from previous execution on
%               older parts of the .movie file if specified, default = null)
%
% NOTES:
%
% OUTPUT:
% pos_list:     list of particle positions
% histostep:   histogram of steps
%
% CREATED: Alessio Caciagli, University of Cambridge, 27/09/2016
% MODIFIED: Alessio Caciagli, University of Cambridge, 23/12/2016

if nargin < 6
    pos_list=zeros(size(file,1),3);
    frameoffset=0;
else
    frameoffset=pos_list(end,1);
    pos_list=[pos_list;zeros(size(file,1),3)];
end

if nargin < 5
    step_max=5;
end

if nargin < 4
    feature=0;
end

filter=filter_list{1,1};
norm=filter./sum(sum(abs(filter)));



parfor frame=1:size(file,1)
    
    
    % read in file
    image =squeeze(file(frame,:,:));
    
    % Background substraction (comment one of the two)
    
    %Non-uniform
    %background = imopen(image,strel('disk',feature));
    %Uniform
    background = mean2(double(image));
    
    image_ref=(double(image) - background)./double(max(image(:))); %Background substraction
        
    % Cross correlation
    c = normxcorr2_general(norm,image_ref,size(norm,1)*size(norm,2)*0.2); %Consider an overlap of minimum 10% of norm (zero-pad the boundaries)
    [ypeak, xpeak] = find(c==max(c(:))); %Maximum of cross correlation
    peakint=max(c(:));
    
    if peakint < 0.6
        ypeak=0;
        xpeak=0;
    end
    
    
    % If point is valid, refine location estimates using polynomial fitting
    if xpeak*ypeak~=0
        subc = c((ypeak-2:ypeak+2),(xpeak-2:xpeak+2));
        fit=ParabFitting(subc);
        xpeak=xpeak+fit(1)-3-floor(size(filter,2)/2);
        ypeak=ypeak+fit(2)-3-floor(size(filter,1)/2);
        cnt = [frame+frameoffset,xpeak,ypeak];

        % Write the new frame to the poslist table
        pos_list(frame+frameoffset,:) = cnt;
        
        % If point is not valid, outputs a 0
    else
        pos_list(frame+frameoffset,:) = [frame+frameoffset,0,0];
    end
    if mod(frame,progressStep)==0
        ppm.increment();
    end
end
disp(['Number of lost frames: ', num2str(length(pos_list(arrayfun(@(x) isequal(pos_list(x,2:3),[0 0]),1:size(pos_list,1))')))])
pos_list=sortrows(pos_list,1);
histostep=pos_list(:,2:3)~=0;
histostep=histostep(:,1)&histostep(:,2);
histostep=diff(pos_list(histostep,2:3)).^2;
histostep=sqrt(histostep(:,1)+histostep(:,2));
% It's better to separate the particle tracking from the trajectory
% analysis. Once the particles are located in each frame, we can re-run the
% trajectory anslysis using different parameters.
disp('Partial tracking complete.')
end



