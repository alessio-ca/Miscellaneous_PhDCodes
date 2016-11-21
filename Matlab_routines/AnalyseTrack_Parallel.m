function [pos_list,histostep]=AnalyseTrack_Parallel(filehead,first,last,filter_list,feature)
% 
% USAGE:   AnalyseTrack_Parallel(filehead,first,last,filter,feature)
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
% filehead: A string common to the image file names
%           i.e. we assume that the filename is of the form 'framexxxx.tif'
%
% start:    First frame to read
% end:      Final frame to read
%
% filter:   Mask used for cross correlation
% feature:  Maximum diameter of the particles (for uneven background substraction) (optional)
%
% NOTES:
%
% OUTPUT:  list of particle positions
%          histogram of steps
%
% CREATED: Alessio Caciagli, University of Cambridge, 11/08/2016

if nargin < 5
    feature=0;
end 
pos_list=zeros(last-first+1,3);
filter=filter_list{1};
norm=filter./(sum(filter(:)));

%Enable progress bar for parallel pool
try
    parpool;
catch ME
    if ~strcmp(ME.identifier,'parallel:convenience:ConnectionOpen')
        rethrow(ME)
    end
end
warning('off','MATLAB:Java:DuplicateClass')
pctRunOnAll javaaddpath java
progressStep=5000;
ppm = ParforProgMon('Progress: ', last-first+1, progressStep, 300, 80);
disp('Tracking in progress...');

parfor frame=1:last-first+1
    
    % read in file
    image = double(imread([filehead, num2str(frame+first-1,'%06u'),'.tiff']));
    
    
    
    % Background substraction (comment one of the two)
    
    %Non-uniform
    %background = imopen(image,strel('disk',feature)); 
    %Uniform
    background = mean2(image);
    
    image_ref=mat2gray(image - background); %Background substraction
    
    % Cross correlation

    c = normxcorr2(norm,image_ref); 
    [ypeak, xpeak] = find(c==max(c(:))); 
    peakint=max(c(:));
    
    % If peakint is too faint or peak location is too close to edges (less than filter/2 length),
    % skip the point
    if all([ypeak,xpeak]-[floor(size(norm,2)/2),floor(size(norm,1)/2)]>0)==0 || all([size(c,1),size(c,2)]-[ypeak,xpeak]>[floor(size(norm,2)/2),floor(size(norm,1)/2)])==0
        ypeak=0;
        xpeak=0;
    end
    
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
        cnt = [frame+first-1,xpeak,ypeak];
        cornercheck=[sum((cnt(2:3)-[1,1]).^2),sum((cnt(2:3)-[1,size(image,1)]).^2),sum((cnt(2:3)-[size(image,2),1]).^2),sum((cnt(2:3)-[size(image,2),size(image,1)]).^2)];
        
        % If centroid is detected in corners, zero it
        if any(cornercheck<(max(size(filter,1),size(filter,2))/2).^2)
            cnt(2:3)= [0,0];
        end
        % Write the new frame to the poslist table
        pos_list(frame,:) = cnt;
        
        % If point is not valid, outputs a 0
    else
        pos_list(frame,:) = [frame+first-1,0,0];
    end

    if mod(frame,progressStep)==0
        ppm.increment();
    end
end
ppm.delete()

disp(['Number of lost frames: ', num2str(length(pos_list(arrayfun(@(x) isequal(pos_list(x,2:3),[0 0]),1:size(pos_list,1))')))])
pos_list=sortrows(pos_list,1);
histostep=pos_list(:,2:3)~=0;
histostep=histostep(:,1)&histostep(:,2);
histostep=diff(pos_list(histostep,2:3)).^2;
histostep=sqrt(histostep(:,1)+histostep(:,2));


% It's better to separate the particle tracking from the trajectory
% analysis. Once the particles are located in each frame, we can re-run the
% trajectory anslysis using different parameters.
%

disp('Tracking complete!')
end


