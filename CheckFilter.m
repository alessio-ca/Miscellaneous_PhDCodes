function out=CheckFilter(image,filter_list)
% USAGE:    CheckFilter
% PURPOSE:  Check that a reprise with previous
%           filters is OK.
%           
% 
% INPUT: Image (assumed to be already loaded via imread)
%
% OUTPUT: Filter mask appropriate for current frame
%
% CREATED: Alessio Caciagli, University of Cambridge, 14/08/2017

out=-1;
count = size(filter_list,1);
filter_indic = count;
filter_MY = filter_list{count,1};
while out == -1

    % Background substraction (comment one of the two)
    
    %Non-uniform 
    %background = imopen(image,strel('disk',feature)); 
    
    %Uniform
    background = mean2(double(image));
    
    %Background substraction
    image_ref=(double(image) - background)./double(max(image(:)));
    norm=filter_MY./sum(sum(abs(filter_MY)));
    
    %Cross correlation between mask and background substracted image
    c = normxcorr2_general(norm,image_ref,size(norm,1)*size(norm,2)*0.2); %Consider an overlap of minimum 20% of norm (zero-pad the boundaries)
    [ypeak, xpeak] = find(c==max(c(:))); 
    peakint=max(c(:));
    
    if (ypeak-length(filter_MY)+1 <= 0) || (xpeak-length(filter_MY)+1 <= 0) || (peakint < 0.6)
        %Filter is not good. Scale the counter
        if (filter_indic<size(filter_list,1) && filter_indic>0)
            count = size(filter_list,1);
        else
            count=count-1;
        end
        
        if count==0
            %No filter is OK
            out = 0;
        else
            %Still filters to check
            filter_MY=filter_list{count,1};
        end
    else
        %Last checked filter is good. Go on.
        out = count;
    end
end
        
            
