function [filter,r_max,searchSize]=CreateFilterDouble(image)
% USAGE:        CreateFilterDouble
% PURPOSE:      Overlays a mask over an image for empirically checking whether
%               the mask is appropriate (version for 2-particle tracking)
%           
% 
% INPUT:        Image (assumed to be already loaded via imread)
%
% OUTPUT:
% filter:       appropriate mask for usage in particle tracking
% r_max:        max dimension of the mask
% searchSize:   length of the square in which corr maximum is most probably
%               located
%
% CREATED: Alessio Caciagli, University of Cambridge, 04/04/2016

out=0;
%image=imread(image);
close all
figure
imshow(imadjust(image))
while out == 0
    param=[-1,-1,-1];
    %Set security commands for wrong user input
    while any(param(1:3) < 0) || any(diff(param(1:3))>0) || any(mod(param(1:3),2) == 1)
        param=input('Insert filter parameters [d_max d_mid d_min] ');
        if (length(param) <= 3)
            filter_MY=create_filter_2(param(1),param(2),param(3)); %Create filter with specified parameters
        else
            filter_MY=create_filter_2(param(1),param(2),param(3),param(4),param(5),param(6)); %Create filter with specified parameters
        end
    end
    
    %Definition of the size of the square where max relative to one particle might be
    maxArea=param(3);
    if maxArea==0
        maxArea=param(2);
        if maxArea==0
            maxArea=param(1);
        end
    end
    
    X_Lim=get(gca,'XLim'); %Commands to keep the zoom in imshow during filter creation
    Y_Lim=get(gca,'YLim');
    NORMA=filter_MY./(sum(filter_MY(:)));
    c = normxcorr2(NORMA,image); %Cross correlation between mask and background substracted image
    [~,ix]=sort(c(:),'descend');
    [ypeak,xpeak]=ind2sub(size(c),ix(1:maxArea*maxArea)); % position in the matrix
    peakint=peakint(1:maxArea*maxArea);
    cond1=ypeak>length(filter_MY)+1 & (ypeak<length(c)-length(filter_MY)-1);
    cond2=xpeak>length(filter_MY)+1 & (xpeak<length(c)-length(filter_MY)-1);
    ypeak=ypeak(cond1&cond2);
    xpeak=xpeak(cond1&cond2);
    peakint=peakint(cond1&cond2);
    peakdist=sqrt(sum(abs([ypeak,xpeak] - [ones(length(xpeak),1)*ypeak(1),ones(length(xpeak),1)*xpeak(1)]).^2,2));
    
    ypeak = [ypeak(1);ypeak(find(peakdist>maxArea,1))];
    xpeak = [xpeak(1);xpeak(find(peakdist>maxArea,1))]; %Find maximum of cross correlation
    peakint=[peakint(1);peakint(find(peakdist>searchSize,1))];%Find maximum of cross correlation

    
    if (any(ypeak-length(filter_MY)+1 <= 0) || any(length(c)-ypeak-length(filter_MY)+1 <= 0) || any(xpeak-length(filter_MY)+1 <= 0) || any(length(c)-xpeak-length(filter_MY)+1 <= 0))
        disp('Error in correlation matrix maximum. Try redefining the mask.')
        continue
    end
    
    xpeak(1)=xpeak(1)-floor(size(filter_MY,2)/2); %Manipulation for ShowFilter (which automatically calculates the offset)
    ypeak(1)=ypeak(1)-floor(size(filter_MY,1)/2);
    xpeak(2)=xpeak(2)-floor(size(filter_MY,2)/2);
    ypeak(2)=ypeak(2)-floor(size(filter_MY,1)/2);
    
    %Refining of filter appearence in imshow as a green mask
    ShowFilter(image,filter_MY,[xpeak(1),ypeak(1);xpeak(2),ypeak(2)],[0 1 0])
    if any(peakint(1:2)<0.3)
        disp('Error in peak maximum. Position might be right, but correlation with image is poor. Try slightly redefining the mask.')
        continue
    end
    
    %Input to accept the filter or try another one
    loop='A';
    while (loop ~= 'Y') && (loop ~= 'N')
        loop=input('Are you satisfied? [Y/N] ','s');
        if loop=='Y'
            out=1;
            filter=filter_MY;
            r_max=param(1);
            searchSize=maxArea;
        elseif loop=='N'
            out=0;
        else
            disp('Incorrect input from user.')
            loop = 'A';
        end
    end
end