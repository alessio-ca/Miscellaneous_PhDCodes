function [filter,filter_pres,r_max]=CreateFilter(image)
% USAGE:    CreateFilter
% PURPOSE:  Overlays a mask over an image for empirically checking whether
%           the mask is appropriate
%           
% 
% INPUT: Image (assumed to be already loaded via imread)
%
% OUTPUT: Filter mask appropriate for current frame
%
% CREATED: Alessio Caciagli, University of Cambridge, 05/12/2015
% MODIFIED: Alessio Caciagli, University of Cambridge, 21/12/2016

out=0;
figure(1)
imshow(mat2gray(image))
drawnow
commandwindow
while out == 0
    param=[-1,-1,-1];
    %Set security commands for wrong user input
    while any(param(1:3) < 0) || any(diff(param(1:3))>0) || any(mod(param(1:3),2) == 1)
        param=input('Insert filter parameters [d_max d_mid d_min] ');
        if (length(param) <= 3)
            filter_MY=create_filter_3(param(1),param(2),param(3)); %Create filter with specified parameters
        else
            filter_MY=create_filter_3(param(1),param(2),param(3),param(4),param(5),param(6)); %Create filter with specified parameters
        end
    end
    X_Lim=get(gca,'XLim'); %Commands to keep the zoom in imshow during filter creation
    Y_Lim=get(gca,'YLim');
    
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
    
    if (ypeak-length(filter_MY)+1 <= 0) || (xpeak-length(filter_MY)+1 <= 0)
        disp('Error in correlation matrix maximum. Try redefining the mask.')
        continue
    end
    
    %Manipulation for ShowFilter (which automatically calculates the offset)
    xpeak=xpeak-floor(size(filter_MY,2)/2); 
    ypeak=ypeak-floor(size(filter_MY,1)/2);
    
    if (length(param) <= 3)
        filter_pres=create_filter_3(param(1),param(2),param(3),0.5,0,1); %Create filter with specified parameters for visualization
    else
        A=[0,0.5,1];
        [~,I]=sort(param(4:6));
        filter_pres=create_filter_3(param(1),param(2),param(3),A(I==1),A(I==2),A(I==3));
        
    end
   
    
    
    %Refining of filter appearence in imshow as a green mask
    ShowFilter(mat2gray(image),filter_pres,[xpeak,ypeak],[0 1 0])
    drawnow
    commandwindow
    if (peakint < 0.6)
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
        elseif loop=='N'
            out=0;
        else
            disp('Incorrect input from user.')
            loop = 'A';
        end
    end
end