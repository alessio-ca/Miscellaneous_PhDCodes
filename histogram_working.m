 function h = histogram(h,data)
    
      % If this is the first time this histogram is being accessed,
      % initialize all the properties of this instacne
    
    if (h.count == 0)
        % Determine the number of bins by evaluating the histogram''s range and increment size  
        nBins = ceil((h.range(2)-h.range(1))/h.increment);
        
        % Adjust the histogram''s range
        % Useful if the total range is not an exact multiple of the increment size
        h.range(2) = h.range(1) + nBins * h.increment;
     
        % Set all bins to zero
        h.histo = zeros(1,nBins);
    
        % Set the values vector to be in the center of each bin
        h.values = 1:nBins;
        h.values = h.range(1) + h.increment*(h.values-0.5);
    end
    
    % Now that the histogram is initialized, add the data in the right bin
    
    if (data > h.range(1) && data <= h.range(2)) % Make sure the data fits the range
        % Find the right bin position
        binIndex = ceil((data-h.range(1))/h.increment);
        
        % Add 1 to the bin
        h.histo(binIndex) = h.histo(binIndex)+1;
        
        % Increment the count by 1
        h.count = h.count+1;
    
    else
       % Display error message
        disp(strvcat('histogram- Value out of range:',num2str(data)));
        return
    end
    
 end