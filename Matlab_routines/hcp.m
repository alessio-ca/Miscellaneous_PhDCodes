    function [coords Lx Ly] = hcp(nPart,density,radius)
    
    % Make sure nPart is a perfect square
    intRoot = floor(sqrt(nPart));
    if (sqrt(nPart) - intRoot) > 1e-7
        % Display an error message
        disp('Number of particles should be a perfect square');
        coords = zeros(2,nPart);
        Lx = 0.0;
        Ly = 0.0;
        return 
    end
    
    % Calculate the seperation length between disk centers
    sepDist = sqrt(2*pi*radius*radius/(density*sqrt(3)));
    
    % Find the box size
    Lx = sepDist * sqrt(nPart);
    Ly = Lx*sqrt(3)/2;
    
    % Create a vector of linearly spaced points along the x-direction
    xPos = linspace(sepDist/2, Lx-sepDist/2, sqrt(nPart));
    % And find the corresponsing y-direction increments
    yPos = (sqrt(3)/2)*xPos;
    
    % Create a matrix with all combinations of x and y
    [X,Y] = meshgrid(xPos,yPos);
    % Shift coordinates to the be at the center of each particle
    X(1:2:end,:) = X(1:2:end,:) + sepDist/2;
    
    % Reshape the matrix to be 1D in X and 1D in Y
    % (numel returns the number of elements in a given array)
    coords = [reshape(X,1,numel(X));reshape(Y,1,numel(Y))];
    
    end