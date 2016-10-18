    function plotCircle(center,radius,ptsNumber,color)
    
    
    %
    % Set default values
    %
    if nargin < 4    % If the user did not supply the 4th value
        color = 'k'; % Set it to default
    end
    
    if nargin < 3       % If the user did not supply the 4th value
        ptsNumber = 20; % Set it to default
    end
    
    
    % Make sure the to draw on the existing figure
    hold on;
    
    % Define the set of angles points should appear in
    angles = 0:(2*pi/ptsNumber):2*pi;
    
    % Define the array of x and y for those angles
    % defined by the center point and the radius given
    x = center(1) + radius*cos(angles);
    y = center(2) + radius*sin(angles);
    
    % Plot the circle
    plot(x,y,color);
    
    end
    

