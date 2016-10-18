    function drawConfig(coords,Lx,Ly,radius)
    
        %
        % Set default values
        %
        if nargin < 4 
            radius = 0.25; % set a defualt value for particle radius
        end
    
        % Prepare the figure
        cla;
        hold on;
    
        % Draw the box outline
        boxX = [0 Lx Lx 0 0];
        boxY = [0 0 Ly Ly 0];
        plot(boxX,boxY);
    
    
        % Find the number of particles
        nPart = size(coords,2);
    
        % Plot the particles 
        for i=1:nPart
            plotCircle(coords(:,i),radius);
        end
    
        % Some more figure adjustments
        axis equal;
        axis([-0.5 Lx+0.5 -0.5 Ly+0.5]); % Set the lower and upper limits of the axes
    
        drawnow;
        figure(gcf);
    end
    

