function gR = radialDistribution2D(switchVal,gR,coords,Lx,Ly,NumOfBins)
        % before running code create gR with: 
        % gR = struct;
        
%         input:
%         ~~~~~~
%         switchVal - operation of the function
%         initialize=0;
%         sample=1;
%         results=2;
%         plot = 3;
%           
%         gR - the output struct contaning the histogram
%         (should be created with "gR = struct;" before initialization
%
%         coords - x,y coordinates of particle centers
%         first row x, second row y
%         
%         Lx - x size of area to calculate
%         Ly - y size of area to calculate
%         
%         NumOfBins - number of bins in the final histogram
%
%         output:
%         ~~~~~~~
%         gR is a histogram struct contaning:
%                 gR.count - total number of object counted
%                 gR.range - x axis range for the histogram
%                 gR.increment - increment of x axis  
%                 gR.saveFileName - name of saved dat file
%   

%       method:
%       ~~~~~~~
%       the RDF is calculated by binnig all pair partical distances into 
%       a histogram, and normalizing each bin with it's Ideal gas number of
%       particles. 
%       when binning the pair distances, we take into account Periodic
%       Boudary Conditions (PBC)
%       finally, to ensure that when r->infinity : RDF ->1 , we
%       normalize by multiplying RDF*2/(N-1) where N is the number of
%       particles. 
%       for more information http://www2.msm.ctw.utwente.nl/sluding/TEACHING/APiE_Script_v2011.pdf
%       page 48 - "Radial distribution function"
%       
%       -------
%     |important|
%       -------
%       this function uses the functions "histogram", "distPBC2D"
%       "histogram" can be found here: http://www.cchem.berkeley.edu/chem195/histogram_8m.html#aedd379efd57ae78820ad8787bfab0cce
%       "distPBC2D" can be found here: http://www.mathworks.com/matlabcentral/fileexchange/46575-distpbc2d-m
%       
%       


        L = sqrt(Lx^2+Ly^2);
        % Set three operation cases for this function
        initialize=0;
        sample=1;
        results=2;
        plot = 3;
    
        % Choose the operation method according to the providerd switchVal
        switch switchVal
            
            case initialize
                % Initialize a histogram to hold the radial distribution
                
                
                gR.count = 0;
                gR.countWeight = 0;
                gR.range = [0 0.3*L];
                gR.increment = 0.3*L/NumOfBins;
                gR.saveFileName = 'radialDist.dat';
                gR.nBins = ceil((gR.range(2)-gR.range(1))/gR.increment);
                gR.nPartTot = 0;
                gR.histoTot = zeros(1,gR.nBins);


                
            case sample
                % Loop over pairs and determine the distribution of distances
                gR.histoAlt = zeros(1,gR.nBins);
                nPart = size(coords,1);
                for partA = 1:(nPart-1)
                    for partB = (partA+1):nPart
                        % Calculate particle-particle distance
                        % Account for PBC (assuming 2D)                               
                        dr = coords(partA,:) - coords(partB,:);                         
                        
                        % Calculate the half box size in each direction
                        hLx = Lx/2.0;
                        hLy = Ly/2.0;
                        
                        % Distance vector should be in the range -hLx -> hLx and -hLy -> hLy
                        % Therefore, we need to apply the following changes if it's not in this range:
                        if dr(1) > hLx
                            dr(1) = dr(1) - Lx;
                        elseif dr(1) < -hLx
                            dr(1) = dr(1) + Lx;
                        end
                        
                        if dr(2) > hLy
                            dr(2) = dr(2) - Ly;
                        elseif dr(2) < -hLy
                            dr(2) = dr(2) + Ly;
                        end
                        
                        
                        % Get the size of this distance vector
                        r = sqrt(dr*dr');
                        
                        % Add to g(r) if r is in the right range [0 0.3*L]
                        if (r < 0.3*L)
                            gR = histogram_GR(gR,r);
                        end
                    end
                end
                %Weighted mean
                gR.histoTot = gR.histoTot + nPart*gR.histoAlt;
                %Weight sum
                gR.nPartTot = gR.nPartTot + nPart;
                gR.countWeight = gR.countWeight + 1;
                
                
            case results
                
                % each bin should be normalized according to its volume
                nBins = size(gR.values,2);
                nPart = size(coords,1);
                rho = gR.nPartTot/(Lx*Ly*gR.countWeight); % Density of the particles
                disp(rho)
                
                for bin = 1:nBins
                    % rVal is the number of cells in some layer of area
                    % da(r)=2 pi * r * dr, a distance r from the central cell
                    rVal = gR.values(bin);
                    next_rVal = gR.increment + rVal;
                    % Calculate the area of the bin (a ring of radii r,
                    % r+dr)
                    AreaBin = pi*next_rVal^2 - pi*rVal^2;
                    % Calculate the number of particles expected in this bin in
                    % the ideal case
                    nIdeal = AreaBin*rho;
                    % Normalize the bin
                    gR.histo(bin) = 2*gR.histo(bin) / nIdeal;
                    gR.histoTot(bin) = 2*gR.histoTot(bin) / nIdeal;
                end
                
                %Normalize the histogram
                gR.histo = gR.histo/(gR.increment*gR.count);
                gR.histoTot = 2*gR.countWeight*gR.histoTot/(gR.increment*gR.count*gR.nPartTot);
                
                
            case plot
                
                plot(gR.values*0.1,gR.histo);
                xlabel('r (\mum)');
                ylabel('g(R)');
                
            otherwise
                % Wrong switch
                disp('radialDistribution : You have entered an illegal switch value');
        end
        
end
