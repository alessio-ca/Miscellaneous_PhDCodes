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
                gR.range = [0 0.3*L];
                gR.increment = 0.3*L/NumOfBins;
                gR.saveFileName = 'radialDist.dat';
                gR.nPartTot = 0;


                
            case sample
                % Loop over pairs and determine the distribution of distances
                nPart = size(coords,1);
                for partA = 1:(nPart-1)
                    for partB = (partA+1):nPart
                        % Calculate particle-particle distance
                        % Account for PBC (assuming 2D)                               
                        dr = coords(partA,:) - coords(partB,:); 
                        dr = distPBC2D(dr,Lx,Ly);
                        % Get the size of this distance vector
                        r = sqrt(sum(dot(dr,dr)));
                        
                        % Add to g(r) if r is in the right range [0 0.3*L]
                        if (r < 0.3*L)
                            gR = histogram(gR,r);
                        end
                    end
                end
%                 %Weighted mean
%                 gR.histoTot = gR.histoTot + nPart*gR.histo;
%                 %Weight sum
%                 gR.nPartTot = gR.nPartTot + nPart;
%                 gR.count = gR.count+1;

                
            case results
              
                % each bin should be normalized according to its volume
                gR.histo = gR.histo/(gR.count*gR.increment);
                nBins = size(gR.values,2);
                nPart = size(coords,1);
                rho = nPart/(Lx*Ly); % Density of the particles
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
                end

            case plot
                
                plot(gR.values,gR.histo);
                xlabel('distance');
                ylabel('g(R)');
                
            otherwise
                % Wrong switch
                disp('radialDistribution : You have entered an illegal switch value');
        end
    
    end
    