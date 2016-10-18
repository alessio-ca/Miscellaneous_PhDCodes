function msd = MSD_INDEPPOINTS(trackdata)
% USAGE:    msd = MSD(trackdata)
% PURPOSE:  
%           Calcuate the MSD from the output of track.m
%           
% 
% INPUT:
% trackdata 
%           Trajectory data of the 4-column format:
%           x-position     y-position    frame     particle i.d.
%           
% NOTES: The particle i.d. is assumed to increment from 1 to a maximum
% number of particles. For each particle i.d. the frame number is assumed
% to increase. Trajectories can include skipped frames.
%   
% OUTPUT:  
% MSD: will have three columns:
%       Column 1: lag time (in frames)
%       Column 2: MSD (in pixels)
%       Column 3: number of observations in average
%
% CREATED: Eric M. Furst, University of Delaware, September 8, 2013
%  Modifications:

MSD_list = [];

% Step through all particles in the tracked data set
for(particleid=1:max(trackdata(:,4)))
    
    % Find the starting row of the particleid in trackdata
    % This is the offset
    index=find(trackdata(:,4)==particleid,1);
    
    % Find the number of frames for the particle
    % NOTE that this is NOT equivalent to time steps of 1 if the data has
    % "skipped frames"
    totalframe = sum(trackdata(:,4)==particleid);
    
    % The maximum frame separation is totalframe - 1
    max_step = totalframe - 1;
    
    % Only analyze if there is more than one frame to calculate
    % displacement
    if max_step >= 1
        disp(['Particle ', num2str(particleid),' of ',...
          num2str(max(trackdata(:,4))),'. Total frames: ',...
                                                num2str(totalframe)]);
        % Step through all frame separations starting from 1 up to max_step
        for step=1:max_step
%        Remove comment for debugging        
%        disp(['Number of steps of length ', ...
%           num2str(step), ': ',num2str(fix(totalframe/step))]);
            for j=1:(fix(max_step/step))
           
           % Caculate lag time and mean-squared displacement between steps
              delta_t = trackdata(index+j*step,3)...
                  - trackdata(index+(j-1)*step,3);
              delta_r2 = (trackdata(index+j*step,1)...
                  - trackdata(index+(j-1)*step,1))^2 ...
                  + (trackdata(index+j*step,2)...
                  -  trackdata(index+(j-1)*step,2))^2;
%           Remove comment for debugging           
%           disp(['  Step ', num2str(j),' of ',...
%               num2str(fix(totalframe/step)), ' Delta t = ',...
%               num2str(delta_t), ' Delta r^2 = ',num2str(delta_r2)]);

             MSD_list = [MSD_list, [delta_t,delta_r2]'];
            end
        end
    end
end

MSD_list = MSD_list';


% Build the MSD from the MSD_list
% the MSD will have three columns:
%   Column 1: lag time (in frames)
%   Column 2: MSD (in pixels)
%   Column 3: number of observations in average

msd = [];
min_lag = min(MSD_list(:,1));
max_lag = max(MSD_list(:,1));
for lag=min_lag:max_lag
    number_obs = sum(MSD_list(:,1)==lag);
    if(number_obs>=2)
        ind = find(MSD_list(:,1)==lag);
        msd = [msd, [lag mean(MSD_list(ind,2)) number_obs]'];
    end
end

msd = msd';

% Run with 
% >>msd = MSD(result)
% Plot with 
% >>loglog(msd(:,1),msd(:,2:3))


end

