function msd = MSD_single_track(trackdata)
% USAGE:    msd = MSD_single_track(trackdata)
% PURPOSE:  
%           Calcuate the MSD from single particle trajectory
%           
% 
% INPUT:
% trackdata 
%           Trajectory data of the 3-column format:
%           frame      x-position     y-position
%           
%   
% OUTPUT:  
% MSD: will have three columns:
%       Column 1: lag time (in frames)
%       Column 2: MSD (in pixels)
%       Column 3: number of observations in average
%       Column 4: std of MSD (in pixels)
%
% CREATED: Eric M. Furst, University of Delaware, September 8, 2013
%  Modifications:

MSD_list = [];

% Step through all particles in the tracked data set
% The maximum frame separation is totalframe - 1
    totalframe = length(trackdata);
    max_step = totalframe - 1;
    index=1;
    
    % Only analyze if there is more than one frame to calculate
    % displacement
    if max_step >= 1
          disp(['Total frames: ',...
                                num2str(totalframe)]);
        % Step through all frame separations starting from 1 up to max_step
        for step=1:max_step
            for j=1:(fix(max_step/step))
           % Caculate lag time and mean-squared displacement between steps
              delta_t = trackdata(index+j*step,1)...
                  - trackdata(index+(j-1)*step,1);
              delta_r2 = (trackdata(index+j*step,2)...
                  - trackdata(index+(j-1)*step,2))^2 ...
                  + (trackdata(index+j*step,3)...
                  -  trackdata(index+(j-1)*step,3))^2;
             MSD_list = [MSD_list, [delta_t,delta_r2]'];
            end
        end
    end

MSD_list = MSD_list';

% Build the MSD from the MSD_list
% the MSD will have three columns:
%   Column 1: lag time (in frames)
%   Column 2: MSD (in pixels)
%   Column 3: number of observations in average
%   Column 4: standard dev on the mean

msd = [];
min_lag = min(MSD_list(:,1));
max_lag = max(MSD_list(:,1));
for lag=min_lag:max_lag
    number_obs = sum(MSD_list(:,1)==lag);
    if(number_obs>=2)
        ind = find(MSD_list(:,1)==lag);
        mean_MSD=mean(MSD_list(ind,2));
        sum_mean=(MSD_list(ind,2)-mean_MSD);
        std_MSD=sqrt(sum(sum_mean.^2)./number_obs);
        %msd = [msd, [lag mean(MSD_list(ind,2)) number_obs]'];
        msd = [msd, [lag mean_MSD number_obs std_MSD]'];
    end
end

msd = msd';

end

