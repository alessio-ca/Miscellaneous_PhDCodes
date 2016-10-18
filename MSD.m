function msd = MSD(trackdata)
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

MSD_list_tot = [];

% Step through all particles in the tracked data set
for particleid=1:max(trackdata(:,4))
    
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
            MSD_list=zeros(max_step,1);
        % Step through all frame separations starting from 1 up to max_step
        for step=1:max_step
            deltacoords=trackdata(index+step:index+max_step,1:2)-trackdata(index:index+max_step-step,1:2);
            squaredisp=sum(deltacoords.^2,2);
            if(length(squaredisp)>1)
            MSD_list(step,1) = step; 
            MSD_list(step,2) = mean(squaredisp); %# average
            MSD_list(step,3) = length(squaredisp); %# n
            MSD_list(step,4) = std(squaredisp); %# std
            end
        end
    end
    MSD_list_tot=[MSD_list_tot;MSD_list];
end
% Build the MSD from the MSD_list
% the MSD will have three columns:
%   Column 1: lag time (in frames)
%   Column 2: MSD (in pixels)
%   Column 3: number of observations in average

msd = [];
min_lag = min(MSD_list_tot(:,1));
max_lag = max(MSD_list_tot(:,1));
for lag=min_lag:max_lag
    ind = find(MSD_list_tot(:,1)==lag & MSD_list_tot(:,4)>0);
    number_obs=sum(MSD_list_tot(ind,3));
    if(number_obs>1)
        weight_sum = sum(MSD_list_tot(ind,3));
        %weight_sum = sum(1./MSD_list_tot(ind,4).^2);
        mean_MSD=sum(MSD_list_tot(ind,2).*MSD_list_tot(ind,3))./weight_sum;
        %mean_MSD=sum(MSD_list_tot(ind,2)./(MSD_list_tot(ind,4).^2))./weight_sum;
        std_MSD=sqrt(sum(MSD_list_tot(ind,2).*MSD_list_tot(ind,2).*MSD_list_tot(ind,3))./weight_sum - mean_MSD.^2);
        msd = [msd, [lag mean_MSD number_obs std_MSD]'];
    end
end

msd = msd';

% Run with 
% >>msd = MSD(result)
% Plot with 
% >>loglog(msd(:,1),msd(:,2:3))


end

