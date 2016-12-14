function MSD_list = MSD_single_track(trackdata)
% USAGE:    msd = MSD_single_track(trackdata)
% PURPOSE:
%           Calcuate the MSD from single particle trajectory
%
%
% INPUT:
% trackdata
%           Trajectory data of the 4-column format:
%           x-position     y-position       frame#      bead#
%
%
% OUTPUT:
% MSD: will have three columns:
%       Column 1: lag time (in frames)
%       Column 2: MSD (in pixels)
%       Column 3: number of observations in average
%       Column 4: std error of MSD (in pixels)
%
% CREATED: Eric M. Furst, University of Delaware, September 8, 2013
%  Modifications:


% Step through all particles in the tracked data set
% The maximum frame separation is totalframe - 1
totalframe = length(trackdata);
max_step = totalframe - 1;
x=trackdata(:,1);
y=trackdata(:,2);

% Only analyze if there is more than one frame to calculate
% displacement
if max_step >= 1
    disp(['Total frames: ',...
        num2str(totalframe)]);
    % Step through all frame separations starting from 1 up to max_step
    t=zeros(max_step,1);
    mu=zeros(max_step,1);
    n=zeros(max_step,1);
    s=zeros(max_step,1);
    
    %Enable progress bar for parallel pool
    try
        parpool;
    catch ME
       if ~strcmp(ME.identifier,'parallel:convenience:ConnectionOpen')
           rethrow(ME)
       end
    end
    warning('off','MATLAB:Java:DuplicateClass')
    pctRunOnAll javaaddpath java
    progressStep=10000;
    ppm = ParforProgMon('Progress: ', max_step, progressStep, 300, 80);
    
    parfor step=1:max_step
        %deltacoords=trackdata(1+step:end,1:2)-trackdata(1:end-step,1:2);
        deltacoordsx=x(1+step:end)-x(1:end-step);
        deltacoordsy=y(1+step:end)-y(1:end-step);
        squaredisp=deltacoordsx.^2+deltacoordsy.^2;
        if(length(squaredisp)>1)
            t(step) = step;
            mu(step) = mean(squaredisp); %# average
            n(step) = length(squaredisp); %# n
            s(step) = std(squaredisp)/sqrt(n(step)); %# std error of the mean (defined as std/sqrt(N))
        end
        if mod(step,progressStep)==0
            ppm.increment();
        end
    end
    ppm.delete()
    
    t=t(t~=0);
    MSD_list = [t,mu(1:length(t)),n(1:length(t)),s(1:length(t))];
end