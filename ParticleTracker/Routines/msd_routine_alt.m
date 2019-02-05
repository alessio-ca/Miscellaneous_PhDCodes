function [msd,tau,msd_err,msd_vec] = msd_routine_alt(dt,Vx,varargin)
% msd_routine   Mean Square Displacement
%
% [MSD,TAU,MSDerr,MSDvec] = msd_routine_alt(DT,VX) calculates the mean square displacement MSD
%   of signal VX sampled at time intervals DT. It can correct some
%   shortcomings of the official routine (e.g. bias in the results due to
%   correlated data at large times), but it is MUCH slower (scales with O(N^2)).
%   VX can be a matrix containing several signals (one per row, all the same length).
%   TAU are the values of the time delays.
%   MSDerr is the error calculated as standard deviation of the MSD value.
%   MSDvec is the set of all MSD calculated
%
% [MSD,TAU,MSDerr,MSDvec] = msd_routine_alt(DT,VX,'PropertyName',PropertyValue) permits
%   to set the value of PropertyName to PropertyValue.
%   Admissible Properties are:
%       TauMax      -   Maximum delay to be calcualted (default = +Inf)

% CREATED: Alessio Caciagli, University of Cambridge, January 2017


% Maximum delay
taumax = +Inf;
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'taumax')
        taumax = varargin{n+1};
    end
end
if isfinite(taumax)
    maxlags = ceil(taumax/dt);
else
    maxlags = size(Vx,1);
end

% Analysis (parallel is enabled only if sample size is larger than 10000)
msd_vec = zeros(maxlags,size(Vx,2));
if (size(Vx,1) > 1e4)
    parfor m = 1:1:maxlags-1
        msd_vec(m+1,:) = mean((Vx(m+1:end,:)-Vx(1:end-m,:)).^2 , 1);
    end
else
    for m = 1:1:maxlags-1
        msd_vec(m+1,:) = mean((Vx(m+1:end,:)-Vx(1:end-m,:)).^2 , 1);
    end
end

tau = dt*(0:1:maxlags-1)';

msd = mean(msd_vec,2);
msd_err = std(msd_vec,0,2);