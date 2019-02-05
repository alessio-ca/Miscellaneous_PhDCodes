function [msd,tau,msd_err,msd_vec] = msd_routine(dt,Vx,varargin)
% msd_routine   Mean Square Displacement
%
% [MSD,TAU,MSDerr,MSDvec] = msd_routine(DT,VX) calculates the mean square displacement MSD
%   of signal VX sampled at time intervals DT.
%   VX can be a matrix containing several signals (one per row, all the same length).
%   TAU are the values of the time delays.
%   MSDerr is the error calculated as standard deviation of the MSD value.
%   MSDvec is the set of all MSD calculated
%
% [MSD,TAU,MSDerr,MSDvec] = msd_routine(DT,VX,'PropertyName',PropertyValue) permits
%   to set the value of PropertyName to PropertyValue.
%   Admissible Properties are:
%       TauMax      -   Maximum delay to be calcualted (default = +Inf)

% CREATED: Alessio Caciagli, University of Cambridge, February 2017


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

% Analysis

%First component, S1
D = zeros(size(Vx,1)+2,size(Vx,2));
D(2:size(Vx,1)+1,:) = (Vx-repmat(mean(Vx,1),[size(Vx,1),1])).^2;
Q = 2*sum(D,1);

S1 = zeros(maxlags,size(Vx,2));
for n = 1:1:maxlags
    Q = Q - D(n,:) - D(size(Vx,1)+3-n,:);
    S1(n,:) = Q;
end

%Second component, S2 (autocorrelation function)
tmp = fft(padarray(Vx-repmat(mean(Vx,1),[size(Vx,1),1]),[size(Vx,1),0],'post'));
tmp = real(ifft(tmp.*conj(tmp)));
S2 = tmp(1:1:maxlags,:);

%Total, sum of terms (with unbiased normalization)
scale = size(Vx,1)*ones(maxlags,1) - (0:maxlags-1)';    
msd_vec = (S1 - 2*S2)./repmat(scale,[1,size(Vx,2)]);
msd_vec(1,:) = 0; %Correct round-off errors on first point

tau = dt*(0:1:maxlags-1)';
msd = mean(msd_vec,2);
msd_err = std(msd_vec,0,2);