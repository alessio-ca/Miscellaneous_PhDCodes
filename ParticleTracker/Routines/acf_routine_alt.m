function [acf,tau,acf_err,acf_vec] = acf_routine_alt(dt,Vx,varargin)
% acf_routine   Autocorrelation function
%
% [ACF,TAU,ACFerr,ACFvec] = acf_routine_alt(DT,VX) calculates the autocorrelation function ACF
%   of signal VX sampled at time intervals DT. It can correct some
%   shortcomings of the official routine (e.g. bias in the results due to
%   correlated data at large times), but it is MUCH slower (scales with O(N^2)).
%   VX can be a matrix containing several signals (one per row, all the same length).
%   TAU are the values of the time delays.
%   ACFerr is the error calculated as standard deviation of the ACF value.
%   ACFvec is the set of all ACF calculated
%
% [ACF,TAU,ACFerr,ACFvec] = acf_routine_alt(DT,VX,'PropertyName',PropertyValue) permits
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
acf_vec = zeros(maxlags,size(Vx,2));
mean_vec = repmat(mean(Vx,1),[size(Vx,1),1]);
if (size(Vx,1)*size(Vx,2) > 1e4)
    parfor i = 1:maxlags
        acf_vec(i,:) = mean((Vx(1:end-i+1,:) - mean_vec(1:end-i+1,:)).*(Vx(i:end,:)-mean_vec(i:end,:)),1);
    end
else
    for i = 1:maxlags
        acf_vec(i,:) = mean((Vx(1:end-i+1,:) - mean_vec(1:end-i+1,:)).*(Vx(i:end,:)-mean_vec(i:end,:)),1);
    end
end
tau = dt*(0:1:maxlags-1)';

acf = mean(acf_vec,2);
%Normalization
acf = acf./acf(1); 
acf_err = std(acf_vec,0,2);
