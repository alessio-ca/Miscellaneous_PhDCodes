function [acf,tau,acf_err,acf_vec] = acf_routine(dt,Vx,varargin)
% acf_routine   Autocorrelation function
%
% [ACF,TAU,ACFerr,ACFvec] = acf_routine(DT,VX) calculates the autocorrelation function ACF
%   of signal VX sampled at time intervals DT.
%   VX can be a matrix containing several signals (one per row, all the same length).
%   TAU are the values of the time delays.
%   ACFerr is the error calculated as standard deviation of the ACF value.
%   ACFvec is the set of all ACF calculated
%
% [ACF,TAU,ACFerr,ACFvec] = acf_routine(DT,VX,'PropertyName',PropertyValue) permits
%   to set the value of PropertyName to PropertyValue.
%   Admissible Properties are:
%       TauMax      -   Maximum delay to be calcualted (default = +Inf)

% CREATED: Alessio Caciagli, University of Cambridge, January 2017


% Maximum delay & normalization
norm = 1;
taumax = +Inf;
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'taumax')
        taumax = varargin{n+1};
    elseif strcmpi(varargin{n},'norm')
        norm = varargin{n+1};
        if norm ~=0 && norm ~= 1
            error('Wrong variable type (expected logical 0 or 1)')
        end
    end
end
if isfinite(taumax)
    maxlags = ceil(taumax/dt);
else
    maxlags = size(Vx,1);
end

% Analysis

%FFT method based on zero padding & unbiased scaling
%tmp = fft(padarray(Vx-repmat(mean(Vx,1),[size(Vx,1),1]),[size(Vx,1),0],'post'));
tmp = fft(padarray(Vx,[size(Vx,1),0],'post'));
tmp = real(ifft(tmp.*conj(tmp)));
scale = size(Vx,1)*ones(maxlags,1) - (0:maxlags-1)';    
acf_vec = tmp(1:1:maxlags,:)./repmat(scale,[1,size(Vx,2)]);
    

%Normalization
tau = dt*(0:1:maxlags-1)';
%acf_vec = acf_vec./repmat(acf_vec(1,:),[size(acf_vec,1),1]);
if norm == 1
    acf_vec=acf_vec./repmat(acf_vec(1,:),size(Vx,1),1);
end
acf = mean(acf_vec,2);
acf_err = std(acf_vec,0,2);
