function [psd,f,psd_err,psd_vec] = psd_routine(dt,Vx,varargin)
% PSD   Power spectral density
%
% [PSD,F,PSDerr,PSDvec,Fvec] = PSD(DT,VX) calculates the power spectral density PSD
%   of signal VX sampled at time intervals DT.
%   VX can be a matrix containing several signals (one per row, all the same length).
%   F are the values of the frequencies.
%   PSDerr is the PSD error.
%   PSDvec is the set of all PSD calculated
%
% [PSD,F,PSDerr,PSDvec,Fvec] = PSD(DT,VX,'PropertyName',PropertyValue) permits
%   to set the value of PropertyName to PropertyValue.
%   Admissible Properties are:
%       Fmin            -   minimum frequency (default = 1/total acquisition time)
%       Fmax            -   maximum frequency (defauls = 1/DT)
%       RunAvg          -   running average (default = 'none' | 'lin' | 'log')
%       NumBins         -   number of bins for the running average (default = 100)

% CREATED: Alessio Caciagli, University of Cambridge, February 2017


% Minimum frequency
fmin = 1/(dt*size(Vx,1));
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'fmin')
        fmin = varargin{n+1};
    end
end

% Maximum frequency
fmax = 1/dt;
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'fmax')
        fmax = varargin{n+1};
    end
end

% Blocking: none (default) | lin | log
RunAvg = 'none';
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'RunAvg')
        RunAvg = varargin{n+1};
    end
end

% Number of bins for blocking (N.B ! It should be sufficiently large so
% that central limit theorem assures the blocked Pk are gaussian
% distributed)
NumBins = 100;
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'NumBins')
        NumBins = varargin{n+1};
    end
end

% Analysis
fs = 1/dt;
T = dt*size(Vx,1); % measurement time

ft_vec = dt*fft(Vx);
psd_vec = ft_vec.*conj(ft_vec)/T;

fNyq = fs/2; % Nyquist frequency
f = (1:1:size(Vx))'/T;

% Only data up to the Nyquist frequency
% and between fmin and fmax
ind = find(f<=fNyq & f>=fmin & f<=fmax);
f = f(ind);
%Matlab indexing: disregard 0 component & take up to Nyquist
psd_vec = psd_vec(ind+1,:);

psd = mean(psd_vec,2);
psd_err = std(psd_vec,0,2);

% Blocking
if strcmpi(RunAvg,'lin')
    tic
    %Window size
    WinSize=round(length(f)/NumBins);
    % Centers of bins
    ind_min=ind(ind-floor(WinSize)>0);
    ind_min=ind_min(1);
    
    ind_max=ind(ind+ceil(WinSize)-length(f)<0);
    ind_max=ind_max(end);
    
    binscenters = [1,(ind_min:WinSize:ind_max),length(f)];
    
    psd_b = movmean(psd,WinSize);
    psd_err_b = movstd(psd,WinSize);
    f_b = movmean(f,WinSize);
    
    f = f_b(binscenters);
    psd = psd_b(binscenters);
    psd_err = psd_err_b(binscenters)/sqrt(WinSize/2);
    psd_err(2:end-1)=psd_err(2:end-1)/sqrt(2);
    toc
elseif strcmpi(RunAvg,'log')

    %Window size
    WinSize=round(length(f)/NumBins);
    % Centers of bins
    ind_min=ind(ind-floor(WinSize/2)>0);
    ind_min=ind_min(1);
    
    ind_max=ind(ind+ceil(WinSize/2)-length(f)<0);
    ind_max=ind_max(end);
    
    binscenters = (log(ind_min):log(length(f))/NumBins:log(ind_max));
    binscenters = unique(floor(exp(binscenters)));
    
    psd_b = movmean(psd,WinSize);
    psd_err_b = movstd(psd,WinSize);
    f_b = movmean(f,WinSize);
    
    f = f_b(binscenters);
    psd = psd_b(binscenters);
    psd_err = psd_err_b(binscenters)/sqrt(WinSize);

end