function [D,I,I_vec]=DLS_DistroFit(t,data,varargin)
% DLS_DistroFit Estimation of Intensity Size Distribution from DLS autocorrelation
% data
%
% [D,I,I_vec]=DLS_DistroFit(t,data) calculates the intensity size distribution 
%  parameters I of a sample from the autocorrelation data obtained by a
%  DLS experiment, using a non-negative least squares fit scheme based on CONTIN.
%  The measurement setup is assumed to be a Malvern Zetasizer Nano APS.
%  
%  DATA can be a matrix containing several signals (one per column, all the same length).
%  I_vec is the set of all the calculated I.


% [D,I,I_vec]=DLS_CumFit(t,data,'PropertyName',PropertyValue)
%  permits to set the value of PropertyName to PropertyValue.
%  Admissible Properties are:
%       T       -   Temperature (default = 298 K)
%       eta     -   Solvent viscosity (default = 0.8872 cP)
%       n       -   Solvent refractive index (default = 1.33)


% CREATED: Alessio Caciagli, University of Cambridge, October 2017

T = 298;
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'T')
        T = varargin{n+1};
    end
end
eta = 0.0008872;
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'eta')
        eta = varargin{n+1};
    end
end
n = 1.33;
for i = 1:2:length(varargin)
    if strcmpi(varargin{i},'n')
        n = varargin{i+1};
    end
end
%%
%Constants
kB = 1.38*10^(-23);
lambda = 633*10^(-9);
theta = 173*2*pi/360;
q = 4*pi*n*sin(theta/2)/lambda;

alpha = 1e-2;

I_vec = zeros(100,size(data,2));


%%Main loop
for j = 1:size(data,2)
    
    %Delete first 2 points (usually noise)
    Time_ToFit = t(3:end);
    ACF_ToFit = data(3:end,j);
    
    %Least-Squares fit
    CritPts = find(ACF_ToFit< 0.1*ACF_ToFit(1));
    tTemp = Time_ToFit(1:CritPts(1));
    dataTempLog = log(ACF_ToFit(1:CritPts));
    A = [ones(length(tTemp),1),tTemp,tTemp.^2];
    Coeff = A\dataTempLog;
    
    %g1(tau) calculation
    beta = exp(Coeff(1));
    dataTemp = sqrt(ACF_ToFit(1:CritPts))/sqrt(beta);

    %CONTIN loop run
    coarse_s = logspace(0,5,10)';
    coarse_g = ones(size(coarse_s));
    [coarse_g,~,~] = rilt(tTemp,dataTemp,coarse_s,coarse_g,alpha,'logarithmic',[],[],[],{'g>0'},[],[]);
    for i = 2:10
        s = logspace(0,5,i*10)';
        g0 = interp1(coarse_s,coarse_g,s,'linear');
        [g,~,~] = rilt(tTemp,dataTemp,s,g0,alpha,'logarithmic',[],[],[],{'g>0'},[],[]);
        coarse_g = g;
        coarse_s = s;
    end
    I_vec(:,j) = g;
end
close all

%D (in um) and I (in unity)
D = s.*q^2*kB*T./(3*pi*eta);
I = mean(I_vec,2);
I = I./sum(I);

