function [Z,PDI,Z_vec,PDI_vec,Beta_vec]=DLS_CumFit(t,data,varargin)
% DLS_CumFit Estimation of Z-Ave and PDI parameters from DLS autocorrelation
% data
%
% [Z,PDI,Z_vec,PDI_vec,Beta_vec]=DLS_CumFit(t,data) calculates the size distribution parameters
%  (Z and PDI) of a sample from the autocorrelation data obtained by a
%  DLS experiment, using the cumulants fit scheme as in the ISO 22412
%  (2008). The measurement setup is assumed to be a Malvern Zetasizer Nano
%  APS.
%  
%  DATA can be a matrix containing several signals (one per column, all the same length).
%  Z_vec is the set of all the calculated Z.
%  PDI_vec is the set of all the calculated PDI.
%  BETA is the y-intercept of the autocorrelation function (obtained by the fit)
%  Beta_vec is the set of all the calculated Beta.

% [Z,PDI,Z_vec,PDI_vec,Beta_vec]=DLS_CumFit(t,data,'PropertyName',PropertyValue)
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


%Delete first 2 points (usually noise)
Time_ToFit = t(3:end);
ACF_ToFit = data(3:end,:);

%%Main loop
CritPts = zeros(1,size(ACF_ToFit,2));
Z_vec = zeros(1,size(ACF_ToFit,2));
PDI_vec = zeros(1,size(ACF_ToFit,2));
Beta_vec = zeros(1,size(ACF_ToFit,2));

for i = 1:size(ACF_ToFit,2)
    %Find far end for fit (<0.1*intercept)
    CritPts = find(ACF_ToFit(:,i) < 0.1*ACF_ToFit(1,i));
    CritPts = CritPts(1);
    %Linearize for fit
    TimeLog = Time_ToFit(1:CritPts);
    ACFLog = log(ACF_ToFit(1:CritPts,i));

    %Least-Squares fit
    A = [ones(length(TimeLog),1),TimeLog,TimeLog.^2];
    Coeff = A\ACFLog;

    %Data extraction
    Beta_vec(i) = exp(Coeff(1));
    Gamma = -Coeff(2)/2;
    k2 = Coeff(3);
    D = Gamma/q^2;
    Z_vec(i) = kB*T/(3*pi*eta*D);
    PDI_vec(i) = k2/Gamma^2;
end

Z = mean(Z_vec);
PDI = mean(PDI_vec);







