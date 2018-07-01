function [tau,MSD,ACF,MSD_vec,ACF_vec,I_vec]=DLS_CRIEAnalysis(t,data,varargin)
% DLS_CRIEAnalysis Estimation of g1 and MSD from DLS autocorrelation
% data using a Constrained Regularization of Integral Equations (CRIE)
% method
%
% [tau,MSD,ACF,MSD_vec,ACF_vec]=DLS_CRIEAnalysis(t,data) calculates the
%   autocorrelation function g1 and the MSD from the autocorrelation
%   function ICF obtained in a DLS experiment, using the inversion relation.
%   The measurement setup is assumed to be a ZetaSizer APS Nano (Malvern).
%
%  DATA can be a matrix containing several signals (one per column, all the same length).
%  MSD_vec is the set of all the calculated MSDs.
%  ACF_vec is the set of all the calculated ACFs.
%  I_vec is the set of all the calculated g(s) coefficients. 

% [tau,MSD,ACF,MSD_vec,ACF_vec]=DLS_CRIEAnalysis(t,data,'PropertyName',PropertyValue)
%  permits to set the value of PropertyName to PropertyValue.
%  Admissible Properties are:
%       T       -   Temperature (default = 298 K)
%       eta     -   Solvent viscosity (default = 0.8872 cP)
%       n       -   Solvent refractive index (default = 1.33)
%       R       -   Bead radius (default = 115 nm)
%       beta    -   Calculation of coherence factor (default = 1)
%       Ng      -   Number of grid points for CRIE routine (default = 161)
%       cut     -   ACF cut value for fit (0<cut<1) (default = 0.1)
%       tail    -   Fit times for tail characterization (default =
%                   [0,0] i.e. no tail correction)
%       batch   -   Switch on/off batch calculation (no plots)


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
R = 115e-9;
for i = 1:2:length(varargin)
    if strcmpi(varargin{i},'R')
        R = varargin{i+1};
    end
end
beta = 1;
for i = 1:2:length(varargin)
    if strcmpi(varargin{i},'beta')
        beta = varargin{i+1};
    end
end
Ng = 161;
for i = 1:2:length(varargin)
    if strcmpi(varargin{i},'Ng')
        Ng = varargin{i+1};
    end
end
cut = 0.1;
for i = 1:2:length(varargin)
    if strcmpi(varargin{i},'cut')
        cut = varargin{i+1};
    end
end
tail = [0,0];
for i = 1:2:length(varargin)
    if strcmpi(varargin{i},'tail')
        tail = varargin{i+1};
    end
end
batch = 1;
for i = 1:2:length(varargin)
    if strcmpi(varargin{i},'batch')
        batch = varargin{i+1};
    end
end
%% Constants
kB = 1.38*10^(-23);
lambda = 633*10^(-9);
theta = 173*2*pi/360;
q = 4*pi*n*sin(theta/2)/lambda;

%% g1(tau) calc. and fit to g1(tau) (for smoothness)

I_vec = ones(Ng,size(data,2));
for j = 1:size(data,2)
    
    %Delete first 6 points (usually noise) & format data
    Time_ToFit = t(7:end);
    ACF_ToFit = data(7:end,j);
    if all(tail)
        tailFit = ACF_ToFit(Time_ToFit > tail(1) & Time_ToFit < tail(2));
        tailFit = mean(tailFit);
        ACF_ToFit = ACF_ToFit - tailFit;
    end
    CritPts = find(ACF_ToFit< cut*ACF_ToFit(1));
    tTemp = Time_ToFit(1:CritPts(1));
    
    %g1(tau) calculation with  optional head treatment (least-Squares)
    switch beta
        case 1
            dataTempLog = log(ACF_ToFit(1:10));
            A = [ones(length(tTemp(1:10)),1),tTemp(1:10),tTemp(1:10).^2];
            Coeff = A\dataTempLog;
            beta = exp(Coeff(1));
            
        case 0
            beta = 1;
            
        otherwise
            error('Error: Unrecognized value for beta.');
            
    end
    dataTemp = sqrt(ACF_ToFit(1:CritPts(1)))/sqrt(beta);
    
    if batch == 1
        %Plotting section
        f=figure(1);
        semilogx(t,data(:,1),'o')
        hold on
        semilogx(tTemp,dataTemp,'s')
        semilogx(t,cut*ones(length(t),1),'--k');
        semilogx(t,ones(length(t),1),'-k');
        legend('g2 (Raw)','g1 (Raw)','Location','northeast')
        xlabel('Time (us)')
        ylabel('ACF')
        hold off
        drawnow
        
        disp(' ')
        disp('Check your g2: the initial point MUST be < 1! If you are satisfied...')
        disp('Press the RETURN key to GO ON (make sure the Figure is active!).')
        disp('Press any other key to STOP.')
        
        
        pause; % wait for a keypress
        k = waitforbuttonpress;
        currkey=get(gcf,'CurrentKey');
        if strcmp(currkey, 'return')
            currkey=1;
        else
            return;
        end
    end
    
    %Input parameters for CRIE run
    DLS_Input.Iquad = 2;
    DLS_Input.Igrid = 2;
    DLS_Input.g_lims = [1e-4,1e5];
    DLS_Input.Kernel = 1;
    DLS_Input.Nnq = 0;
    DLS_Input.Anq = 0;
    DLS_Input.Neq = 0;
    DLS_Input.Aeq = 0;
    DLS_Input.Nneg = 1;
    DLS_Input.Ny0 = 1;
    DLS_Input.iwt = 1;
    DLS_Input.Nalpha = 40;
    DLS_Input.alpha_lims = [0.01,100];
    
    %CRIE run
    s = CRIE(tTemp*1e-6,dataTemp,Ng,0,DLS_Input);
    I_vec(:,j) = s;
end
%% Re-fitting to expand range of fitted g1
I = find(sqrt(ACF_ToFit)/sqrt(beta) > 0.5);
deltaRel = - log10(t(1)*1e-6) + log10(t(I(end))*1e-6);
tFinal = 10^(log10(t(1)*1e-6) + 2*deltaRel);
logSpacing = mean(diff(log10(tTemp*1e-6)));
Nexpand = floor((log10(tFinal) - log10(tTemp(end)*1e-6))/logSpacing);
tExpand = logspace(log10(tTemp(end)*1e-6),log10(tFinal),Nexpand)';
tExpand = tExpand(2:end);
tFit = [tTemp*1e-6;tExpand];

%Kernel generation
gFit_min = log10(DLS_Input.g_lims(1));
gFit_max = log10(DLS_Input.g_lims(2));
gFit = linspace(gFit_min,gFit_max,Ng);
cFit = [4*ones(1,(Ng-3)/2);2*ones(1,(Ng-3)/2)]; %Simpson rule for coeff
cFit = cFit(:);
hFit = (gFit(end) - gFit(1))/(Ng-1);
cFit = (hFit/3)*[1;cFit;4;2];
[gFitM,~] = meshgrid(gFit,1e-6*tTemp);
[cFitM,tFitM] = meshgrid(cFit,1e-6*tTemp);
A = cFitM.*(log(10)*10.^(gFitM)).*exp(-tFitM.*10.^(gFitM));

ACF_vec = zeros(length(tFit),size(data,2));
for i = 1:size(data,2)
    ACF_vec(1:length(tTemp),i) = A*I_vec(:,i);
    Afit = [1,-tTemp(end-1)*1e-6;1,-tTemp(end)*1e-6];
    Bfit = [log(ACF_vec(length(tTemp)-1,i));log(ACF_vec(length(tTemp),i))];
    xfit = Afit\Bfit;
    ACF_vec(length(tTemp)+1:end,i) = exp(xfit(1)).*exp(-xfit(2)*tExpand);
    %Correct 0-intercept by fitting the first 3 points
    if ACF_vec(1,i) > 1
        tCrit = tFit(1:3,:);
        dataTempLog = log(ACF_vec(1:3,i));
        B = [ones(length(tCrit),1),tCrit,tCrit.^2];
        Coeff = B\dataTempLog;
        beta = exp(Coeff(1));
    else
        beta = 1;
    end
    ACF_vec(:,i) = ACF_vec(:,i)/beta;
end
ACF = mean(ACF_vec,2);

%% MSD calculation
tau = tFit;
MSD_vec = (6/q^2)*(-log(ACF_vec));
MSD = mean(MSD_vec,2);

%Plotting
close all
subplot(1,2,1)
semilogx(tTemp,dataTemp,tau*1e6,ACF,'*')
xlabel('t [us]')
ylabel('g1(t)')
title('ACF')
subplot(1,2,2)
loglog(tau*1e6,MSD*1e6,'o',tau*1e6,6*(kB*T/(6*pi*eta*R))*tau*1e6,'k--')
xlabel('Time [us]')
ylabel('MSD [um^2]')
legend('Data','Pure diffusive')
title('MSD')
