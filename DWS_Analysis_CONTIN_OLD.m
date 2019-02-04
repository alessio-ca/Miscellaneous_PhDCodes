function [tau,MSD,ACF,MSD_vec,ACF_vec]=DWS_Analysis_CONTIN_OLD(t,data,varargin)
% DLS_Analysis Estimation of MSD and G* from DWS autocorrelation
% data
%
% [tau,MSD,ACF,MSD_vec,ACF_vec]=DWS_Analysis(t,data) calculates the MSD and complex
%   modulus G* from the autocorrelation function ICF obtained in a DWS
%   experiment, using the inversion relation and Evans scheme for inverting
%   the MSD into the complex modulus G. The measurement setup is assumed to 
%   be a DWS Rheolab (LS instruments)
%  
%  DATA can be a matrix containing several signals (one per column, all the same length).
%  MSD_vec is the set of all the calculated MSDs.
%  ACF_vec is the set of all the calculated ACFs.

% [tau,MSD,ACF,MSD_vec,ACF_vec]=DWS_Analysis(t,data,'PropertyName',PropertyValue)
%  permits to set the value of PropertyName to PropertyValue.
%  Admissible Properties are:
%       T       -   Temperature (default = 298 K)
%       eta     -   Solvent viscosity (default = 0.8872 cP)
%       n       -   Solvent refractive index (default = 1.33)
%       R       -   Bead radius (default = 115 nm)
%       tail    -   Fit times for tail characterization (default =
%                   [0,0] i.e. no tail correction)
%       fitmeth -   Fit method used (default = "CON")
%                   Allowed methods: CON (CONTIN), RAT (rational fit), SPL (spline)


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
tail = [0,0];
for i = 1:2:length(varargin)
    if strcmpi(varargin{i},'tail')
        tail = varargin{i+1};
    end
end
fitmeth = 'CON';
for i = 1:2:length(varargin)
    if strcmpi(varargin{i},'fit')
        fitmeth = varargin{i+1};
    end
end
%% Constants
kB = 1.38*10^(-23);
lambda = 633*10^(-9);
theta = 173*2*pi/360;
q = 4*pi*n*sin(theta/2)/lambda;

%% g1(tau) calc. and fit to g1(tau) (for smoothness)

alpha = 1e-2;
I_vec = ones(100,size(data,2));
fit_vec = zeros(7,size(data,2));

for j = 1:size(data,2)
    
    %Delete first 2 points (usually noise) & format data
    Time_ToFit = t(3:end);
    ACF_ToFit = data(3:end,j);
    if all(tail)
        tailFit = ACF_ToFit(Time_ToFit > tail(1) & Time_ToFit < tail(2));
        tailFit = mean(tailFit);
        ACF_ToFit = ACF_ToFit - tailFit;
    end
    CritPts = find(ACF_ToFit< 0.3*ACF_ToFit(1));
    tTemp = Time_ToFit(1:CritPts(1));
    
    %g1(tau) calculation with  optional tail treatment (least-Squares)
    if all(tail)
        dataTempLog = log(ACF_ToFit(1:CritPts));
        A = [ones(length(tTemp),1),tTemp,tTemp.^2];
        Coeff = A\dataTempLog;
        beta = exp(Coeff(1));
    else
        beta = 1;
    end
    dataTemp = sqrt(ACF_ToFit(1:CritPts))/sqrt(beta);
    
    switch fitmeth
        case 'CON'
            %CONTIN loop run
            g_inf = -6;
            g_sup = -1;
            check = 1;
            while(check>0)
                check=0;
                coarse_s = logspace(g_inf,g_sup,10)';
                coarse_g = ones(size(coarse_s));
                [coarse_g,~,~] = CONTIN_Rilt(tTemp,dataTemp,coarse_s,coarse_g,alpha,'logarithmic',100,[],[],[],[],[]);
                % s refinement
                g_peak = find(coarse_g>0.01);
                if g_peak(1)<4
                    g_inf = g_inf-1;
                    check=1;
                end
                if g_peak(end)>7
                    g_sup = g_sup+1;
                    g_inf = g_inf+1;
                    check=1;
                end
            end
            for i = 2:10
                disp(['Iteration: ',num2str(i),'/10'])
                s = logspace(g_inf,g_sup,i*10)';
                g0 = interp1(coarse_s,coarse_g,s,'linear');
                [g,~,~] = CONTIN_Rilt(tTemp,dataTemp,s,g0,alpha,'logarithmic',100,[],[],[],[],[]);
                coarse_g = g;
                coarse_s = s;
                disp(' ')
            end
            I_vec(:,j) = g;
            I_vec = I_vec./repmat(sum(I_vec,1),[length(I_vec),1]);
            
        case 'RAT'
            %Rational fit loop run
            U = tTemp;
            Z = dataTemp.*tTemp;
            V = dataTemp.*tTemp.^2;
            F = dataTemp.*tTemp.^3;
            K = dataTemp;
            A = [U,ones(length(U),1),-V,-Z,-K];
            Coeff = A\F;
            fo = fitoptions('rat33');
            fo.StartPoint = [0,0,Coeff'];
            my_fit = fit(tTemp,dataTemp,'rat33',fo);
            fit_vec(:,j) = coeffvalues(my_fit);
            
        case 'SPL'
            %Spline fitting
            data_spline=csape(tTemp',dataTemp');
            
        otherwise
            error('Error: Unrecognized fit method.');
    end
end
%%
switch fitmeth
    case 'CON'
        tFit = logspace(floor(log10(tTemp(1))),1+log10(tTemp(end)))'; %Extrapolated
        ACF_vec = zeros(length(tFit),size(data,2));
        for i = 1:length(tFit)
            ACF_vec(i,:) = sum(I_vec.*exp(-tFit(i)./s),1);
        end
        
    case 'RAT'
        tFit = logspace(floor(log10(tTemp(1))),1+log10(tTemp(end)))'; %Extrapolated
        ACF_vec = zeros(length(tFit),size(data,2));
        fit_33 = @(x,i) (fit_vec(1,i).*x.^3 + fit_vec(2,i).*x.^2 + fit_vec(3,i).*x + fit_vec(4,i))./(x.^3 + fit_vec(5,i).*x.^2 + fit_vec(6,i).*x + fit_vec(7,i));
        for i = 1:size(data,2)
            ACF_vec(:,i) = fit_33(tFit,i)./(fit_vec(4,i)./fit_vec(7,i));
            if ACF_vec(end,i) < 0
                ACF_vec(:,i) = (ACF_vec(:,i) - 2*min(ACF_vec(:,i)))./(1 -  2*min(ACF_vec(:,i)));
            end
        end
        
    case 'SPL'
        tFit = logspace(log10(tTemp(1)),log10(tTemp(end)))'; %Non-extrapolated
        ACF_vec=fnval(tFit,data_spline)';
        if size(ACF_vec,1)==1
            ACF_vec=ACF_vec';
        end
        ACF_vec = ACF_vec./fnval(0.999*tFit(1),data_spline);
        
    otherwise
        error('Error: Unrecognized fit method.');
end
ACF = mean(ACF_vec,2);

close all
subplot(1,2,1)
semilogy(tTemp,dataTemp,'o',tFit(tFit<=tTemp(end)),ACF(tFit<=tTemp(end)))
xlabel('Time [s]')
ylabel('g1(t)')
title('g1 (raw and smooth fit)')
legend('Raw','Fit')

%% MSD calculation
tau = tFit;
MSD_vec = (6/q^2)*(-log(ACF_vec));
MSD = mean(MSD_vec,2);
subplot(1,2,2)
loglog(tau,MSD,'o',tau,6*(kB*T/(6*pi*eta*R))*tau,'k--')
xlabel('Time [s]')
ylabel('MSD(t)')
legend('Data','Pure diffusive')
title('MSD')
