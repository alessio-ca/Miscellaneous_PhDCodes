function [tau,MSD,ACF,MSD_vec,ACF_vec]=DWS_Analysis(t,data,varargin)
% DWS_Analysis Estimation of MSD and G* from DWS autocorrelation
% data
%
% [MSD,G,MSD_vec,G_vec]=DWS_Analysis(t,data) calculates the MSD and complex
%   modulus G* from the autocorrelation function ICF obtained in a DWS
%   experiment, using the inversion relation and Evans scheme for inverting
%   the MSD into the complex modulus G. The measurement setup is assumed to 
%   be a DWS Rheolab (LS instruments)
%  
%  DATA can be a matrix containing several signals (one per column, all the same length).
%  MSD_vec is the set of all the calculated MSDs.
%  G_vec is the set of all the calculated Gs.

% [MSD,G,MSD_vec,G_vec]=DWS_Analysis(t,data,'PropertyName',PropertyValue)
%  permits to set the value of PropertyName to PropertyValue.
%  Admissible Properties are:
%       T       -   Temperature (default = 298 K)
%       eta     -   Solvent viscosity (default = 0.8872 cP)
%       n       -   Solvent refractive index (default = 1.33)
%       tail    -   Fit times for tail characterization (default =
%       [0,0] i.e. no tail correction)
%       fitmeth -   Fit method used (default = "CON")
%       Allowed methods: CON (CONTIN), RAT (rational fit), SPL (spline)


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
lambda = 685*10^(-9); %Laser wavelength
L = 2e-3; %Cuvette thickness
l_star = 401.87e-6; %Transport mean free path
k0 = 2*pi*n/lambda;
R = 115e-9; %Bead radius

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
    CritPts = find(ACF_ToFit< 0.1*ACF_ToFit(1));
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
            coarse_s = logspace(-6,-1,10)';
            coarse_g = ones(size(coarse_s));
            [coarse_g,~,~] = CONTIN_Rilt(tTemp,dataTemp,coarse_s,coarse_g,alpha,'logarithmic',100,[],[],[],[],[]);
            for i = 2:10
                disp(['Iteration: ',num2str(i)])
                s = logspace(-6,-1,i*10)';
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
close all
%%
switch fitmeth
    case 'CON'
        tFit = logspace(-6,-2)'; %Extrapolated
        ACF_vec = zeros(length(tFit),size(data,2));
        for i = 1:length(tFit)
            ACF_vec(i,:) = sum(I_vec.*exp(-tFit(i)./s),1);
        end
        
    case 'RAT'
        tFit = logspace(-6,-2)'; %Extrapolated
        ACF_vec = zeros(length(tFit),size(data,2));
        fit_33 = @(x,i) (fit_vec(1,i).*x.^3 + fit_vec(2,i).*x.^2 + fit_vec(3,i).*x + fit_vec(4,i))./(x.^3 + fit_vec(5,i).*x.^2 + fit_vec(6,i).*x + fit_vec(7,i));
        for i = 1:size(data,2)
            ACF_vec(:,i) = fit_33(tFit,i)./(fit_vec(4,i)./fit_vec(7,i));
            if ACF_vec(end,i) < 0
                ACF_vec(:,i) = (ACF_vec(:,i) - 2*min(ACF_vec(:,i)))./(1 -  2*min(ACF_vec(:,i)));
            end
        end
        
    case 'SPL'
        tFit = logspace(-6,log10(tTemp(end)))'; %Non-extrapolated
        ACF_vec = zeros(length(tFit),size(data,2));
        ACF_vec=fnval(tFit,data_spline)';
        if size(ACF_vec,1)==1
            ACF_vec=ACF_vec';
        end
        
    otherwise
        error('Error: Unrecognized fit method.');
end
ACF = mean(ACF_vec,2);
semilogy(tTemp,dataTemp,'o',tFit(tFit<=tTemp(end)),ACF(tFit<=tTemp(end)))
%% MSD calculation
g1_an = @(x) (L/l_star + 4/3)/(5/3) * (sinh(x) + (2/3)*x*cosh(x)) / ...
    ((1 + (4/9)*x^2)*sinh((L/l_star)*x) + (4/3)*x*cosh((L/l_star)*x));
x0 = sqrt(-3*log(ACF_vec))./(L/l_star);
xAcc = zeros(size(x0));

for i = 1:length(x0)
    for j = 1:size(x0,2)
        fun = @(x) g1_an(x) - ACF_vec(i,j);
        xAcc(i,j) = fzero(fun,x0(i,j));
    end
end

tau = tFit;
%MSD_vec = 0.33 * (xAcc + x0).^2/k0^2;
MSD_vec = xAcc.^2/k0^2;
MSD = mean(MSD_vec,2);
%MSDFit_Spl = csape(tFit',MSD_vec');
% %% Microrheology
% t_Micro = (1e-6:1e-6:1)';
% MSD_Micro = fnval(t_Micro,MSDFit_Spl)';
% if size(MSD_Micro,1)==1
%     MSD_Micro=MSD_Micro';
% end
% Jfactor = 1 / (kB*T/(pi*R));
% [omega,G,~,G_vec]=MSDtoG_Evans_oversampling(t_Micro,MSD_Micro,1e+6,'Jfactor',Jfactor);



