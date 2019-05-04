%% DLS_Gmoduli

%Script to extract G* from the data from DLS_Analysis 
%Input data are variables with names
% - tau
% - ACF
% - MSD

%Assumes decreasing temperatures with Trend = -1
%Assumes increasing temperatures with Trend = +1

%Recalculates MSDs according to a renormalization

clear variables
close all

NT = 11; %Number of temperatures
Tin = 60;
Tstep = 2;
Id = 'SNR';

Trend = -1;

PathName = uigetdir;
PathName = [PathName,'/'];
FileList = dir([PathName,Id,'_T*']);
if size(FileList,1) ~= NT
    error('Wrong datafile length!');
end

%% Constants 
n = 1.33;
lambda = 633*10^(-9);
theta = 173*2*pi/360;
q = 4*pi*n*sin(theta/2)/lambda;
kB = 1.38*10^(-23);
A = 2.414e-5;
B = 247.8;
C = 140;
eta = @(x) A * 10.^(B/((x+273.15) - C));
D = @(x) kB*(x + 273.15)/(3*pi*230e-9*eta(x));

%% Preliminary: renormalization
fP = ones(NT,1);
fM = ones(NT,1);

for i=1:NT
    load([PathName,Id,'_T_',num2str(Tin + Trend*(i-1)*Tstep)]);
    %fP(i) = 6*D(Tin + Trend*(i-1)*Tstep)*tau(1)/(sqrt(i));
    %fP(i) = -fP(i)*(q^2/6)/log(ACF(1));
end
%fP(2) = fP(2)*1.2;
%fP(4) = fP(4)*0.85;
%fP = 0.5*fP;
%% G* calculation
for i=1:NT
    load([PathName,Id,'_T_',num2str(Tin + Trend*(i-1)*Tstep)]);
    tau(ACF==0) = [];
    ACF(ACF==0) = []; %Remove zeros
    MSD = (6/q^2)*(-log(fM(i).*ACF.^fP(i)));
    [omega,G]=MSDtoG_Mason(tau*1e-6,MSD,'R',115e-9,'CG',1.01,'T',273.15 + Tin + Trend*(i-1)*Tstep);
    %[omega,G]=MSDtoG_Evans_oversampling(tau,MSD,1/tau(1),'R',115e-9,'CG',1.01,'T',273.15 + Tin - (i-1)*Tstep,'Nf0',5,'Nfinf',5);
    save([PathName,'T_',num2str(Tin + Trend*(i-1)*Tstep),'_G'],'omega','G');
end





