%Script to plot the ACFs and MSDs data from DLS_Analysis 
%Input data are variables with names
% - tau
% - ACF
% - MSD

%Assumes decreasing temperatures with Trend = -1
%Assumes increasing temperatures with Trend = +1

clear variables
close all

NM = 3; %Number of measurements for each temperature
NT = 2; %Number of temperatures
Tin = 20;
Tstep = 2;
Id = 'SNR';

Trend = +1;

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

%% Preliminary: renormalization & colormaps
fP = ones(NT,1);
fM = ones(NT,1);
for i=1:NT
    load([PathName,Id,'_T_',num2str(Tin + Trend*(i-1)*Tstep)]);
end
initialColorOrder = get(gca,'ColorOrder');
newDefaultColors = flip(cool(NT));
%% ACF plot
figure(1)
subplot(1,2,1)
set(gca, 'ColorOrder', newDefaultColors, 'NextPlot', 'replacechildren');
hold on

for i=1:NT
    load([PathName,Id,'_T_',num2str(Tin + Trend*(i-1)*Tstep)]);
    semilogx(tau,fM(i)*ACF.^fP(i));
end
hold off
h=gca;
h.XScale='log';
xlabel('Time [s]')
ylabel('g_1 (\tau)')
legend(cellstr(int2str(linspace(Tin,Tin + Trend*(NT-1)*Tstep,NT)')))

%% MSD plot
figure(1)
subplot(1,2,2)
set(gca, 'ColorOrder', newDefaultColors, 'NextPlot', 'replacechildren');
hold on
for i=1:NT
    load([PathName,Id,'_T_',num2str(Tin + Trend*(i-1)*Tstep)]);
    loglog(tau,1e+12*(6/q^2)*(-log(fM(i).*ACF.^fP(i))),'*');
end
loglog(tau,1e+12*6*D(Tin)*tau,'k--')
hold off
h=gca;
h.XScale='log';
h.YScale='log';
xlabel('Time [s]')
ylabel('MSD [um^2]')
legend(cellstr(int2str(linspace(Tin,Tin + Trend*(NT-1)*Tstep,NT)')),'Location','southeast')
