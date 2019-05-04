%Script to plot the Gmoduli data from DLS_Analysis 
%Input data are variables with names
% - tau
% - ACF
% - MSD

%Assumes decreasing temperatures with Trend = -1
%Assumes increasing temperatures with Trend = +1

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

%% Preliminary: colormaps
fP = ones(NT,1);
fM = ones(NT,1);
initialColorOrder = get(gca,'ColorOrder');
newDefaultColors = flip(cool(NT));

%% G' plot
figure(1)
subplot(1,2,1)
set(gca, 'ColorOrder', newDefaultColors, 'NextPlot', 'replacechildren');
hold on
for i=1:NT
    load([PathName,'T_',num2str(Tin + Trend*(i-1)*Tstep),'_G.mat']);
    loglog(omega,real(G))
end
hold off
h=gca;
h.XScale='log';
h.YScale='log';
xlabel('\omega (rad/s)')
ylabel('G'' (\omega) [Pa]')
legend(cellstr(int2str(linspace(Tin,Tin + Trend*(NT-1)*Tstep,NT)')),'Location','northwest')

%% G'' plot
figure(1)
subplot(1,2,2)
set(gca, 'ColorOrder', newDefaultColors, 'NextPlot', 'replacechildren');
hold on

for i=1:NT
    load([PathName,'T_',num2str(Tin + Trend*(i-1)*Tstep),'_G.mat']);
    loglog(omega,imag(G),'*')
end
loglog(omega,eta(Tin)*omega,'k--')
hold off
h=gca;
h.XScale='log';
h.YScale='log';
xlabel('\omega (rad/s)')
ylabel('G'''' (\omega) [Pa]')
legend(cellstr(int2str(linspace(Tin,Tin + Trend*(NT-1)*Tstep,NT)')),'Location','southeast')