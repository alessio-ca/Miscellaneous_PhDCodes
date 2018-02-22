%Script to plot the data from DLS_Analysis 
%Input data are variables with names
% - tau
% - ACF
% - MSD

%You need to be in the folder where the "mat" files are
%Assumes decreasing T

clear variables
close all

NM = 5; %Number of measurements for each temperature
NT = 9; %Number of temperatures
Tin = 50;
Tstep = 3;
picks = [1,4,1,2,4,2,4,1,2];


%% Constants 
n = 1.33;
lambda = 633*10^(-9);
theta = 173*2*pi/360;
q = 4*pi*n*sin(theta/2)/lambda;


%% Preliminary: renormalization & colormaps
f0 = 0;
f1 = 1;
for i=1:NT
    j = picks(i);
    load(['T_',num2str(Tin - (i-1)*Tstep),'_',num2str(j)]);
    f0 = max([ACF(1,1),f0]);
    f1 = min([ACF(1,1),f1]);
end
initialColorOrder = get(gca,'ColorOrder'); 
newDefaultColors = flip(cool(NT));
f = linspace(f1,f0,NT);

%% ACF plot
figure(1)
set(gca, 'ColorOrder', newDefaultColors, 'NextPlot', 'replacechildren');
hold on

for i=1:NT
    j = picks(i);
    load(['T_',num2str(Tin - (i-1)*Tstep),'_',num2str(j)]);
    semilogx(tau,f(i)*ACF./ACF(1,1));
end
hold off
h=gca;
h.XScale='log';
xlabel('Time (us)')
ylabel('g_1 (\tau)')
legend(cellstr(int2str(linspace(Tin,Tin - (NT-1)*Tstep,NT)')))

%% MSD plot
figure(2)
set(gca, 'ColorOrder', newDefaultColors, 'NextPlot', 'replacechildren');
hold on
for i=1:NT
    j = picks(i);
    load(['T_',num2str(Tin - (i-1)*Tstep),'_',num2str(j)]);
    f0 = ACF(1,1);
    loglog(tau,(6/q^2)*(-log(f(i)*ACF./ACF(1,1))));
end
hold off
h=gca;
h.XScale='log';
h.YScale='log';
xlabel('Time (us)')
ylabel('MSD (m^2)')
legend(cellstr(int2str(linspace(Tin,Tin - (NT-1)*Tstep,NT)')),'Location','southeast')



