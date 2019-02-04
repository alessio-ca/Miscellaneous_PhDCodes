%Script to plot the data from DLS_Analysis 
%Input data are variables with names
% - tau
% - ACF
% - MSD

%Assumes decreasing temperatures with Trend = -1
%Assumes increasing temperatures with Trend = +1

clear variables
close all

NM = 2; %Number of measurements for each temperature
NT = 5; %Number of temperatures
Tin = 54;
Tstep = 2;

Trend = -1;

[FileName,PathName] = uigetfile;
load([PathName,FileName]);
if size(data,1) ~= NM*NT
    error('Wrong datafile length!');
end

Nid = 1;
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
%% ACF plot
figure(1)
subplot(1,2,1)
hold on
for j = 1:NM
    load([PathName,'T_',num2str(Tin + Trend*(Nid-1)*Tstep),'_',num2str(j)]);
    semilogx(tau,ACF);
end

hold off
h=gca;
h.XScale='log';
xlabel('Time [s]')
ylabel('g_1 (\tau)')
%% MSD plot
figure(1)
subplot(1,2,2)
hold on
for j = 1:NM
    load([PathName,'T_',num2str(Tin + Trend*(Nid-1)*Tstep),'_',num2str(j)]);
    loglog(tau,(6/q^2)*(-log(ACF)),'*');
end

loglog(tau,6*D(Tin)*tau,'k--')
hold off
h=gca;
h.XScale='log';
h.YScale='log';
xlabel('Time [s]')
ylabel('MSD [m^2]')

