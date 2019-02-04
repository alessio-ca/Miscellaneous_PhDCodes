% Calibration using mean square displacement
%
% Calibrates an optical tweezers using the mean square displacement

% Clear everything
clear all;
close all;
clc;

%%

% Data matrix
data = dlmread('L:\OT_ParticleWall\z=1d_G=0_0025.csv');

rx = data(:,1:2:end-1);
ry = data(:,2:2:end);

% Normalize signals to the mean 
rx = rx - repmat(mean(rx,1),size(rx,1),1);
ry = ry - repmat(mean(ry,1),size(ry,1),1);

% Assign 2 calib series: 2D plane (x & y)
Vx=1e-6*rx/10; %pixel-to-um conversion
Vy=1e-6*ry/10;

% Coefficients
R=0.75e-6;
eta=8.94*1e-4;
T=298;
dt=1/1127;
gamma=6*pi*eta*R;
kB=1.38e-23; %Boltzmann constant
kBT=kB*T;
D=kBT/gamma;
%% MSD inspection: no fitting
[MSD_X,T_X,~,MSD_X_VEC] = msd_routine(dt,Vx,'TauMax',5);
[MSD_Y,T_Y,~,MSD_Y_VEC] = msd_routine(dt,Vy,'TauMax',5);

plot(T_X,MSD_X*1e+12,'b')
hold on
plot(T_Y,MSD_Y*1e+12,'r')
hold off
legend('X','Y')
xlabel('t [s]')
ylabel('MSD [um^2/s]')
set(gca,'XScale','Log','YScale','Log')

%% MSD Fitting
% Fitting
g = @(a, b, x) a + log10(1 - exp(-x/b)); %Fitting Model

%MSD tail fitting
MSD_X_END = mean(MSD_X(1000:end)); 
MSD_Y_END = mean(MSD_Y(1000:end));

%Fit: start points estimation
kx_0 = 2*kB*T/MSD_X_END; 
ky_0 = 2*kB*T/MSD_Y_END;
tauc_x_0 = 6*pi*eta*R/kx_0;
tauc_y_0 = 6*pi*eta*R/ky_0;

W = length(T_X) - (10:1:length(MSD_X)); %Fitting weights

%Fitting
fitresult_x = fit(T_X(10:end),log10(MSD_X(10:end)),fittype(g),'StartPoint',[log10(2*kB*T/kx_0) tauc_x_0],'Weight',W.^4);
fitresult_y = fit(T_Y(10:end),log10(MSD_Y(10:end)),fittype(g),'StartPoint',[log10(2*kB*T/ky_0) tauc_y_0],'Weight',W.^4);

%Fit coefficients
ax = fitresult_x.a; 
ay = fitresult_y.a;
bx = fitresult_x.b;
by = fitresult_y.b;

%Fit result
tauc_x = bx;
tauc_y = by;
kx = gamma/tauc_x;
ky = gamma/tauc_y;

D_fit_x = 10^(ax)/(2*tauc_x);
D_fit_y = 10^(ay)/(2*tauc_y);


plot(T_X,MSD_X*1e+12,'b')
hold on
plot(T_X,10.^g(ax,bx,T_X)*1e+12,'k--')
plot(T_Y,MSD_Y*1e+12,'r')
plot(T_Y,10.^g(ay,by,T_Y)*1e+12,'k-.')
hold off
legend('X','X-Fit','Y','Y-Fit')
xlabel('t [s]')
ylabel('MSD [um^2/s]')
set(gca,'XScale','Log','YScale','Log')
          
%% MSD Fitting Errors
tauc_x_vec = zeros(size(rx,2),1);
kx_vec = zeros(size(rx,2),1);
Dx_vec = zeros(size(rx,2),1);

%Build fit coefficient distribution from data series
for n = 1:size(rx,2)
    fitresult = fit(T_X(2:end),log10(MSD_X_VEC(2:end,n)),fittype(g),'StartPoint',[log10(D_fit_x*tauc_x) tauc_x]);
    a = fitresult.a;
    b = fitresult.b;
    tauc_x_vec(n) = b;
    kx_vec(n) = gamma/tauc_x_vec(n);
    Dx_vec(n) = 10^a/(2*tauc_x_vec(n));
end
%Estimate std of coefficient distro
tauc_x_err = std(tauc_x_vec);
kx_err = std(kx_vec);
Dx_err = std(Dx_vec);

tauc_y_vec = zeros(size(ry,2),1);
ky_vec = zeros(size(ry,2),1);
Dy_vec = zeros(size(ry,2),1);
%Build fit coefficient distribution from data series
for n = 1:size(ry,2)
    fitresult = fit(T_Y(2:end),log10(MSD_Y_VEC(2:end,n)),fittype(g),'StartPoint',[log10(D_fit_y*tauc_y) tauc_y]);
    b = fitresult.b;
    tauc_y_vec(n) = b;
    ky_vec(n) = gamma/tauc_y_vec(n);
    Dy_vec(n) = 10^a/(2*tauc_y_vec(n));
end
%Estimate std of coefficient distro
tauc_y_err = std(tauc_y_vec);
ky_err = std(ky_vec);
Dy_err = std(Dy_vec);

% Print results
txt = ['Mean squared displacement analysis: \n' ...
    int2str(size(Vx,2)) ' signals with ' int2str(size(Vx,1)) ' samples each\n' ...
    int2str(size(Vy,2)) ' signals with ' int2str(size(Vy,1)) ' samples each\n' ...
    '\n' ...
    'tauc_x : ' num2str(tauc_x) ' +/- ' num2str(tauc_x_err) ' s\n' ...
    'tauc_y : ' num2str(tauc_y) ' +/- ' num2str(tauc_y_err) ' s\n' ...
    'kx : ' num2str(kx*1e+6) ' +/- ' num2str(kx_err*1e+6) ' fN/nm\n' ...
    'ky : ' num2str(ky*1e+6) ' +/- ' num2str(ky_err*1e+6) ' fN/nm\n' ...
    'D_fit_x : ' num2str(D_fit_x*1e+12) ' +/- ' num2str(Dx_err*1e+12) ' um^2/s\n' ...
    'D_fit_y : ' num2str(D_fit_y*1e+12) ' +/- ' num2str(Dy_err*1e+12) ' um^2/s\n' ...
    '\n' 
    ];
fprintf(txt)
