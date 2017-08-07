% Calibration using power spectral density method
%
% Calibrates an optical tweezers using the power spectral density method.

% Clear everything
clear all;
close all;
clc;

%%

% Data matrix
%rx=dataTot(:,1:5);
%ry=dataTot(:,6:10);
rx = r(:,1:10);
ry = r(:,11:20);

% Normalize signals to the mean 
rx = rx - repmat(mean(rx,1),size(rx,1),1);
ry = ry - repmat(mean(ry,1),size(ry,1),1);

% Assign 2 calib series: 2D plane (x & y)
Vx=rx;
Vy=ry;


% Coefficients
R=0.75e-6;
eta=1e-3;
T=300;
dt=1/1e4;
gamma=6*pi*eta*R;
kB=1.38e-23; %Boltzmann constant
kBT=kB*T;
D=kBT/gamma;
%% Power spectrum inspection: no fitting
[PSD_X,F_X] = psd_routine(dt,Vx);
[PSD_Y,F_Y] = psd_routine(dt,Vy);

plot(F_X,PSD_X,'b')
hold on
plot(F_Y,PSD_Y,'r')
hold off
legend('X','Y')
xlabel('f [Hz]')
ylabel('PSD [mV^2/Hz]')
set(gca,'XScale','Log','YScale','Log')

%%
% Power spectrum
[PSD_X,F_X,PSD_ERR_X] = psd_routine(dt,Vx, 'fmin',.1,'fmax',1e3,'RunAvg','lin','NumBins', 500);
[PSD_Y,F_Y,PSD_ERR_Y] = psd_routine(dt,Vy, 'fmin',.1,'fmax',1e3,'RunAvg','lin','NumBins', 500);


% PSD Fitting
for p = 0:1:2
    for q = 0:1:2
        eval(['sx' num2str(p) num2str(q) ' = sum ((F_X.^ (2*' num2str(p) ')) .* (PSD_X.^ (' num2str(q) ')));']);
        eval(['sy' num2str(p) num2str(q) ' = sum ((F_Y.^ (2*' num2str(p) ')) .* (PSD_Y.^ (' num2str(q) ')));']);
    end
end

ax = (sx01 * sx22 - sx11 * sx12) / (sx02 * sx22 - sx12.^2);
bx = (sx11 * sx02 - sx01 * sx12) / (sx02 * sx22 - sx12.^2);
ay = (sy01 * sy22 - sy11 * sy12) / (sy02 * sy22 - sy12.^2);
by = (sy11 * sy02 - sy01 * sy12) / (sy02 * sy22 - sy12.^2);

fcx = sqrt(ax/bx);
fcy = sqrt(ay/by);

kx = 2*pi*gamma*fcx;
ky=  2*pi*gamma*fcy;

Times_x = dt*size(Vx,1);
Times_y = dt*size(Vy,1);

D_fit_x = (1/bx) * 2 * (pi.^2);
D_fit_y = (1/by) * 2 * (pi.^2);

Chimin_x = sx00 - (sx01*sx01*sx22 + sx11*sx11*sx02 - 2*sx01*sx11*sx12)/(sx02*sx22 - sx12*sx12);
Chimin_y = sy00 - (sy01*sy01*sy22 + sy11*sy11*sy02 - 2*sy01*sy11*sy12)/(sy02*sy22 - sy12*sy12);

S_fit_x = sqrt(D/(D_fit_x*1e-6));
S_fit_y = sqrt(D/(D_fit_y*1e-6));

psd_fit_x = 1 ./ (ax + bx .* F_X.^2);
psd_fit_y = 1 ./ (ay + by .* F_Y.^2);

% PSD Errors

minx   = min(F_X) / fcx;
maxx   = max(F_X) / fcx;
miny   = min(F_Y) / fcy;
maxy   = max(F_Y) / fcy;

u = @(x1,x2) 2 * x2 / (1 + x2^2) - 2 * x1 / (1 + x1^2) + 2 * atan((x2 - x1)/(1 + x1 * x2));
v = @(x1,x2) 4/(x2 - x1) * atan((x2 - x1)/(1 + x1 * x2))^2;
sf = @(x1,x2) sqrt(pi / (u(x1,x2) - v(x1,x2)));
sd = @(x1,x2) sqrt(u(x1,x2) / ((1 + pi/2) * (x2 - x1)));

sfx = sf(minx,maxx);
sfy = sf(miny,maxy);
sdx = sd(minx,maxx);
sdy = sd(miny,maxy);

fc_err_x = sfx * fcx / sqrt(pi * fcx * Times_x);
k_err_x = abs(fc_err_x/fcx)*kx;
fc_err_y = sfy * fcy / sqrt(pi * fcy * Times_y);
k_err_y = abs(fc_err_y/fcy)*ky;

D_fit_err_x  = D_fit_x * sqrt( (1 + pi/2) / (pi * fcx * Times_x) )*sdx;
D_fit_err_y  = D_fit_y * sqrt( (1 + pi/2) / (pi * fcy * Times_y) )*sdy;
S_fit_err_x = 0.5 * S_fit_x * D_fit_err_x / D_fit_x;
S_fit_err_y = 0.5 * S_fit_y * D_fit_err_y / D_fit_y;



%Print results
txt = ['Power spectrum analysis: \n' ...
    int2str(size(Vx,2)) ' signals with ' int2str(size(Vx,1)) ' samples each\n' ...
    int2str(size(Vy,2)) ' signals with ' int2str(size(Vy,1)) ' samples each\n' ...
    '\n' ...
    'fcx : ' num2str(fcx) ' +/- ' num2str(fc_err_x) ' Hz\n' ...
    'fcy : ' num2str(fcy) ' +/- ' num2str(fc_err_y) ' Hz\n' ...
    'kx : ' num2str(kx*1e+6) ' +/- ' num2str(k_err_x*1e+6) ' fN/nm\n' ...
    'ky : ' num2str(ky*1e+6) ' +/- ' num2str(k_err_y*1e+6) ' fN/nm\n' ...
    'D_fit_x : ' num2str(D_fit_x) ' +/- ' num2str(D_fit_err_x) ' mV^2/s\n' ...
    'D_fit_y : ' num2str(D_fit_y) ' +/- ' num2str(D_fit_err_y) ' mV^2/s\n' ...
    'S_fit_x : ' num2str(S_fit_x) ' +/- ' num2str(S_fit_err_x) ' m/V\n' ...
    'S_fit_y : ' num2str(S_fit_y) ' +/- ' num2str(S_fit_err_y) ' m/V\n' ...
    '\n'
    ];
fprintf(txt)
%%
%Conversion from mV to nm
factor_x = S_fit_x *(1e+9/1e+3);
factor_y = S_fit_y *(1e+9/1e+3);

errorbar(F_X,factor_x^2 * PSD_X, factor_x^2 * PSD_ERR_X,'.')
hold on
errorbar(F_Y,factor_y^2 * PSD_Y, factor_y^2 * PSD_ERR_Y,'.')
legend('Focal plane (X)','Focal plane (Y)')

plot(F_X,psd_fit_x*factor_x^2,'k')
plot(F_Y,psd_fit_y*factor_y^2,'k')

hold off
xlabel('f [Hz]')
ylabel('PSD [nm^2/Hz]')
set(gca,'XScale','Log','YScale','Log')

%% Optional: export data (change names!)
F_X_syn = F_X;
F_Y_syn = F_Y;
PSD_X_syn = factor_x^2 * PSD_X;
PSD_Y_syn = factor_y^2 * PSD_Y;
PSD_X_fit_syn = psd_fit_x*factor_x^2;
PSD_Y_fit_syn = psd_fit_y*factor_y^2;
PSD_ERR_X_syn = factor_x^2 * PSD_ERR_X;
PSD_ERR_Y_syn = factor_y^2 * PSD_ERR_Y;
save('fit_PSD_10kHz_syn.mat',...
    'F_X_syn', 'F_Y_syn', ...
    'PSD_X_syn','PSD_Y_syn', ...
    'PSD_X_fit_syn','PSD_Y_fit_syn',...
    'PSD_ERR_X_syn','PSD_ERR_Y_syn')

