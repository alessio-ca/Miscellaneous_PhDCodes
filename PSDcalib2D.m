% Calibration using power spectral density method
%
% Calibrates an optical tweezers using the power spectral density method.

% Clear everything
clear all;
close all;
clc;

%%

% Data matrix
rx=test_dx_x;
ry=test_dx_y;

% Normalize signals to the mean 
rx = rx - repmat(mean(rx,1),size(rx,1),1);
ry = ry - repmat(mean(ry,1),size(ry,1),1);

% Assign 2 calib series: 2D plane (x & y)
Vx=rx;
Vy=ry;


% Coefficients
R=1.47e-6;
eta=1e-3;
T=300;
dt=1/150;
gamma=6*pi*eta*R;
kB=1.38e-23; %Boltzmann constant
kBT=kB*T;
D=kBT/gamma;

% Power spectrum
[PSD_X,F_X,PSD_ERR_X] = psd_routine(dt,Vx, 'fmin',1e-2,'fmax',30,'RunAvg','lin','NumBins', 5000);
[PSD_Y,F_Y,PSD_ERR_Y] = psd_routine(dt,Vy, 'fmin',1e-2,'fmax',30,'RunAvg','lin','NumBins', 5000);


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

D_fit_x = (1/bx) * 2 * (pi.^2);
D_fit_y = (1/by) * 2 * (pi.^2);

Xmin_x = sx00 - (sx01*sx01*sx22 + sx11*sx11*sx02 - 2*sx01*sx11*sx12)/(sx02*sx22 - sx12*sx12);
Xmin_y = sy00 - (sy01*sy01*sy22 + sy11*sy11*sy02 - 2*sy01*sy11*sy12)/(sy02*sy22 - sy12*sy12);

S_fit_x = sqrt(D/D_fit_x);
S_fit_y = sqrt(D/D_fit_y);

psd_fit_x = 1 ./ (ax + bx .* F_X.^2);
psd_fit_y = 1 ./ (ay + by .* F_Y.^2);

% PSD Errors
Times_x = dt*size(Vx,1);
Times_y = dt*size(Vy,1);

xx   = min(F_X) / fcx;
yx   = max(F_X) / fcx;
xy   = min(F_Y) / fcy;
yy   = max(F_Y) / fcy;

sx   = sqrt(pi) * ( (2*yx) / (1 + yx^2) - (2*xx) / (1 + xx^2) + 2 * atan((yx - xx) / (1 + xx*yx)) - ...
    (4/(yx - xx)) * (atan((yx - xx) / (1 + xx*yx)))^2) ^ (-1/2);
sy   = sqrt(pi) * ( (2*yy) / (1 + yy^2) - (2*xy) / (1 + xy^2) + 2 * atan((yy - xy) / (1 + xy*yy)) - ...
    (4/(yy - xy)) * (atan((yy - xy) / (1 + xy*yy)))^2) ^ (-1/2);

fc_err_x = sx * fcx / sqrt(pi * fcx * Times_x);
k_err_x = abs(fc_err_x/fcx)*kx;
fc_err_y = sy * fcy / sqrt(pi * fcy * Times_y);
k_err_y = abs(fc_err_y/fcy)*ky;

gx   = sqrt( ((2*yx)/(1 + yx^2)-(2*xx)/(1 + xx^2) + 2*atan((yx - xx) / (1 + xx*yx)) )/((1 + pi/2)*(yx - xx)) );
gy   = sqrt( ((2*yy)/(1 + yy^2)-(2*xy)/(1 + xy^2) + 2*atan((yy - xy) / (1 + xy*yy)) )/((1 + pi/2)*(yy - xy)) );


D_fit_err_x  = D_fit_x * sqrt( (1 + pi/2) / (pi * fcx * Times_x) )*gx*sx;
D_fit_err_y  = D_fit_y * sqrt( (1 + pi/2) / (pi * fcy * Times_y) )*gy*sy;
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
    'D_fit_x : ' num2str(D_fit_x) ' +/- ' num2str(D_fit_err_x) ' V^2/s\n' ...
    'D_fit_y : ' num2str(D_fit_y) ' +/- ' num2str(D_fit_err_y) ' V^2/s\n' ...
    'S_fit_x : ' num2str(S_fit_x) ' +/- ' num2str(S_fit_err_x) ' V/m\n' ...
    'S_fit_y : ' num2str(S_fit_y) ' +/- ' num2str(S_fit_err_y) ' V/m\n' ...
    '\n'
    ];
fprintf(txt)
%%
errorbar(F_X,PSD_X*1e+18,PSD_ERR_X*1e+18,'.')
hold on
errorbar(F_Y,PSD_Y*1e+18,PSD_ERR_Y*1e+18,'.')
legend('Focal plane (X)','Focal plane (Y)')

plot(F_X,psd_fit_x*1e+18,'k')
plot(F_Y,psd_fit_y*1e+18,'k')

hold off
xlabel('f [Hz]')
ylabel('PSD [nm^2/Hz]')
set(gca,'XScale','Log','YScale','Log')