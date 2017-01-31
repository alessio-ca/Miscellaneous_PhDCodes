% Calibration using power spectral density method
%
% Calibrates an optical tweezers using the power spectral density method.

clear all;
close all;
clc;

% Load data file
load('data.mat')

% Normalize signals to the mean 
rx = rx - repmat(mean(rx),size(rx,1),1);
ry = ry - repmat(mean(ry),size(ry,1),1);
rz = rz - repmat(mean(rz),size(rz,1),1);

% Assign 2 calib series: 2D plane (x & y) and z-direction
Vx=[rx ry];
Vz=rz;


% Coefficients
gamma=6*pi*eta*R;
kB=1.38e-23; %Boltzmann constant
kBT=kB*T;
D=kBT/gamma;

% Power spectrum
[PSD_X,F_X,PSD_ERR_X] = psd_routine(dt,Vx, 'fmin',1,'fmax',1e+4,'blocking','lin','binsnumber', 600);
[PSD_Z,F_Z,PSD_ERR_Z] = psd_routine(dt,Vz, 'fmin',1,'fmax',1e+4,'blocking','lin','binsnumber', 600);


% PSD Fitting
for p = 0:1:2
    for q = 0:1:2
        eval(['sx' num2str(p) num2str(q) ' = sum ((F_X.^ (2*' num2str(p) ')) .* (PSD_X.^ (' num2str(q) ')));']);
        eval(['sz' num2str(p) num2str(q) ' = sum ((F_Z.^ (2*' num2str(p) ')) .* (PSD_Z.^ (' num2str(q) ')));']);
    end
end

ax = (sx01 * sx22 - sx11 * sx12) / (sx02 * sx22 - sx12.^2);
bx = (sx11 * sx02 - sx01 * sx12) / (sx02 * sx22 - sx12.^2);
az = (sz01 * sz22 - sz11 * sz12) / (sz02 * sz22 - sz12.^2);
bz = (sz11 * sz02 - sz01 * sz12) / (sz02 * sz22 - sz12.^2);

fcx = sqrt(ax/bx);
fcz = sqrt(az/bz);

kx = 2*pi*gamma*fcx;
kz=  2*pi*gamma*fcz;

D_fit_x = (1/bx) * 2 * (pi.^2);
D_fit_z = (1/bz) * 2 * (pi.^2);

S_fit_x = sqrt(D/D_fit_x);
S_fit_z = sqrt(D/D_fit_z);

psd_fit_x = 1 ./ (ax + bx .* F_X.^2);
psd_fit_z = 1 ./ (az + bz .* F_Z.^2);

% PSD Errors
Times_x = dt*size(Vx,1);
Times_z = dt*size(Vz,1);

xx   = min(F_X) / fcx;
yx   = max(F_X) / fcx;
xz   = min(F_Z) / fcz;
yz   = max(F_Z) / fcz;
sx   = sqrt(pi) * ( (2*yx) / (1 + yx^2) - (2*xx) / (1 + xx^2) + 2 * atan((yx - xx) / (1 + xx*yx)) - ...
    (4/(yx - xx)) * (atan((yx - xx) / (1 + xx*yx)))^2) ^ (-1/2);
sz   = sqrt(pi) * ( (2*yz) / (1 + yz^2) - (2*xz) / (1 + xz^2) + 2 * atan((yz - xz) / (1 + xz*yz)) - ...
    (4/(yz - xz)) * (atan((yz - xz) / (1 + xz*yz)))^2) ^ (-1/2);

fc_err_x = sx * fcx / sqrt(pi * fcx * Times_x);
k_err_x = abs(fc_err_x/fcx)*kx;
fc_err_z = sz * fcz / sqrt(pi * fcz * Times_z);
k_err_z = abs(fc_err_z/fcz)*kz;

gx   = sqrt( ((2*yx)/(1 + yx^2)-(2*xx)/(1 + xx^2) + 2*atan((yx - xx) / (1 + xx*yx)) )/((1 + pi/2)*(yx - xx)) );
gz   = sqrt( ((2*yz)/(1 + yz^2)-(2*xz)/(1 + xz^2) + 2*atan((yz - xz) / (1 + xz*yz)) )/((1 + pi/2)*(yz - xz)) );


D_fit_err_x  = D_fit_x * sqrt( (1 + pi/2) / (pi * fcx * Times_x) )*gx*sx;
D_fit_err_z  = D_fit_z * sqrt( (1 + pi/2) / (pi * fcz * Times_z) )*gz*sz;


%Print results
txt = ['Power spectrum analysis: \n' ...
    int2str(size(Vx,2)) ' signals with ' int2str(size(Vx,1)) ' samples each\n' ...
    int2str(size(Vz,2)) ' signals with ' int2str(size(Vz,1)) ' samples each\n' ...
    '\n' ...
    'fcx : ' num2str(fcx) ' +/- ' num2str(fc_err_x) ' Hz\n' ...
    'fcz : ' num2str(fcz) ' +/- ' num2str(fc_err_z) ' Hz\n' ...
    'kx : ' num2str(kx*1e+6) ' +/- ' num2str(k_err_x*1e+6) ' fN/nm\n' ...
    'kz : ' num2str(kz*1e+6) ' +/- ' num2str(k_err_z*1e+6) ' fN/nm\n' ...
    'D_fit_x : ' num2str(D_fit_x*1e+12) ' +/- ' num2str(D_fit_err_x*1e+12) ' um^2/s\n' ...
    'D_fit_z : ' num2str(D_fit_z*1e+12) ' +/- ' num2str(D_fit_err_z*1e+12) ' um^2/s\n' ...
    '\n'
    ];
fprintf(txt)
%%
hold on
errorbar(F_X,PSD_X*1e+18,PSD_ERR_X*1e+18,'.')
errorbar(F_Z,PSD_Z*1e+18,PSD_ERR_Z*1e+18,'.')
legend('Focal plane (XY)','Beam axis (Z)')

plot(F_X,psd_fit_x*1e+18,'k')
plot(F_Z,psd_fit_z*1e+18,'k')

hold off
xlabel('f [Hz]')
ylabel('PSD [nm^2/Hz]')
set(gca,'XScale','Log','YScale','Log')

