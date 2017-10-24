%% Constants
T = 298;
eta = 0.0008872;
n = 1.33;
tail = [1e-1,1];
kB = 1.38*10^(-23);
lambda = 685*10^(-9); %Laser wavelength
L = 2e-3; %Cuvette thickness
l_star = 401.87e-6; %Mean free path
k0 = 2*pi*n/lambda;
R = 0.5e-6; %Bead radius
D = kB*T/(6*pi*eta*R);
k = 1e-6; %Trap Stiffness
lambda_rel = k/(6*pi*eta*R);
%% TEST ON WATER: MICRORHEOLOGY ROUTINE
t = (1e-4:1e-4:1)';
MSD = 6*D*t;
Jfactor = ((kB*T)/(pi*R))^(-1);
[omega_H20,G_H20]=MSDtoG_Evans_oversampling(t,MSD,1e+4,'Jfactor',Jfactor,'BM',1);
A = [ones(length(omega_H20(1:60)),1),omega_H20(1:60)]; %First 60 pts for fit
Coeff = A\imag(G_H20(1:60));
disp(['ETA = ',num2str(Coeff(2)*1e+3),' mPa/s (cf with eta fixed at beginning)']);
%% TEST ON WATER + TRAP: MICRORHEOLOGY ROUTINE
PACF = exp(-t*lambda_rel);
MSD_T = (6*kB*T/k)*(1 - PACF);
[omega_H20_T,G_H20_T]=MSDtoG_Evans_oversampling(t,MSD_T,1e+4,'Jfactor',Jfactor,'BM',1);
A = [ones(length(omega_H20_T(1:60)),1),omega_H20_T(1:60)]; %First 60 pts for fit of G2
Coeff = A\imag(G_H20_T(1:60));
disp(['ETA = ',num2str(Coeff(2)*1e+3),' mPa/s (cf with eta fixed at beginning)']);
A = [ones(length(omega_H20_T(1:30)),1),omega_H20_T(1:30)]; %First 30 pts for fit of G1
Coeff = A\real(G_H20_T(1:30));
disp(['k = ',num2str(Coeff(1)*(6*pi*R)),' m/N (cf with k fixed at beginning)']);

%% TEST ON WATER: DWS ROUTINE
g1_an = @(x) (L/l_star + 4/3)/(5/3) .* (sinh(x) + (2/3).*x.*cosh(x)) ./ ...
    ((1 + (4/9).*x.^2).*sinh((L/l_star).*x) + (4/3).*x.*cosh((L/l_star).*x));
tDWS = logspace(-6,1,300)';
g2DWS = g1_an(sqrt(k0^2 * 6 * D * tDWS)).^2;
[tau,MSD]=DWS_Analysis(tDWS,g2DWS,'fit','RAT');