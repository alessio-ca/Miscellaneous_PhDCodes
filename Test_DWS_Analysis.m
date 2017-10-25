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
rho = 1.04; %Bead density (in g/cm^3)
m = 1e3*rho*(4/3)*pi*R^3; %Bead mass (in Kg)
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
[tau_CON,MSD_CON]=DWS_Analysis(tDWS,g2DWS,'fit','CON');
[tau_RAT,MSD_RAT]=DWS_Analysis(tDWS,g2DWS,'fit','RAT');
[tau_SPL,MSD_SPL]=DWS_Analysis(tDWS,g2DWS,'fit','SPL');
%%
loglog(tau_CON,MSD_CON,tau_RAT,MSD_RAT,tau_SPL,MSD_SPL,tau_CON,6*D*tau_CON,'k--')

%% TEST ON MAXWELL: DWS ROUTINE
g1_an = @(x) (L/l_star + 4/3)/(5/3) .* (sinh(x) + (2/3).*x.*cosh(x)) ./ ...
    ((1 + (4/9).*x.^2).*sinh((L/l_star).*x) + (4/3).*x.*cosh((L/l_star).*x));
tDWS = logspace(-10,1,300)';
tau_Max = eta/1e2;
tau_Br = m/(6*pi*eta*R);
ACF_DWS = (3*kB*T/m)*exp(-tDWS/(2*tau_Max)) .* ...
    (cos(sqrt(4*(tau_Max/tau_Br) - 1)/(2*tau_Max)*tDWS) ...
    + 1/sqrt(4*(tau_Max/tau_Br) - 1) * sin(sqrt(4*(tau_Max/tau_Br) - 1)/(2*tau_Max)*tDWS) ...
    );
MSD_DWS = (6*kB*T/m)*(...
    tau_Br*tDWS + tau_Br*(tau_Br - tau_Max) ...
    * (exp(-tDWS/(2*tau_Max)) .* cos(sqrt(4*(tau_Max/tau_Br) - 1)/(2*tau_Max)*tDWS) - 1) ...
    + 1/4 * tau_Br^2 * exp(-tDWS/(2*tau_Max)) .* sin(sqrt(4*(tau_Max/tau_Br) - 1)/(2*tau_Max)*tDWS) ...
    * (1/sqrt(4*(tau_Max/tau_Br) - 1) - 3*sqrt(4*(tau_Max/tau_Br) - 1)) ...
    );
%g2DWS = g1_an(sqrt(k0^2 * (kB*T)/(pi*R) * JDWS)).^2;
%[tau_CON,MSDMAX_CON]=DWS_Analysis(tDWS,g2DWS,'fit','CON');
%[tau_RAT,MSDMAX_RAT]=DWS_Analysis(tDWS,g2DWS,'fit','RAT');
%[tau_SPL,MSDMAX_SPL]=DWS_Analysis(tDWS,g2DWS,'fit','SPL');
%%
loglog(tau_CON,MSDMAX_CON,tau_RAT,MSDMAX_RAT,tau_SPL,MSDMAX_SPL,tau_CON,(kB*T)/(pi*R) * (tau_CON/eta + 1/1e2),'k--')
