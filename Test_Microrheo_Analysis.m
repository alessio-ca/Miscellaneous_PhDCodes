%% Constants
T = 298;
eta = 0.0008872; %Fluid viscosity (in Pa.s)
n = 1.33;
tail = [1e-1,1];
kB = 1.38*10^(-23);
lambda = 685*10^(-9); %Laser wavelength
L = 2e-3; %Cuvette thickness
l_star = 401.87e-6; %Mean free path
k0 = 2*pi*n/lambda;
R = 0.1e-6; %Bead radius
rho_p = 1.04; %Bead density (in g/cm^3)
rho_f = 1.00; %Fluid intensity (in g/cm^3)
m_p = 1e3*rho_p*(4/3)*pi*R^3; %Bead mass (in Kg)
m_f = 1e3*rho_f*(4/3)*pi*R^3; %Displaced fluid mass (in Kg)
D = kB*T/(6*pi*eta*R);
k = 1e-6; %Trap Stiffness
%% TEST ON NEWTONIAN: MICRORHEOLOGY ROUTINE
dt=1e-4;
t = (dt:dt:1e4*dt)';
MSD = 6*D*t;
Jfactor = ((kB*T)/(pi*R))^(-1);
[omega_H20,G_H20]=MSDtoG_Evans_oversampling(t,MSD,1/dt,'Jfactor',Jfactor,'BM',0);
A = [ones(length(omega_H20(1:60)),1),omega_H20(1:60)]; %First 60 pts for fit
Coeff = A\imag(G_H20(1:60));
disp(['ETA = ',num2str(Coeff(2)*1e+3),' mPa/s (cf with eta fixed at beginning)']);
loglog(omega_H20,real(G_H20)+.1,'*-',omega_H20,imag(G_H20),'x-')
legend('G'' + 0.1','G''''','Location','northwest')
xlabel('\omega (rad/s)')
ylabel('Complex Modulus (Pa)')

%% TEST ON SINGLE RELAXATION MAXWELL: MICRORHEOLOGY ROUTINE
dt=1e-5;
t = (dt:dt:1e4*dt)';
Jfactor = ((kB*T)/(pi*R))^(-1);
eta_Max = 20*eta;
E_Max = 1e2;
MSD = (1/Jfactor)*(1/E_Max + t/eta_Max);
[omega_Max,G_Max]=MSDtoG_Evans_oversampling(t,MSD,1/dt,'Jfactor',Jfactor,'BM',0);
Coeff = omega_Max(real(G_Max)>imag(G_Max));
Coeff = Coeff(1); %Crossing point
disp(['TAU_MAX = ',num2str(Coeff),' s (cf with E_Max/eta_Max fixed at beginning)']);
loglog(omega_Max,real(G_Max),'*-',omega_Max,imag(G_Max),'x-')
legend('G''','G''''','Location','northwest')
xlabel('\omega (rad/s)')
ylabel('Complex Modulus (Pa)')

%% TEST ON SINGLE RELAXATION MAXWELL WITH SOLVENT DISSIPATION: MICRORHEOLOGY ROUTINE
dt=5e-7;
t = (dt:dt:1e4*dt)';
Jfactor = ((kB*T)/(pi*R))^(-1);
eta_Max = 10*eta;
G_Max = 2e2;
m_star = m_p + 0.5*m_f;
tau_p = m_star/(6*pi*eta*R);
tau_Max = eta_Max/G_Max;
Omega=sqrt(6*pi*G_Max*R/m_star - 0.25*(1/tau_Max - 1/tau_p)^2);
Damp_const = 0.5*(1/tau_Max + 1/tau_p); 
MSD = (3*kB*T/m_star)*Gen_Maxwell_fun(t,Omega,Damp_const,tau_Max,tau_p);
[omega_Max_Sol,G_Max_Sol,~,~,DeltaA]=MSDtoG_Evans_oversampling(t,MSD,1/dt,'Jfactor',Jfactor,'eps',1.5e-8,'BM',1);
loglog(omega_Max_Sol,real(G_Max_Sol),'*-',omega_Max_Sol,imag(G_Max_Sol),'x-')
legend('G''','G''''','Location','northwest')
xlabel('\omega (rad/s)')
ylabel('Complex Modulus (Pa)')

%% TEST ON WATER + TRAP: MICRORHEOLOGY ROUTINE
t = (1e-5:1e-5:1e-1)';
lambda_rel = k/(6*pi*eta*R);
Gfactor = eta*lambda_rel;
PACF = exp(-t*lambda_rel);
[omega_H20_T,G_H20_T]=PACFtoG_Evans_oversampling(t,PACF,1e+4,'Gfactor',Gfactor);
A = [ones(length(omega_H20_T(1:60)),1),omega_H20_T(1:60)]; %First 60 pts for fit of G2
Coeff = A\imag(G_H20_T(1:60));
disp(['ETA = ',num2str(Coeff(2)*1e+3),' mPa/s (cf with eta fixed at beginning)']);
A = [ones(length(omega_H20_T(1:30)),1),omega_H20_T(1:30)]; %First 30 pts for fit of G1
Coeff = A\real(G_H20_T(1:30));
disp(['k = ',num2str(Coeff(1)*(6*pi*R)),' m/N (cf with k fixed at beginning)']);
loglog(omega_H20_T,real(G_H20_T),'*-',omega_H20_T,imag(G_H20_T),'x-')
legend('G''','G''''','Location','northwest')
xlabel('\omega (rad/s)')
ylabel('Complex Modulus (Pa)')