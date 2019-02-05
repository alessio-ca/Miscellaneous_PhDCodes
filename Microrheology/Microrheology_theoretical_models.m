R = 0.5e-6; %Radius of bead
rho = 2650; %Density of the bead
m = (4./3.)*pi*R.^3*rho; 
eta = 0.001;
k = 1000e-6;
T = 300; %Temperature
kB=1.38e-23; %Boltzmann constant
kBT=kB*T;

tau_rel = 1e-15; %Relaxation time

tau_rel=m/(6*pi*eta*R);
Gfactor=k/(6*pi*R);

Omega=sqrt(k(1)/m - 0.25/tau_rel^2);
CosImag = @(x) (exp(-x*(imag(Omega) + 1/(2*tau_rel))) + exp(x*(imag(Omega) - 1/(2*tau_rel))))/2;
SinImag = @(x) (-exp(-x*(imag(Omega) + 1/(2*tau_rel))) + exp(x*(imag(Omega) - 1/(2*tau_rel))))/2;

dt=1e-12;
tau=(0:dt:1e6*dt)';
%pacf_th=(kBT/m)*(CosImag(tau) + 1/(2*tau_rel*imag(Omega))*SinImag(tau));
vacf_th=(kBT/m)*(CosImag(tau) - 1/(2*tau_rel*imag(Omega))*SinImag(tau));
[omega_VACF_th,GVACF_th]=VACFtoG_Evans(tau,vacf_th,1/dt);
%%
loglog(omega_VACF_th,real(GVACF_th),'-x','Color',[214 180 13] ./ 255)
hold on
loglog(omega_VACF_th,imag(GVACF_th),'-o','Color',[214 180 13] ./ 255)

loglog(omega_VACF_th,Gfactor*ones(length(omega_VACF_th),1),'k-.',omega_VACF_th,eta*omega_VACF_th,'k-')
hold off
