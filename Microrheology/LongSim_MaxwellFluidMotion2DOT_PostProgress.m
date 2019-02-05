%f_sample=1e5;
f_sample=1e7;
R = 0.5e-6; %Radius of bead
rho = 2650; %Density of the bead
m = (4./3.)*pi*R.^3*rho; 
eta = 0.001;
k = 1e-3;
T = 300; %Temperature
kB=1.38e-23; %Boltzmann constant
kBT=kB*T;

tau_rel = 1e-15; %Relaxation time



%%

[pacf,tau]=acf_routine(1/f_sample,r_sample,'taumax',1e-3);
[~,GPACF]=PACFtoG_Evans(tau,pacf,f_sample);
[vacf,tau]=acf_routine(1/f_sample,v_sample,'taumax',1e-3,'norm',0);
[~,GVACF]=VACFtoG_Evans(tau,vacf,f_sample);
[vacf_r,tau]=acf_routine(1/f_sample,diff(r_sample)*f_sample,'taumax',1e-3,'norm',0);
[omega,GVACF_R]=VACFtoG_Evans(tau,vacf_r,f_sample);

%%
tau_rel=m/(6*pi*eta*R);
Omega=sqrt(k(1)/m - 0.25/tau_rel^2);
CosImag = @(x) (exp(-x*(imag(Omega) + 1/(2*tau_rel))) + exp(x*(imag(Omega) - 1/(2*tau_rel))))/2;
SinImag = @(x) (-exp(-x*(imag(Omega) + 1/(2*tau_rel))) + exp(x*(imag(Omega) - 1/(2*tau_rel))))/2;


figure(1)
semilogx(tau,pacf)
hold on
%semilogx(tau,(CosImag(tau) + 1/(2*tau_rel*imag(Omega))*SinImag(tau)))
%hold off
%figure(2)
semilogx(tau,vacf./vacf(1,:))

hold on
%semilogx(tau,(CosImag(tau) - 1/(2*tau_rel*imag(Omega))*SinImag(tau)))
%semilogx(tau,exp(-tau./tau_rel))
hold off
%%
pacf_th=(CosImag(tau) + 1/(2*tau_rel*imag(Omega))*SinImag(tau));
[~,GPACF_th]=PACFtoG_Evans(tau,pacf_th,f_sample);
vacf_th=(kBT/m)*(CosImag(tau) - 1/(2*tau_rel*imag(Omega))*SinImag(tau));
%vacf_th=GLE_VACF(tau(1:1:500),sqrt(100e-6/m),0.8,6.5e+6);
[omega_VACF_th,GVACF_th]=VACFtoG_Evans(tau,vacf_th,f_sample);
%[omega_VACF_th,GVACF_th]=VACFtoG_Evans(tau(1:1:500),vacf(1)*vacf_th,f_sample/10);
%%
figure(3)
Gfactor=k/(6*pi*R);
loglog(omega,real(GVACF),'x','Color',[142 42 111] ./ 255)
hold on
loglog(omega,imag(GVACF),'o','Color',[142 42 111] ./ 255)

loglog(omega,real(GVACF_R),'x','Color',[246 127 41] ./ 255)
loglog(omega,imag(GVACF_R),'o','Color',[246 127 41] ./ 255)

loglog(omega_VACF_th,real(GVACF_th),'-x','Color',[214 180 13] ./ 255)
loglog(omega_VACF_th,imag(GVACF_th),'-o','Color',[214 180 13] ./ 255)

loglog(omega,Gfactor*ones(length(omega),1),'k-.',omega,eta*omega,'k-')
hold off
%%
figure(4)
loglog(omega,Gfactor*real(GPACF),'bx',omega,Gfactor*imag(GPACF),'bo')
hold on
loglog(omega,Gfactor*real(GPACF_th),'-x','Color',[214 180 13] ./ 255)
loglog(omega,Gfactor*imag(GPACF_th),'-o','Color',[214 180 13] ./ 255)
loglog(omega,Gfactor*ones(length(omega),1),'k-.',omega,eta*omega,'k-')
hold off
