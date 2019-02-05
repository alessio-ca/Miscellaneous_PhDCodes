Fs=1000;
T=1/Fs;

%%
%FFT section
tau=(T:T:20*pi);
L=2*floor(length(tau)/2);
df=Fs/L;
lambda=20;
ustep=@(x) x>=0;
autocorr=exp(-lambda*tau);
autocorrFFT=fft(autocorr,L)*T;
autocorrFFT=autocorrFFT(1:floor(L/2));


f=df*(-L/2 : L/2);
f=f(L/2+2:end);
omegaFFT=2*pi*f;
figure(1)
semilogx(omegaFFT,real(autocorrFFT))
hold on
semilogx(omegaFFT,imag(autocorrFFT))
hold off
%%
%FFT of Pi section
autocorrP=1./(1i.*omegaFFT) - autocorrFFT;
figure(1)
semilogx(omegaFFT,real(autocorrP))
hold on
semilogx(omegaFFT,imag(autocorrP))
hold off

figure(2)
loglog(omegaFFT,abs(autocorrFFT))
hold on
loglog(omegaFFT,abs(autocorrP))
hold off
%%

%G modulus
G=1./(1 - 1i.*omegaFFT.*autocorrFFT);
Ga=1./(1i.*omegaFFT.*autocorrP);

figure(3)
loglog(omegaFFT,real(G),omegaFFT,real(Ga))
hold on
loglog(omegaFFT,imag(G),omegaFFT,imag(Ga))
hold off

%%
%Evans section (no oversampling)
tau=(T:T:2*pi);
L=length(tau); 
omega=2*pi*(Fs/L)*linspace(1,floor(L/2));
autocorr=exp(-lambda*tau);
Evansfft=1i*omega;
Evansfft=Evansfft + (1-exp(-1i*omega*tau(1)))*(autocorr(1)-1)/tau(1);
%Enable progress bar for parallel pool
try
    parpool;
catch ME
    if ~strcmp(ME.identifier,'parallel:convenience:ConnectionOpen')
        rethrow(ME)
    end
end
warning('off','MATLAB:Java:DuplicateClass')
pctRunOnAll javaaddpath java
progressStep=ceil(length(omega)/100);
ppm = ParforProgMon('Progress: ', length(omega), progressStep, 300, 80);

parfor ii=1:length(omega)
    Evansfft(ii)=Evansfft(ii)-sum(diff(exp(-1i*omega(ii)*tau)).*diff(autocorr)./diff(tau));
    if mod(ii,progressStep)==0
        ppm.increment();
    end
end
ppm.delete()

Evansfft=Evansfft./(-omega.^2);

figure(1)
semilogx(omega,real(Evansfft))
hold on
semilogx(omega,imag(Evansfft))
hold off

Gev=1./(1 - 1i*omega.*Evansfft);

figure(2)
loglog(omega/lambda,real(Gev))
hold on
loglog(omega/lambda,imag(Gev))
hold off

%%
%Evans section (oversampling)
tau=(T:T:4*pi);
L=length(tau); 
omega=2*pi*(Fs/L)*linspace(1,floor(L/2));
autocorr=exp(-lambda*tau);
acf_spline=csape(tau,autocorr);

Beta=20;
oversample_tau=(1/(Beta*Fs):1/(Beta*Fs):tau(end));
oversample_acf=fnval(oversample_tau,acf_spline);

figure(1)
semilogx(oversample_tau,oversample_acf)
hold on
semilogx(tau,autocorr)
hold off


Evansfft=1i*omega;
Evansfft=Evansfft + (1-exp(-1i*omega*oversample_tau(1)))*(oversample_acf(1)-1)/oversample_tau(1);
%Enable progress bar for parallel pool
try
    parpool;
catch ME
    if ~strcmp(ME.identifier,'parallel:convenience:ConnectionOpen')
        rethrow(ME)
    end
end
warning('off','MATLAB:Java:DuplicateClass')
pctRunOnAll javaaddpath java
progressStep=ceil(length(omega)/100);
ppm = ParforProgMon('Progress: ', length(omega), progressStep, 300, 80);

parfor ii=1:length(omega)
    Evansfft(ii)=Evansfft(ii)-sum(diff(exp(-1i*omega(ii)*oversample_tau)).*diff(oversample_acf)./diff(oversample_tau));
    if mod(ii,progressStep)==0
        ppm.increment();
    end
end
ppm.delete()

Evansfft=Evansfft./(-omega.^2);
Gev=1./(1 - 1i*omega.*Evansfft);

figure(2)
loglog(omega/lambda,real(Gev))
hold on
loglog(omega/lambda,imag(Gev))
hold off


