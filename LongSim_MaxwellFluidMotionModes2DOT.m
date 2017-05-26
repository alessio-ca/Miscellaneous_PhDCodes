%The second order integrator of the Generalized Langevin Equation (GLE) is
%based on Baczweski and Bond (J. Chem. Phys. 129, 044107 (2013))
tic

close all; clc;


%Parameter declaration
dt = 1e-8; %Time step
f_sample=1e8; %Sampling frequency (for simulating real data acquisition)
R = 0.5e-6; %Radius of bead
rho = 2650; %Density of the bead
%rho = 1500; %Density of the bead
%eta = 0.0002; %Viscosity (fictituous)
eta = 0.001;
T = 300; %Temperature
k = [.01 .01 .01]; %Trap elastic constant
%k=[0 0 0];
r_eq = [0 0 0]; %Trap Equilibrium position
N = 1e+8; %Number of steps
run = 10; %Number of trajectories
r0 = [0 0 0]; %Initial position

Nsample=round(N*dt*f_sample);
outfreq=round(1/(dt*f_sample));

%Power law kernel
lambda = 0.8;
gamma_l = 1;
K = @(x) gamma_l/gamma(1-lambda)*x.^(-lambda);

%% Power law kernel, Mode Study
%tau_min=dt/10;
tau_min = dt/100;
%tau_max = N*dt;
tau_max=1e5*tau_min;
t = linspace(tau_min,tau_max,1e6);
figure(1)
plot(log10(t),log10(K(t)),'o','MarkerSize',1)
hold on
y = K(t)';
t = t';
modes = [2 4 6 8];
for i=1:length(modes)
    tau_k = logspace(log10(tau_min),log10(tau_max),modes(i));
    X = exp(-t./tau_k);
    a = X\y;
    plot(log10(t),log10(X*a))
end
hold off
legend('Exact','2 Modes','4 Modes','6 Modes','8 Modes')
xlabel('Time (Log)')
ylabel('K(t) (Log)')
drawnow
%%

% Pre-calculation coefficients
m = (4./3.)*pi*R.^3*rho; 
kB=1.38e-23; %Boltzmann constant
kBT=kB*T;
D = kBT/(6*pi*eta*R);
gamma_kern = kBT/D;

nmodes=2;
gamma_l=gamma_kern;
clear gamma
K = @(x) gamma_l/gamma(1-lambda)*x.^(-lambda);
y = K(t);
tau_rel = logspace(log10(tau_min),log10(tau_max),nmodes);
tau_rel(1)=1e-15;
X = exp(-t./tau_rel);
gamma_kern = (X\y).*tau_rel';
gamma_kern(1)= 10*gamma_kern(end);
gamma_kern(end)=gamma_kern(end)*500;
theta = exp(-dt./tau_rel)';
alpha = (1-theta)./sqrt(dt)';
%%

% Inizialization

r = repmat(r0,1,run);
r_sample=zeros(Nsample,3*run); 
v_sample=zeros(Nsample,3*run); 
v = zeros(1,3*run);
S = zeros(nmodes,3*run);

%%Equilibrate%%
for n = 1:100
    vhalf = v + dt/(2*m)*(-repmat(k,1,run).*(r-repmat(r_eq,1,run)) + sum(S));
    r = r + dt*vhalf;
    S = theta.*S - (1 - theta).*gamma_kern.*repmat(vhalf,[nmodes,1]) + alpha.*sqrt(2*kBT*gamma_kern).*randn(nmodes,3*run);
    Favg = mean(sum(S)-repmat(k,1,run).*(r-repmat(r_eq,1,run))); %Center of mass acceleration, to be substracted
    v = vhalf + dt/(2*m)*(-repmat(k,1,run).*(r-repmat(r_eq,1,run)) + sum(S) - Favg);
end

%%Sample%%
tic
OnTheFly_CoarseGrainACF(r,dt,30,0);
for n = 1:N
    vhalf = v + dt/(2*m)*(-repmat(k,1,run).*(r-repmat(r_eq,1,run)) + sum(S));
    r = r + dt*vhalf;
    S = theta.*S - (1 - theta).*gamma_kern.*repmat(vhalf,[nmodes,1]) + alpha.*sqrt(2*kBT*gamma_kern).*randn(nmodes,3*run);
    Favg = mean(sum(S)-repmat(k,1,run).*(r-repmat(r_eq,1,run))); %Center of mass acceleration, to be substracted
    v = vhalf + dt/(2*m)*(-repmat(k,1,run).*(r-repmat(r_eq,1,run)) + sum(S) - Favg);
    if (mod(n,outfreq)==0)
        r_sample(n/outfreq,:)=r;
        v_sample(n/outfreq,:)=v;
    end
    if (mod(100*n,N)==0)
        disp([num2str(100*n/N),' %'])
    end
    OnTheFly_CoarseGrainACF(r,dt,30,1);
end
[pacfCG,t_vacfCG,pacfCG_err]=OnTheFly_CoarseGrainACF(r,dt,30,2,'taumax',N*dt);