%The second order integrator of the Generalized Langevin Equation (GLE) is
%based on Baczweski and Bond (J. Chem. Phys. 129, 044107 (2013))
tic

close all; clc;

%Parameter declaration

dt = 1e-8; %Time step
f_sample=1e7; %Sampling frequency (for simulating real data acquisition)
R = 0.5e-6; %Radius of bead
rho = 2650; %Density of the bead
%rho = 1500; %Density of the bead
%eta = 0.0002; %Viscosity (fictituous)
eta = 0.001;
tau_rel = 1e-15; %Relaxation time
T = 300; %Temperature
k = [1000e-6 1000e-6 1000e-6]; %Trap elastic constant
%k=[0 0 0];
r_eq = [0 0 0]; %Trap Equilibrium position
N = 1e+7; %Number of steps
run = 10; %Number of trajectories
r0 = [0 0 0]; %Initial position

Nsample=round(N*dt*f_sample);
outfreq=round(1/(dt*f_sample));



% Pre-calculation coefficients
m = (4./3.)*pi*R.^3*rho; 
kB=1.38e-23; %Boltzmann constant
kBT=kB*T;
D = kBT/(6*pi*eta*R);
gamma = kBT/D;
theta = exp(-dt/tau_rel);
alpha = (1-theta)/sqrt(dt);
disp(['Langevin Time Scale: ',num2str(m/(6*pi*eta*R))])
%%
% Inizialization

r = repmat(r0,1,run);
r_sample=zeros(Nsample,3*run); 
v_sample=zeros(Nsample,3*run); 
v = zeros(1,3*run);
S = zeros(1,3*run);
Omega = sqrt(gamma/(m*tau_rel) - 1/(4*tau_rel.^2));
disp(['Omega = ', num2str(Omega)])
disp(['Tau_p = ', num2str(m/gamma)])
disp(['Tau_M = ', num2str(1/(tau_rel.^(-1) + gamma/m))])

%%Equilibrate%%
for n = 1:100
    vhalf = v + dt/(2*m)*(-repmat(k,1,run).*(r-repmat(r_eq,1,run)) + S);
    r = r + dt*vhalf;
    S = theta*S - (1 - theta)*gamma*vhalf + alpha*sqrt(2*kBT*gamma)*randn(1,3*run);
    Favg = mean(S-repmat(k,1,run).*(r-repmat(r_eq,1,run))); %Center of mass acceleration, to be substracted
    v = vhalf + dt/(2*m)*(-repmat(k,1,run).*(r-repmat(r_eq,1,run)) + S - Favg);
end

%%Sample%%
tic
%OnTheFly_CoarseGrainACF(v(1,:),dt,20,0);
for n = 1:N
    vhalf = v + dt/(2*m)*(-repmat(k,1,run).*(r-repmat(r_eq,1,run)) + S);
    r = r + dt*vhalf;
    S = theta*S - (1 - theta)*gamma*vhalf + alpha*sqrt(2*kBT*gamma)*randn(1,3*run);
    Favg = mean(S-repmat(k,1,run).*(r-repmat(r_eq,1,run))); %Center of mass acceleration, to be substracted
    v = vhalf + dt/(2*m)*(-repmat(k,1,run).*(r-repmat(r_eq,1,run)) + S - Favg);
    if (mod(n,outfreq)==0)
        r_sample(n/outfreq,:)=r;
        v_sample(n/outfreq,:)=v;
    end
    if (mod(100*n,N)==0)
        disp([num2str(100*n/N),' %'])
    end
    %OnTheFly_CoarseGrainACF(v(n,:),dt,20,1);
end
