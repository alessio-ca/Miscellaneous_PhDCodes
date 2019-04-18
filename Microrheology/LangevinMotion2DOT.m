tic
%The second order integrator of the Langevin Equation is based on
%Vanden-Eijden and Ciccotti (Chem. Phys. Lett. 429, 310316 (2006))

close all; clc;

%Parameter declaration

dt = 1e-8; %Time step
dt0_5 = sqrt(dt);
dt1_5 = dt.^(1.5);
dt2 = dt.^2;
R = 1e-6; %Radius of bead
rho = 2650; %Density of the bead
eta = 0.001; %Viscosity
T = 300; %Temperature
k = [1e-6 1e-6]; %Trap elastic constant
r_eq = [0 0]; %Trap Equilibrium position
N = 5e5; %Number of steps

r0 = [0 0]; %Initial position

% Pre-calculation coefficients
m = (4./3.)*pi*R.^3*rho; 
kB=1.38e-23; %Boltzmann constant
kBT=kB*T;
D = kBT/(6*pi*eta*R);
gamma = kBT/D;
gamma = gamma/m; %scaled friction coefficient
k = k/m; %scaled elastic constant
sigma = sqrt(2*kBT*gamma/m);
C1 = k.*r_eq*(0.5*dt.^2);
eps = randn(N,2); % noise, epsilon
th = randn(N,2); % noise, theta

% Inizialization
r = zeros(N,2);
v = zeros(N,2);
r(1,:) = r0;
v(1,:) = [0 0];
t = (0:1:N-1)'*dt;

% Main loop
% It simulates a Brownian Particle in a solvent according to the
% non-overdamped Langevin equation. 
for n = 2:1:N
    C = (0.5*dt2)*(-k.*r(n-1,:)- gamma*v(n-1,:)) + C1 + sigma*dt1_5*(0.5*eps(n,:) + (0.5/sqrt(3.))*th(n,:));
    r(n,:) = r(n-1,:) + dt*v(n-1,:) + C;
    v(n,:) = v(n-1,:) + 0.5*dt*(-k.*r(n,:) - k.*r(n-1,:)) + 2*C1/dt - dt*gamma*v(n-1,:) + sigma*dt0_5*eps(n,:) - gamma*C;
    if mod(n,100000)==0
        disp(n)
    end
end
    

[msd,t_msd]=msd_routine(dt,r,'TauMax',dt*length(t)/10);
[vacf,t_vacf]=acf_routine(dt,v,'TauMax',dt*length(t)/100);

%%
figure(1)
loglog(t_msd,msd)
hold on
loglog(t_msd,2*D*t_msd)
hold off
title('Mean Square Displacement')
xlabel('Time (s)')
ylabel('MSD (m^2/s)')

figure(2)
semilogx(t_vacf,vacf)
hold on
semilogx(t_vacf,exp(-gamma*t_vacf))
hold off
title('Velocity Autocorrelation Function')
xlabel('Time (s)')
ylabel('VACF')

toc